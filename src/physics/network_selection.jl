# Selecting which reactions belong in a network from seed species.

# Is `name` a plausible H-Ca nova network species: a light universal particle
# (n/p/d/t/he3/he4), or an isotope within a bounded band around the valley of
# stability? Bounds are `Z <= zmax`, `A <= amax`, and neutron number
# `N = A - Z` within `[Z - proton_rich_margin, Z + neutron_rich_margin]` of
# `Z` -- i.e. isotopes up to a few nucleons off stability in either
# direction, but not the far proton/neutron drip-line species that would
# never be populated by hydrogen burning at nova temperatures. This keeps
# `select_h_ca_reaction_tables`'s BFS closure (below) from wandering into
# thousands of physically irrelevant exotic-nuclide reactions that a
# library-wide rate compilation (REACLIB) inevitably contains.
function _h_ca_candidate_species(
    name::AbstractString;
    zmax::Int=20,
    amax::Int=46,
    proton_rich_margin::Int=4,
    neutron_rich_margin::Int=3,
)
    normalized = normalize_species_name(name)
    normalized in ("n", "p", "d", "t", "he3", "he4") && return true

    species = try
        species_from_name(normalized)
    catch
        return false
    end

    species.Z < 1 && return false
    species.Z > zmax && return false
    species.A > amax && return false

    neutron_number = species.A - species.Z
    neutron_number < max(0, species.Z - proton_rich_margin) && return false
    neutron_number > species.Z + neutron_rich_margin && return false
    return true
end

"""
    select_h_ca_reaction_tables(tables, seed_species)

Select a broader H-Ca nova post-processing network from a rate-table library
(STARLIB or REACLIB). The selector is deliberately less broad than "all H-Ca
rows": it starts from the supplied seed composition, follows weak, proton,
alpha (and optionally neutron/`he3`/`d`) links, and keeps species inside a
bounded stable/proton-rich isotope band (`_h_ca_candidate_species`). This
keeps single-zone runs tractable while covering the main H-Ca nova flows,
rather than pulling in the library's full complement of far-from-stability
reactions that a nova envelope never reaches.

Algorithm, run in two passes:
1. **Candidacy filter**: keep a table only if its bookkeeping is valid
   (`_valid_nuclear_bookkeeping`), every participant is inside the physical
   band, and it is either a recognized single-reactant weak decay
   (`_weak_source`, when `include_weak=true`) or has at least one of
   `projectiles` among its reactants.
2. **Reachability closure**: starting from `seed_species` (plus `n`, `p`,
   `he4`, always available), repeatedly scan the candidates for any whose
   *every* reactant is already reachable; select it and add its products to
   the reachable set. Repeat to a fixed point. This is what makes the
   network "broader than the seed composition" -- e.g. seeding only ¹²C lets
   the closure pull in ¹³N, ¹³C, ¹⁴N, ... as the proton-capture chain
   unfolds -- while still excluding anything whose precursor is never
   actually produced.

A reaction can therefore be excluded for two structurally different reasons:
it never became a *candidate* at all (outside the physical band, wrong
bookkeeping, no seed projectile), or it was a valid candidate but is
*unreachable* (its reactants never all appeared in the closure). See
`SensitivityStudy/network_audit.jl` for a diagnostic that reports both
reasons explicitly for a given seed composition.
"""
function select_h_ca_reaction_tables(
    tables::AbstractVector{ReactionRateTable},
    seed_species;
    zmax::Int=20,
    amax::Int=46,
    proton_rich_margin::Int=4,
    neutron_rich_margin::Int=3,
    include_weak::Bool=true,
    projectiles=("p", "he4", "he3", "d"),
)
    normalized_seed_species = Set(normalize_species_name.(collect(seed_species)))
    push!(normalized_seed_species, "n")
    push!(normalized_seed_species, "p")
    push!(normalized_seed_species, "he4")

    physical_species(name) = _h_ca_candidate_species(
        name;
        zmax=zmax,
        amax=amax,
        proton_rich_margin=proton_rich_margin,
        neutron_rich_margin=neutron_rich_margin,
    )

    function candidate_table(table)
        _valid_nuclear_bookkeeping(table) || return false
        all(physical_species, vcat(table.reactants, table.products)) || return false
        include_weak && length(table.reactants) == 1 && _weak_source(table.source) && return true
        return any(projectile -> projectile in table.reactants, projectiles)
    end

    candidates = [table for table in tables if candidate_table(table)]
    selected = ReactionRateTable[]
    selected_keys = Set{Any}()
    species = copy(normalized_seed_species)

    changed = true
    while changed
        changed = false
        for table in candidates
            all(reactant -> reactant in species, table.reactants) || continue

            key = _reaction_participant_key(table)
            if !(key in selected_keys)
                push!(selected, table)
                push!(selected_keys, key)
            end

            for product in table.products
                if !(product in species)
                    push!(species, product)
                    changed = true
                end
            end
        end
    end

    return selected
end

# Heuristic recognition of a weak-decay-sourced rate table by its source
# label: STARLIB/REACLIB weak-rate labels conventionally end in "w" (e.g.
# "wc07w"), or spell out "beta"/"ec" (electron capture)/"weak" directly.
# There's no structural field marking a table as weak, so source-label
# pattern matching is the practical way to distinguish "this decays" from
# "this needs a projectile."
function _weak_source(source::AbstractString)
    s = lowercase(strip(source))
    return occursin("w", s) ||
           s == "ec" ||
           occursin("bet", s) ||
           occursin("weak", s)
end

"""
    select_decay_reaction_tables(tables, seed_species)

Select the weak-decay-only sub-network reachable from `seed_species`, for
post-processing decay (`decay_mass_fractions`). Uses the same reachability
closure as `select_h_ca_reaction_tables`, but restricted to single-reactant
weak-decay tables (`_weak_source`) with no physical-band or projectile
filtering -- once burning has stopped and only decay chains matter, there is
no need to bound the isotope band, since decay only ever moves toward
stability (never away from it).
"""
function select_decay_reaction_tables(tables::AbstractVector{ReactionRateTable}, seed_species)
    species = Set(normalize_species_name.(collect(seed_species)))
    candidates = [
        table for table in tables
        if length(table.reactants) == 1 &&
           _valid_nuclear_bookkeeping(table) &&
           _weak_source(table.source)
    ]

    selected = ReactionRateTable[]
    selected_keys = Set{Any}()
    changed = true
    while changed
        changed = false
        for table in candidates
            only(table.reactants) in species || continue

            key = _reaction_participant_key(table)
            if !(key in selected_keys)
                push!(selected, table)
                push!(selected_keys, key)
            end

            for product in table.products
                if !(product in species)
                    push!(species, product)
                    changed = true
                end
            end
        end
    end

    return selected
end

