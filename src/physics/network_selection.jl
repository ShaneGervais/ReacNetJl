# Selecting which reactions belong in a network from seed species.

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

Select a broader H-Ca nova post-processing network from STARLIB. The selector is
deliberately less broad than "all H-Ca rows": it starts from the supplied seed
composition, follows weak, proton, and alpha links, and keeps species inside a
bounded stable/proton-rich isotope band. This keeps single-zone examples
tractable while covering the main H-Ca nova flows.
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

function _weak_source(source::AbstractString)
    s = lowercase(strip(source))
    return occursin("w", s) ||
           s == "ec" ||
           occursin("bet", s) ||
           occursin("weak", s)
end

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

