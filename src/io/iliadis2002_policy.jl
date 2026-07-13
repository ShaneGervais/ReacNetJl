# Iliadis 2002 baseline rate policy: merges REACLIB, label sets, and paper tables.

"""
    iliadis2002_rate_tables(sets; T9_grid=STARLIB_T9_GRID, boundary_A=20,
                            include_weak=true, include_reverse=false)

Build the baseline rate tables of the Iliadis et al. (2002, ApJS 142, 105)
nova sensitivity study from a REACLIB library:

- NACRE (`nacr`, Angulo et al. 1999) for reactions whose heaviest reactant has
  A < `boundary_A`
- Iliadis et al. (2001, ApJS 134, 151; `il01`) for A >= `boundary_A`, with
  `nacr` as fallback

Weak rates, which are outside both compilations, are kept from whatever set
label the library provides when `include_weak=true`. Reactions covered by
neither compilation fall back to their available label and are reported in
`report.fallbacks` so rate provenance stays explicit.

Returns `(tables=..., report=...)` where `report` carries per-reaction chosen
sources (`chosen`), counts by category (`counts`), and the fallback list.
Use `read_reaclib` on the ReaclibV1.0 snapshot (`data/reaclib_v1.0.dat`),
which still contains the complete `nacr` and `il01` sets.
"""
function iliadis2002_rate_tables(
    sets::AbstractVector{ReaclibSet};
    T9_grid::AbstractVector{<:Real}=STARLIB_T9_GRID,
    boundary_A::Int=20,
    include_weak::Bool=true,
    include_reverse::Bool=false,
    paper_tables::AbstractVector{ReactionRateTable}=ReactionRateTable[],
)
    deduped_sets = _dedup_reaclib_sets(sets)
    reaction_order = Any[]
    by_reaction = Dict{Any,Dict{String,Vector{ReaclibSet}}}()
    for set in deduped_sets
        set.reverse && continue
        key = (set.chapter, Tuple(set.reactants), Tuple(set.products))
        if !haskey(by_reaction, key)
            by_reaction[key] = Dict{String,Vector{ReaclibSet}}()
            push!(reaction_order, key)
        end
        label_map = by_reaction[key]
        label = lowercase(set.label)
        haskey(label_map, label) || (label_map[label] = ReaclibSet[])
        push!(label_map[label], set)
    end

    tables = ReactionRateTable[]
    chosen = NamedTuple[]
    counts = Dict{Symbol,Int}()
    chosen_forward_labels = Dict{Any,String}()

    for key in reaction_order
        chapter, reactants, products = key
        label_map = by_reaction[key]
        available = sort!(collect(keys(label_map)))

        target_A = 0
        for name in reactants
            species = try
                species_from_name(name)
            catch
                continue
            end
            target_A = max(target_A, species.A)
        end

        preferred = target_A >= boundary_A ? ("il01", "nacr") : ("nacr", "il01")
        label = nothing
        category = :other
        for wanted in preferred
            if wanted in available
                label = wanted
                category = wanted == "il01" ? :il01 : :nacr
                break
            end
        end

        if label === nothing
            weak_labels = [l for l in available if any(set -> set.resonance == 'w', label_map[l]) || _weak_source(l)]
            if !isempty(weak_labels)
                label = first(weak_labels)
                category = :weak
            else
                label = first(available)
                category = :other
            end
        end

        category == :weak && !include_weak && continue

        table = _reaclib_group_table(label_map[label], T9_grid)
        push!(tables, table)
        counts[category] = get(counts, category, 0) + 1
        chosen_forward_labels[(Tuple(reactants), Tuple(products))] = label
        push!(chosen, (
            reaction=_reaction_display(table.reactants, table.products),
            chapter=chapter,
            target_A=target_A,
            preferred=first(preferred),
            label=label,
            source=table.source,
            category=category,
            available=available,
        ))
    end

    if include_reverse
        keep(set) = set.reverse
        order, groups = _group_reaclib_sets(deduped_sets, keep)
        for key in order
            group = groups[key]
            representative = first(group)
            forward_key = (Tuple(representative.products), Tuple(representative.reactants))
            get(chosen_forward_labels, forward_key, nothing) == lowercase(representative.label) || continue

            table = _reaclib_group_table(group, T9_grid)
            push!(tables, table)
            counts[:reverse] = get(counts, :reverse, 0) + 1
            push!(chosen, (
                reaction=_reaction_display(table.reactants, table.products),
                chapter=representative.chapter,
                target_A=0,
                preferred=lowercase(representative.label),
                label=lowercase(representative.label),
                source=table.source,
                category=:reverse,
                available=[lowercase(representative.label)],
            ))
        end
    end

    # Paper tables (tabulated rates from the original publications) override
    # the REACLIB fits of the same reactions; reactions absent from REACLIB
    # are appended. This closes the compilation-coverage gap of the JINA
    # database with the authoritative published numbers.
    replaced = 0
    added = 0
    if !isempty(paper_tables)
        index_by_participants = Dict{Any,Int}(_reaction_participant_key(t) => i for (i, t) in pairs(tables))
        paper_slots = Set{Int}()
        for paper_table in paper_tables
            key = _reaction_participant_key(paper_table)
            slot = get(index_by_participants, key, nothing)
            if slot === nothing
                push!(tables, paper_table)
                slot = length(tables)
                index_by_participants[key] = slot
                push!(paper_slots, slot)
                added += 1
            elseif slot in paper_slots
                # first paper table wins: il01 keeps proton-induced A>=20
                # reactions that NACRE also tabulates
                continue
            else
                tables[slot] = paper_table
                push!(paper_slots, slot)
                replaced += 1
            end

            target_A = 0
            for name in paper_table.reactants
                species = try
                    species_from_name(name)
                catch
                    continue
                end
                target_A = max(target_A, species.A)
            end
            push!(chosen, (
                reaction=_reaction_display(paper_table.reactants, paper_table.products),
                chapter=paper_table.chapter,
                target_A=target_A,
                preferred=paper_table.source,
                label=paper_table.source,
                source=paper_table.source,
                category=:paper,
                available=[paper_table.source],
            ))
        end
        counts[:paper] = replaced + added
    end

    fallbacks = [entry for entry in chosen if entry.category == :other]
    report = (counts=counts, chosen=chosen, fallbacks=fallbacks, paper_overrides=(replaced=replaced, added=added))
    return (tables=tables, report=report)
end

function iliadis2002_rate_tables(path::AbstractString; kwargs...)
    return iliadis2002_rate_tables(read_reaclib(path); kwargs...)
end

"""
The zero-argument form reads the ReaclibV1.0 snapshot for library-wide
coverage (weak rates, neutron captures, and other channels outside the two
compilations), merges in the complete label-specific `nacr` and `il01`
downloads, and overrides with the tabulated paper rates
(`data/iliadis2001_rates.dat`, `data/nacre_rates.dat`) when those files are
present in `data/`.
"""
function iliadis2002_rate_tables(; kwargs...)
    sets = read_reaclib(_default_reaclib_path())
    for path in (DEFAULT_REACLIB_IL01_PATH, DEFAULT_REACLIB_NACR_PATH)
        isfile(path) && append!(sets, read_reaclib(path))
    end

    if !haskey(kwargs, :paper_tables)
        pf = isfile(DEFAULT_WINVNE_PATH) ? read_winvne() : nothing
        ame = isfile(DEFAULT_AME_PATH) ? read_ame_masses() : nothing
        paper = ReactionRateTable[]
        isfile(DEFAULT_ILIADIS2001_PATH) && append!(paper, read_iliadis2001_rates(; partition_functions=pf, mass_excesses=ame))
        isfile(DEFAULT_NACRE_PATH) && append!(paper, read_nacre_rates(; partition_functions=pf, mass_excesses=ame))
        return iliadis2002_rate_tables(sets; paper_tables=paper, kwargs...)
    end
    return iliadis2002_rate_tables(sets; kwargs...)
end

