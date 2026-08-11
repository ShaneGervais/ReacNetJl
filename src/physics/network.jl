# Reaction and ReactionNetwork: the compiled network structure.

"""
    Reaction

A network reaction built from a `ReactionRateTable`.

The first version stores only the reaction participants and a STARLIB rate table.
Later we can add screening, reverse-rate handling, and uncertainty sampling here
without changing the network RHS interface.
"""
struct Reaction
    reactants::Vector{String}
    products::Vector{String}
    rate_table::ReactionRateTable
end

Reaction(table::ReactionRateTable) = Reaction(table.reactants, table.products, table)

struct CompiledReaction
    reactant_indices::Vector{Int}
    product_indices::Vector{Int}
    reactant_species_indices::Vector{Int}
    reactant_species_counts::Vector{Int}
    product_species_indices::Vector{Int}
    product_species_counts::Vector{Int}
    stoichiometric_delta::Vector{Float64}
    symmetry_factor::Float64
    nreactants::Int
    # Sum of Zi*Zj over charged reactant pairs, precomputed so weak screening
    # needs no species parsing in the solver hot loop.
    charge_pair_sum::Float64
    # Chained (Z1, A1, Z2, A2) screening pairs: one pair for two-body
    # reactions, and sequential compound pairs for three-body reactions
    # (e.g. triple-alpha screens he4+he4, then be8+he4). Pairs involving a
    # neutral particle are omitted.
    screening_pairs::Vector{NTuple{4,Float64}}
end

"""
    ReactionNetwork

A single-zone nuclear reaction network.

`species` fixes the order of the abundance vector `Y`. `species_index` maps each
species name to its position in `Y`, and `reactions` stores the reactions that
contribute to `dY/dt`. `jacobian_sparsity` is the structural (value-free)
nonzero pattern of `d(dY/dt)/dY`, precomputed once here rather than
rediscovered every solver step (feature spec Tier 0 #2); see
`_jacobian_sparsity_pattern` and `sparse_jacobian_prototype`.
"""
struct ReactionNetwork
    species::Vector{String}
    species_info::Vector{Species}
    species_index::Dict{String,Int}
    reactions::Vector{Reaction}
    compiled_reactions::Vector{CompiledReaction}
    jacobian_sparsity::SparseMatrixCSC{Bool,Int}
end

function ReactionNetwork(species::AbstractVector{<:AbstractString}, reactions::AbstractVector{Reaction})
    normalized_species = normalize_species_name.(species)
    length(unique(normalized_species)) == length(normalized_species) || throw(ArgumentError("network species must be unique"))

    species_info = species_from_name.(normalized_species)
    species_index = Dict(name => i for (i, name) in pairs(normalized_species))
    normalized_reactions = Reaction[]

    for reaction in reactions
        reactants = normalize_species_name.(reaction.reactants)
        products = normalize_species_name.(reaction.products)

        for name in Iterators.flatten((reactants, products))
            haskey(species_index, name) || throw(ArgumentError("reaction contains species '$name' that is missing from the network"))
        end

        push!(normalized_reactions, Reaction(reactants, products, reaction.rate_table))
    end

    compiled_reactions = [_compile_reaction(reaction, species_index, length(normalized_species)) for reaction in normalized_reactions]
    jacobian_sparsity = _jacobian_sparsity_pattern(length(normalized_species), compiled_reactions)
    return ReactionNetwork(normalized_species, species_info, species_index, normalized_reactions, compiled_reactions, jacobian_sparsity)
end

# Structural nonzero pattern of the analytic Jacobian d(dY_i/dt)/dY_j: for
# reaction r, every reactant column j (compiled.reactant_species_indices)
# potentially affects every row i the reaction touches, i.e. every reactant
# and product species of that same reaction (see `_cached_network_jacobian!`,
# whose dense loop has exactly this (reaction, j, i) structure). The diagonal
# is always included even for species with no local self-coupling, since the
# backward-Euler Newton matrix `I - dt*J` needs a diagonal entry to represent
# regardless of whether the network Jacobian itself happens to populate it.
function _jacobian_sparsity_pattern(nspecies::Int, compiled_reactions::Vector{CompiledReaction})
    entries = Set{Tuple{Int,Int}}()
    for i in 1:nspecies
        push!(entries, (i, i))
    end

    for compiled in compiled_reactions
        affected_rows = Iterators.flatten((compiled.reactant_species_indices, compiled.product_species_indices))
        for jindex in compiled.reactant_species_indices
            for iindex in affected_rows
                push!(entries, (iindex, jindex))
            end
        end
    end

    rows = [e[1] for e in entries]
    cols = [e[2] for e in entries]
    return sparse(rows, cols, trues(length(entries)), nspecies, nspecies)
end

function _factorial(n::Int)
    value = 1
    for i in 2:n
        value *= i
    end
    return value
end

# Symmetry factor prod_j (nu^react_j)! from the flux equation in
# `reaction_flux`: for `k` copies of an identical reactant species (e.g.
# `p+p`, `he4+he4+he4` for triple-alpha), the same physical collision is
# counted `k!` times by the naive `prod Y_j^nu_j` product, so dividing by
# `k!` per distinct reactant species recovers the correct flux. Reactions
# with all-distinct reactants have every count equal to 1 and this factor is
# 1, matching ordinary two-body reactions with no correction needed.
function _symmetry_factor(species::Vector{String})
    factor = 1
    for i in eachindex(species)
        count = 0
        for j in eachindex(species)
            species[j] == species[i] && (count += 1)
        end

        if findfirst(==(species[i]), species) == i
            factor *= _factorial(count)
        end
    end
    return factor
end

function _species_index(index::AbstractDict, name::String)
    haskey(index, name) && return index[name]
    throw(ArgumentError("species '$name' is missing from the species index"))
end

function _index_counts(indices::Vector{Int})
    unique_indices = Int[]
    counts = Int[]
    for index in indices
        position = findfirst(==(index), unique_indices)
        if position === nothing
            push!(unique_indices, index)
            push!(counts, 1)
        else
            counts[position] += 1
        end
    end
    return unique_indices, counts
end

# Precompute per-reaction bookkeeping once at network-construction time so the
# solver's per-timestep RHS/Jacobian evaluation only does array lookups, not
# repeated species-name parsing or charge/mass arithmetic:
#   - stoichiometric_delta: the dY_i/dt coefficient (nu^prod_i - nu^react_i)
#     for every species i, precomputed as a dense vector (see network_rhs).
#   - charge_pair_sum: sum_{i<j} Z_i Z_j over charged reactant pairs, the
#     Coulomb-barrier input to `weak_screening_multiplier`.
#   - screening_pairs: sequential (Z1,A1,Z2,A2) pairs for Chugunov-style
#     screening of reactions with more than two reactants -- e.g. triple-alpha
#     is screened as he4+he4 forming an effective 8Be-like pair, then that
#     compound charge/mass screened again against the third he4.
function _compile_reaction(reaction::Reaction, species_index::AbstractDict, nspecies::Int)
    reactant_indices = [_species_index(species_index, name) for name in reaction.reactants]
    product_indices = [_species_index(species_index, name) for name in reaction.products]
    reactant_species_indices, reactant_species_counts = _index_counts(reactant_indices)
    product_species_indices, product_species_counts = _index_counts(product_indices)
    stoichiometric_delta = zeros(Float64, nspecies)

    for (index, count) in zip(reactant_species_indices, reactant_species_counts)
        stoichiometric_delta[index] -= count
    end
    for (index, count) in zip(product_species_indices, product_species_counts)
        stoichiometric_delta[index] += count
    end

    charge_pair_sum = 0.0
    screening_pairs = NTuple{4,Float64}[]
    if length(reaction.reactants) >= 2
        charges = Int[]
        masses = Int[]
        for name in reaction.reactants
            species = try
                species_from_name(name)
            catch
                Species(String(name), 0, 0)
            end
            push!(charges, species.Z)
            push!(masses, species.A)
        end

        for i in 1:(length(charges) - 1), j in (i + 1):length(charges)
            if charges[i] > 0 && charges[j] > 0
                charge_pair_sum += charges[i] * charges[j]
            end
        end

        accumulated_Z = charges[1]
        accumulated_A = masses[1]
        for k in 2:length(charges)
            if accumulated_Z > 0 && charges[k] > 0
                push!(screening_pairs, (Float64(accumulated_Z), Float64(accumulated_A), Float64(charges[k]), Float64(masses[k])))
            end
            accumulated_Z += charges[k]
            accumulated_A += masses[k]
        end
    end

    return CompiledReaction(
        reactant_indices,
        product_indices,
        reactant_species_indices,
        reactant_species_counts,
        product_species_indices,
        product_species_counts,
        stoichiometric_delta,
        Float64(_symmetry_factor(reaction.reactants)),
        length(reaction.reactants),
        charge_pair_sum,
        screening_pairs,
    )
end

function _select_rate_table(matches::Vector{ReactionRateTable}, label::AbstractString, on_multiple::Symbol)
    isempty(matches) && throw(ArgumentError("no STARLIB rate found for reaction '$label'"))

    if length(matches) == 1
        return only(matches)
    elseif on_multiple == :first
        return first(matches)
    elseif on_multiple == :error
        sources = join([m.source for m in matches], ", ")
        throw(ArgumentError("multiple STARLIB rates found for reaction '$label' with sources: $sources; pass `source=...` or `on_multiple=:first`"))
    else
        throw(ArgumentError("unsupported on_multiple=$on_multiple; use :error or :first"))
    end
end

"""
    reaction_from_label(tables, label; source=nothing, on_multiple=:error)

Find a STARLIB rate table by reaction label and wrap it as a `Reaction`.
If multiple entries match, choose a `source`, or set `on_multiple=:first`.
"""
function reaction_from_label(
    tables::AbstractVector{ReactionRateTable},
    label::AbstractString;
    source=nothing,
    on_multiple::Symbol=:error,
)
    matches = find_rate(tables, label; source=source)
    return Reaction(_select_rate_table(matches, label, on_multiple))
end

function _infer_species_from_reactions(reactions::AbstractVector{Reaction})
    species = String[]
    seen = Set{String}()

    for reaction in reactions
        for name in Iterators.flatten((reaction.reactants, reaction.products))
            if !in(name, seen)
                push!(species, name)
                push!(seen, name)
            end
        end
    end

    return species
end

"""
    network_from_tables(tables; species=nothing)

Build a `ReactionNetwork` directly from a vector of `ReactionRateTable`s (one
reaction per table, after de-duplicating identical reactant/product/chapter
entries via `_unique_reaction_tables`). This is the path used by `run_ppn` and
the nova-network selectors (`select_h_ca_reaction_tables`,
`add_reverse_reaction_tables`), which already produce a specific table list to
turn into a network. If `species` is not supplied, the species list is
inferred from the union of every reaction's reactants and products.
"""
function network_from_tables(tables::AbstractVector{ReactionRateTable}; species=nothing)
    reactions = [Reaction(table) for table in _unique_reaction_tables(tables)]
    network_species = species === nothing ? _infer_species_from_reactions(reactions) : collect(species)
    return ReactionNetwork(network_species, reactions)
end

"""
    network_from_labels(tables, labels; species=nothing, source=nothing, on_multiple=:error)

Build a `ReactionNetwork` directly from human-readable reaction labels (e.g.
`"18F(p,α)15O"`), looking each one up in `tables` via `reaction_from_label`.
This is the label-driven, interactive-experimentation path (used by
`solve_single_zone`); `network_from_tables` is the corresponding path for
programmatically-selected table lists. If `species` is not supplied, the
network species list is inferred from the selected reactions. If a label
matches more than one rate table (e.g. both a STARLIB and a REACLIB entry),
pass `source` to disambiguate or `on_multiple=:first` to take the first match.
"""
function network_from_labels(
    tables::AbstractVector{ReactionRateTable},
    labels::AbstractVector{<:AbstractString};
    species=nothing,
    source=nothing,
    on_multiple::Symbol=:error,
)
    reactions = [reaction_from_label(tables, label; source=source, on_multiple=on_multiple) for label in labels]
    network_species = species === nothing ? _infer_species_from_reactions(reactions) : collect(species)
    return ReactionNetwork(network_species, reactions)
end

