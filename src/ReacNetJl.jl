module ReacNetJl

using LinearAlgebra
using Random
using PrecompileTools: @setup_workload, @compile_workload

export Species,
    Trajectory,
    ReactionRateTable,
    ReaclibSet,
    Reaction,
    ReactionNetwork,
    species_from_name,
    abundance_from_mass_fraction,
    mass_fraction_from_abundance,
    normalize_species_name,
    parse_reaction_label,
    read_starlib,
    read_reaclib,
    read_winvne,
    read_tabulated_rates,
    read_iliadis2001_rates,
    read_nacre_rates,
    reaction_q_value,
    PartitionFunctionTable,
    reaclib_rate,
    reaclib_rate_tables,
    iliadis2002_rate_tables,
    read_trajectory,
    read_initial_abundances,
    trajectory_profiles,
    first_cooling_threshold_time,
    first_mass_fraction_threshold_crossing,
    starlib_chapter_report,
    find_rate,
    find_reverse_rate,
    reaction_from_label,
    generated_detailed_balance_reverse_table,
    add_reverse_reaction_tables,
    select_h_ca_reaction_tables,
    select_decay_reaction_tables,
    decay_mass_fractions,
    network_from_tables,
    network_from_labels,
    weak_screening_multiplier,
    chugunov_screening_multiplier,
    abundances_from_mass_fractions,
    mass_fractions_from_abundances,
    mass_fraction_history,
    total_mass_fraction,
    total_mass_fraction_history,
    mass_fraction_drift,
    abundance_diagnostics,
    interpolate_rate,
    interpolate_factor_uncertainty,
    sampled_interpolate_rate,
    reaction_string,
    reaction_flux,
    reaction_fluxes,
    reaction_flux_history,
    integrated_fluxes,
    energy_generation_rate,
    energy_generation_history,
    integrated_energy_generation,
    species_flux_balance,
    reaction_edges,
    reaction_conservation,
    network_validation_report,
    network_rhs,
    solve_network,
    solve_network_adaptive,
    solve_network_fbdf,
    run_monte_carlo,
    solve_single_zone,
    run_ppn

"""
    Species

A nuclear species identified by charge number `Z`, mass number `A`, and a display
`name`, such as `Species("f18", 9, 18)`.
"""
struct Species
    name::String
    Z::Int
    A::Int
end

"""
    Trajectory

A single-zone thermodynamic trajectory with time in seconds, temperature `T9` in
GK, and density `rho` in g cm^-3.
"""
struct Trajectory
    time::Vector{Float64}
    T9::Vector{Float64}
    rho::Vector{Float64}
end

"""
    ReactionRateTable

A temperature-dependent thermonuclear reaction rate from STARLIB.

`T9` is the temperature grid in GK. `rate` is the recommended STARLIB rate.
`factor_uncertainty` is the STARLIB multiplicative factor uncertainty.
"""
struct ReactionRateTable
    chapter::Int
    reactants::Vector{String}
    products::Vector{String}
    source::String
    q_value::Float64
    T9::Vector{Float64}
    rate::Vector{Float64}
    factor_uncertainty::Vector{Float64}
end

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
contribute to `dY/dt`.
"""
struct ReactionNetwork
    species::Vector{String}
    species_info::Vector{Species}
    species_index::Dict{String,Int}
    reactions::Vector{Reaction}
    compiled_reactions::Vector{CompiledReaction}
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
    return ReactionNetwork(normalized_species, species_info, species_index, normalized_reactions, compiled_reactions)
end

const DEFAULT_STARLIB_PATH = joinpath(dirname(@__DIR__), "data", "starlib.dat")
const LEGACY_STARLIB_PATH = joinpath(dirname(@__DIR__), "starlib.dat")
const DEFAULT_REACLIB_PATH = joinpath(dirname(@__DIR__), "data", "reaclib_v1.0.dat")
const DEFAULT_REACLIB_IL01_PATH = joinpath(dirname(@__DIR__), "data", "reaclib_il01.dat")
const DEFAULT_REACLIB_NACR_PATH = joinpath(dirname(@__DIR__), "data", "reaclib_nacr.dat")
const DEFAULT_ILIADIS2001_PATH = joinpath(dirname(@__DIR__), "data", "iliadis2001_rates.dat")
const DEFAULT_NACRE_PATH = joinpath(dirname(@__DIR__), "data", "nacre_rates.dat")
const STARLIB_ROWS_PER_REACTION = 60

function _default_starlib_path()
    isfile(DEFAULT_STARLIB_PATH) && return DEFAULT_STARLIB_PATH
    isfile(LEGACY_STARLIB_PATH) && return LEGACY_STARLIB_PATH
    return DEFAULT_STARLIB_PATH
end

function _default_reaclib_path()
    isfile(DEFAULT_REACLIB_PATH) || error(
        "REACLIB library not found at $(DEFAULT_REACLIB_PATH); run data/download_rates.sh to fetch it",
    )
    return DEFAULT_REACLIB_PATH
end
const AVOGADRO = 6.02214076e23
const MEV_TO_ERG = 1.602176634e-6

const _PARTICLE_ALIASES = Dict(
    "p" => "p",
    "proton" => "p",
    "h1" => "p",
    "n" => "n",
    "neutron" => "n",
    "d" => "d",
    "h2" => "d",
    "t" => "t",
    "h3" => "t",
    "a" => "he4",
    "alpha" => "he4",
    "he4" => "he4",
    "g" => "gamma",
    "gamma" => "gamma",
)

const _SPECIAL_SPECIES = Dict(
    "n" => Species("n", 0, 1),
    "p" => Species("p", 1, 1),
    "d" => Species("d", 1, 2),
    "t" => Species("t", 1, 3),
    "gamma" => Species("gamma", 0, 0),
)

const _ELEMENT_SYMBOLS = [
    "h", "he", "li", "be", "b", "c", "n", "o", "f", "ne",
    "na", "mg", "al", "si", "p", "s", "cl", "ar", "k", "ca",
    "sc", "ti", "v", "cr", "mn", "fe", "co", "ni", "cu", "zn",
    "ga", "ge", "as", "se", "br", "kr", "rb", "sr", "y", "zr",
    "nb", "mo", "tc", "ru", "rh", "pd", "ag", "cd", "in", "sn",
    "sb", "te", "i", "xe", "cs", "ba", "la", "ce", "pr", "nd",
    "pm", "sm", "eu", "gd", "tb", "dy", "ho", "er", "tm", "yb",
    "lu", "hf", "ta", "w", "re", "os", "ir", "pt", "au", "hg",
    "tl", "pb", "bi", "po", "at", "rn", "fr", "ra", "ac", "th",
    "pa", "u", "np", "pu", "am", "cm", "bk", "cf", "es", "fm",
    "md", "no", "lr", "rf", "db", "sg", "bh", "hs", "mt", "ds",
    "rg", "cn", "nh", "fl", "mc", "lv", "ts", "og",
]

const _ELEMENT_Z = Dict(symbol => Z for (Z, symbol) in pairs(_ELEMENT_SYMBOLS))

"""
    normalize_species_name(name)

Normalize common isotope notation to STARLIB-style lowercase names.

Examples:
- `"18F" -> "f18"`
- `"4He" -> "he4"`
- `"α" -> "he4"`
- `"γ" -> "gamma"`
"""
function normalize_species_name(name::AbstractString)
    raw = lowercase(strip(name))
    raw == "al-6" && return "al26"

    s = raw
    s = replace(s, " " => "", "-" => "", "_" => "")
    s = replace(s, "α" => "alpha", "γ" => "gamma")

    if haskey(_PARTICLE_ALIASES, s)
        return _PARTICLE_ALIASES[s]
    end

    m = match(r"^(\d+)([a-z]+)\*$", s)
    if m !== nothing
        mass_number = m.captures[1]
        symbol = m.captures[2]
        return symbol * "*" * last(mass_number)
    end

    m = match(r"^(\d+)alm$", s)
    if m !== nothing
        mass_number = m.captures[1]
        return "al*" * last(mass_number)
    end

    m = match(r"^(\d+)alg$", s)
    if m !== nothing
        return "al" * m.captures[1]
    end

    m = match(r"^(\d+)([a-z]+)$", s)
    if m !== nothing
        return m.captures[2] * m.captures[1]
    end

    return s
end

"""
    species_from_name(name)

Create a `Species` object from a normalized or common isotope name.

Examples:
- `species_from_name("18F") == Species("f18", 9, 18)`
- `species_from_name("he4") == Species("he4", 2, 4)`
- `species_from_name("p") == Species("p", 1, 1)`
"""
function species_from_name(name::AbstractString)
    normalized = normalize_species_name(name)

    if haskey(_SPECIAL_SPECIES, normalized)
        return _SPECIAL_SPECIES[normalized]
    end

    m = match(r"^([a-z]+)(\d+)$", normalized)
    if m === nothing
        star_match = match(r"^([a-z]+)\*(\d+)$", normalized)
        if star_match !== nothing
            symbol = star_match.captures[1]
            mass_suffix = star_match.captures[2]
            haskey(_ELEMENT_Z, symbol) || throw(ArgumentError("unknown element symbol '$symbol' in species '$name'"))
            A = symbol == "al" && mass_suffix == "6" ? 26 : parse(Int, mass_suffix)
            return Species(normalized, _ELEMENT_Z[symbol], A)
        end
    end
    m === nothing && throw(ArgumentError("could not infer species information from '$name'"))

    symbol = m.captures[1]
    mass_suffix = m.captures[2]
    A = symbol == "al" && mass_suffix == "6" ? 26 : parse(Int, mass_suffix)
    haskey(_ELEMENT_Z, symbol) || throw(ArgumentError("unknown element symbol '$symbol' in species '$name'"))

    return Species(normalized, _ELEMENT_Z[symbol], A)
end

#=
    parse_reaction_label(label)

Parse a reaction label of the form `target(projectile,ejectile)product`.
Returns `(reactants, products)`, where both entries are normalized species-name
vectors.

Examples:
- `parse_reaction_label("18F(p,α)15O") == (["p", "f18"], ["he4", "o15"])`
- `parse_reaction_label("18F(p,γ)19Ne") == (["p", "f18"], ["ne19"])`
- `parse_reaction_label("p(p,eν)d") == (["p", "p"], ["d"])`
=#
function _is_beta_token(token::AbstractString)
    s = lowercase(strip(token))
    s = replace(s, " " => "", "β" => "beta")
    return s in ("beta+", "b+", "β+", "betaplus", "ec", "electroncapture", "beta-", "b-", "β-", "betaminus")
end

function _parse_particle_multiplicity(token::AbstractString)
    s = lowercase(strip(token))
    s = replace(s, " " => "", "α" => "alpha")
    m = match(r"^([2-9])([a-z\+]+)$", s)
    m === nothing && return nothing

    count = parse(Int, m.captures[1])
    particle = m.captures[2]
    particle in ("p", "n", "d", "t", "alpha") || return nothing
    normalized = normalize_species_name(particle)
    return fill(normalized, count)
end

function _is_non_nuclear_ejectile_token(token::AbstractString)
    s = lowercase(strip(token))
    s = replace(s, " " => "", "ν" => "nu")
    return s in ("enu", "e+nu", "e-nu", "nu", "neutrino")
end

function _parse_reaction_token(token::AbstractString; allow_multiplicity::Bool=false)
    _is_non_nuclear_ejectile_token(token) && return String[]

    if allow_multiplicity
        particles = _parse_particle_multiplicity(token)
        particles !== nothing && return particles
    end

    normalized = normalize_species_name(token)
    normalized == "gamma" && return String[]
    return [normalized]
end

function parse_reaction_label(label::AbstractString)
    decay_match = match(r"^\s*([^\(]+)\(([^,\)]+)\)(.+?)\s*$", label)
    if decay_match !== nothing && _is_beta_token(decay_match.captures[2])
        parent = normalize_species_name(decay_match.captures[1])
        products = _parse_reaction_token(decay_match.captures[3]; allow_multiplicity=true)
        return [parent], products
    end

    m = match(r"^\s*([^\(]+)\(([^,]+),([^\)]+)\)(.+?)\s*$", label)
    m === nothing && throw(ArgumentError("reaction label must look like target(projectile,ejectile)product or parent(β+)daughter"))

    target = normalize_species_name(m.captures[1])
    projectile = normalize_species_name(m.captures[2])
    ejectiles = _parse_reaction_token(m.captures[3]; allow_multiplicity=true)
    product = normalize_species_name(m.captures[4])

    reactants = [projectile, target]
    products = vcat(ejectiles, [product])
    return reactants, products
end

#=
    abundance_from_mass_fraction(X, A)

Convert mass fraction `Xᵢ` to abundance `Yᵢ = Xᵢ / Aᵢ`.
=#
abundance_from_mass_fraction(X::Real, A::Integer) = X / A

#=
    mass_fraction_from_abundance(Y, A)

Convert abundance `Yᵢ` to mass fraction `Xᵢ = Aᵢ Yᵢ`.
=#
mass_fraction_from_abundance(Y::Real, A::Integer) = A * Y

function _normalize_mass_fraction_keys(X::AbstractDict)
    normalized = Dict{String,Float64}()
    for (name, value) in X
        normalized_name = normalize_species_name(string(name))
        haskey(normalized, normalized_name) && throw(ArgumentError("duplicate mass fraction entry for species '$normalized_name'"))
        normalized[normalized_name] = Float64(value)
    end
    return normalized
end

#=
    abundances_from_mass_fractions(network, X; normalize=false, check_sum=false, atol=1e-8)

Convert a dictionary of mass fractions into an abundance vector ordered like
`network.species`.

Missing network species are assigned zero mass fraction. Extra species not in the
network raise an error. If `normalize=true`, mass fractions are divided by their
total before conversion.
=#
function abundances_from_mass_fractions(
    network::ReactionNetwork,
    X::AbstractDict;
    normalize::Bool=false,
    check_sum::Bool=false,
    atol::Real=1.0e-8,
)
    normalized = _normalize_mass_fraction_keys(X)

    for name in keys(normalized)
        haskey(network.species_index, name) || throw(ArgumentError("mass fraction provided for species '$name' that is not in the network"))
    end

    total = sum(values(normalized); init=0.0)
    if normalize
        total > 0.0 || throw(ArgumentError("cannot normalize mass fractions with non-positive total"))
    elseif check_sum && abs(total - 1.0) > atol
        throw(ArgumentError("mass fractions sum to $total, not 1 within atol=$atol"))
    end

    Y = zeros(Float64, length(network.species))
    for (i, species) in pairs(network.species_info)
        species.A > 0 || throw(ArgumentError("cannot convert mass fraction for species '$(species.name)' with A=$(species.A)"))
        X_i = get(normalized, species.name, 0.0)
        normalize && (X_i /= total)
        Y[i] = abundance_from_mass_fraction(X_i, species.A)
    end

    return Y
end

#=
    mass_fractions_from_abundances(network, Y)

Convert an abundance vector ordered like `network.species` into a dictionary of
mass fractions keyed by normalized species name.
=#
function mass_fractions_from_abundances(network::ReactionNetwork, Y::AbstractVector{<:Real})
    length(Y) == length(network.species) || throw(ArgumentError("Y length must match the number of network species"))

    X = Dict{String,Float64}()
    for (i, species) in pairs(network.species_info)
        species.A > 0 || throw(ArgumentError("cannot convert abundance for species '$(species.name)' with A=$(species.A)"))
        X[species.name] = mass_fraction_from_abundance(Y[i], species.A)
    end
    return X
end

#=
    mass_fraction_history(network, history)

Convert an abundance-history matrix into a mass-fraction-history matrix with
the same shape. Columns remain ordered like `network.species`.
=#
function mass_fraction_history(network::ReactionNetwork, history::AbstractMatrix{<:Real})
    size(history, 2) == length(network.species) || throw(ArgumentError("history column count must match the number of network species"))

    X_history = Matrix{Float64}(undef, size(history))
    for j in axes(history, 2)
        A = network.species_info[j].A
        A > 0 || throw(ArgumentError("cannot convert abundance for species '$(network.species_info[j].name)' with A=$A"))
        for i in axes(history, 1)
            X_history[i, j] = mass_fraction_from_abundance(history[i, j], A)
        end
    end
    return X_history
end

function first_mass_fraction_threshold_crossing(
    network::ReactionNetwork,
    times::AbstractVector{<:Real},
    history::AbstractMatrix{<:Real},
    species::AbstractString,
    threshold::Real;
    direction::Symbol=:down,
)
    length(times) == size(history, 1) || throw(ArgumentError("times length must match the number of history rows"))
    size(history, 2) == length(network.species) || throw(ArgumentError("history column count must match the number of network species"))
    direction in (:down, :up) || throw(ArgumentError("direction must be :down or :up"))

    limit = Float64(threshold)
    limit >= 0.0 || throw(ArgumentError("threshold must be non-negative"))
    name = normalize_species_name(species)
    haskey(network.species_index, name) || throw(ArgumentError("species '$species' is not present in the network"))
    index = network.species_index[name]
    A = network.species_info[index].A

    mass_fraction_at(row) = mass_fraction_from_abundance(history[row, index], A)
    initial = mass_fraction_at(firstindex(times))
    if (direction == :down && initial <= limit) || (direction == :up && initial >= limit)
        return (
            time=Float64(first(times)),
            state=copy(history[firstindex(times), :]),
            previous_index=firstindex(times),
            fraction=0.0,
            mass_fraction=initial,
        )
    end

    for i in firstindex(times):(lastindex(times)-1)
        x0 = mass_fraction_at(i)
        x1 = mass_fraction_at(i + 1)
        crossed = direction == :down ? (x0 >= limit && x1 <= limit) : (x0 <= limit && x1 >= limit)
        crossed || continue

        fraction = x0 == x1 ? 0.0 : (limit - x0) / (x1 - x0)
        stop_time = Float64(times[i]) + fraction * Float64(times[i+1] - times[i])
        state = (1.0 - fraction) .* history[i, :] .+ fraction .* history[i+1, :]
        return (
            time=stop_time,
            state=Vector{Float64}(state),
            previous_index=i,
            fraction=fraction,
            mass_fraction=limit,
        )
    end

    return nothing
end

#=
    total_mass_fraction(network, Y)

Return `sum_i A_i Y_i` for one abundance state.
=#
function total_mass_fraction(network::ReactionNetwork, Y::AbstractVector{<:Real})
    length(Y) == length(network.species) || throw(ArgumentError("Y length must match the number of network species"))

    total = 0.0
    for (i, species) in pairs(network.species_info)
        species.A > 0 || throw(ArgumentError("cannot convert abundance for species '$(species.name)' with A=$(species.A)"))
        total += mass_fraction_from_abundance(Y[i], species.A)
    end
    return total
end

#=
    total_mass_fraction_history(network, history)

Return the total mass fraction at every saved abundance-history row.
=#
function total_mass_fraction_history(network::ReactionNetwork, history::AbstractMatrix{<:Real})
    size(history, 2) == length(network.species) || throw(ArgumentError("history column count must match the number of network species"))

    totals = Vector{Float64}(undef, size(history, 1))
    for i in axes(history, 1)
        totals[i] = total_mass_fraction(network, view(history, i, :))
    end
    return totals
end

#=
    mass_fraction_drift(network, history)

Summarize total-mass-fraction drift over a history.
=#
function mass_fraction_drift(network::ReactionNetwork, history::AbstractMatrix{<:Real})
    totals = total_mass_fraction_history(network, history)
    initial = first(totals)
    final = last(totals)
    deviations = abs.(totals .- initial)
    return (
        initial=initial,
        final=final,
        drift=final - initial,
        max_abs_drift=maximum(deviations),
        min_total=minimum(totals),
        max_total=maximum(totals),
        totals=totals,
    )
end

function _matrix_minimum_location(values::AbstractMatrix{<:Real})
    min_value = Inf
    min_row = 0
    min_col = 0
    for row in axes(values, 1)
        for col in axes(values, 2)
            value = Float64(values[row, col])
            if value < min_value
                min_value = value
                min_row = row
                min_col = col
            end
        end
    end
    return min_value, min_row, min_col
end

#=
    abundance_diagnostics(network, history)

Report positivity-oriented diagnostics for an abundance history, including the
minimum abundance and minimum mass fraction with species/time indices.
=#
function abundance_diagnostics(network::ReactionNetwork, history::AbstractMatrix{<:Real})
    size(history, 2) == length(network.species) || throw(ArgumentError("history column count must match the number of network species"))

    min_Y, min_Y_row, min_Y_col = _matrix_minimum_location(history)
    X_history = mass_fraction_history(network, history)
    min_X, min_X_row, min_X_col = _matrix_minimum_location(X_history)

    return (
        min_abundance=min_Y,
        min_abundance_time_index=min_Y_row,
        min_abundance_species=network.species[min_Y_col],
        min_mass_fraction=min_X,
        min_mass_fraction_time_index=min_X_row,
        min_mass_fraction_species=network.species[min_X_col],
        has_negative_abundance=min_Y < 0.0,
        has_negative_mass_fraction=min_X < 0.0,
    )
end

function _split_starlib_species(chapter::Int, species::Vector{String})
    arities = Dict(
        1 => (1, 1),
        2 => (1, 2),
        3 => (1, 3),
        4 => (2, 1),
        5 => (2, 2),
        6 => (2, 3),
        7 => (2, 4),
        8 => (3, 1),
        9 => (3, 2),
        10 => (4, 2),
    )

    if haskey(arities, chapter)
        nreactants, nproducts = arities[chapter]
        if length(species) == nreactants + nproducts
            return species[1:nreactants], species[(nreactants + 1):end]
        end
    end

    # Conservative fallback: keep the raw STARLIB order if this chapter is not
    # supported yet. We will expand this as the network grows.
    return species, String[]
end

function _supported_starlib_layout(chapter::Int, nspecies::Int)
    return (chapter == 1 && nspecies == 2) ||
           (chapter in (2, 4) && nspecies == 3) ||
           (chapter in (3, 5, 8) && nspecies == 4) ||
           (chapter in (6, 9) && nspecies == 5) ||
           (chapter in (7, 10) && nspecies == 6)
end

function _supported_rate_table(table::ReactionRateTable)
    return !isempty(table.reactants) && !isempty(table.products)
end

function _valid_nuclear_bookkeeping(table::ReactionRateTable)
    _supported_rate_table(table) || return false
    return try
        reaction_conservation(Reaction(table)).valid_nuclear_bookkeeping
    catch
        false
    end
end

"""
    starlib_chapter_report(tables)

Summarize which STARLIB chapter layouts were parsed into supported reactant and
product bookkeeping. Unsupported rows are kept by `read_starlib`, but they are
not suitable for network construction until their chapter layout is implemented.
"""
function starlib_chapter_report(tables::AbstractVector{ReactionRateTable})
    supported_by_chapter = Dict{Int,Int}()
    unsupported_by_chapter = Dict{Int,Int}()
    unsupported_examples = ReactionRateTable[]

    for table in tables
        if _supported_rate_table(table)
            supported_by_chapter[table.chapter] = get(supported_by_chapter, table.chapter, 0) + 1
        else
            unsupported_by_chapter[table.chapter] = get(unsupported_by_chapter, table.chapter, 0) + 1
            length(unsupported_examples) < 10 && push!(unsupported_examples, table)
        end
    end

    return (
        total=length(tables),
        supported=sum(values(supported_by_chapter); init=0),
        unsupported=sum(values(unsupported_by_chapter); init=0),
        supported_by_chapter=supported_by_chapter,
        unsupported_by_chapter=unsupported_by_chapter,
        unsupported_examples=unsupported_examples,
    )
end

#=
    read_starlib(path=DEFAULT_STARLIB_PATH)

Read a STARLIB v6-style `.dat` file into `ReactionRateTable` entries.

For now, this implements the file format and common chapter layouts needed for
simple one-zone experiments. More chapter-specific reaction bookkeeping can be
added as we expand the physical network.
=#
function _validate_trajectory(time::Vector{Float64}, T9::Vector{Float64}, rho::Vector{Float64})
    length(time) == length(T9) == length(rho) || throw(ArgumentError("trajectory columns must have the same length"))
    length(time) >= 2 || throw(ArgumentError("trajectory must contain at least two rows"))

    for i in 1:(length(time)-1)
        time[i+1] > time[i] || throw(ArgumentError("trajectory time values must be strictly increasing"))
    end

    any(<=(0.0), T9) && throw(ArgumentError("trajectory T9 values must be positive"))
    any(<=(0.0), rho) && throw(ArgumentError("trajectory rho values must be positive"))
end

#=
    read_trajectory(path)

Read a whitespace- or comma-separated trajectory file with columns:

    time_s  T9  rho

Lines beginning with `#` and blank lines are ignored. Metadata assignments are
also supported before the numeric rows:

    AGEUNIT = SEC | YRS
    TUNIT   = T9K | T8K
    RHOUNIT = CGS | LOG
=#
function read_trajectory(path::AbstractString)
    time = Float64[]
    T9 = Float64[]
    rho = Float64[]
    age_unit = "SEC"
    temperature_unit = "T9K"
    density_unit = "CGS"

    open(path, "r") do io
        for raw_line in eachline(io)
            line = strip(split(raw_line, '#'; limit=2)[1])
            isempty(line) && continue

            if occursin("=", line)
                key_value = split(line, '='; limit=2)
                key = uppercase(strip(key_value[1]))
                value = uppercase(strip(key_value[2]))
                if key == "AGEUNIT"
                    age_unit = value
                elseif key == "TUNIT"
                    temperature_unit = value
                elseif key == "RHOUNIT"
                    density_unit = value
                end
                continue
            end

            fields = split(replace(line, ',' => ' '))
            length(fields) >= 3 || throw(ArgumentError("trajectory row must contain at least three columns: $raw_line"))

            time_value = parse(Float64, fields[1])
            T_value = parse(Float64, fields[2])
            rho_value = parse(Float64, fields[3])

            if age_unit == "SEC"
                push!(time, time_value)
            elseif age_unit == "YRS"
                push!(time, time_value * 365.25 * 24.0 * 60.0 * 60.0)
            else
                throw(ArgumentError("unsupported trajectory AGEUNIT=$age_unit; use SEC or YRS"))
            end

            if temperature_unit == "T9K"
                push!(T9, T_value)
            elseif temperature_unit == "T8K"
                push!(T9, T_value / 10.0)
            else
                throw(ArgumentError("unsupported trajectory TUNIT=$temperature_unit; use T9K or T8K"))
            end

            if density_unit == "CGS"
                push!(rho, rho_value)
            elseif density_unit == "LOG"
                push!(rho, 10.0^rho_value)
            else
                throw(ArgumentError("unsupported trajectory RHOUNIT=$density_unit; use CGS or LOG"))
            end
        end
    end

    _validate_trajectory(time, T9, rho)
    return Trajectory(time, T9, rho)
end

function _parse_initial_abundance_species(fields::Vector{SubString{String}}, raw_line::AbstractString)
    if length(fields) == 3 && uppercase(fields[2]) == "PROT"
        return "p", parse(Float64, fields[3])
    elseif length(fields) == 4
        return normalize_species_name(string(fields[2], fields[3])), parse(Float64, fields[4])
    elseif length(fields) == 3
        return normalize_species_name(fields[2]), parse(Float64, fields[3])
    end

    throw(ArgumentError("initial abundance row must look like `Z sym A X`, `Z symA X`, or `1 PROT X`: $raw_line"))
end

#=
    read_initial_abundances(path; normalize=false)

Read an initial-abundance mass-fraction file with rows such as:

    1 PROT  3.5e-1
    6 c 12  1.0e-2
    26 fe56 6.0e-4

Returns a dictionary keyed by normalized species name. If `normalize=true`, all
mass fractions are divided by the file total.
=#
function read_initial_abundances(path::AbstractString; normalize::Bool=false)
    X = Dict{String,Float64}()

    open(path, "r") do io
        for raw_line in eachline(io)
            line = strip(split(raw_line, '#'; limit=2)[1])
            isempty(line) && continue
            fields = split(line)
            name, value = _parse_initial_abundance_species(fields, raw_line)
            haskey(X, name) && throw(ArgumentError("duplicate initial abundance entry for species '$name'"))
            X[name] = value
        end
    end

    normalize || return X
    total = sum(values(X); init=0.0)
    total > 0.0 || throw(ArgumentError("cannot normalize initial abundances with non-positive total"))
    return Dict(name => value / total for (name, value) in X)
end

function _linear_profile(x::Vector{Float64}, y::Vector{Float64}, value::Real)
    v = Float64(value)
    first(x) <= v <= last(x) || throw(ArgumentError("profile value $v is outside the trajectory range $(first(x))–$(last(x))"))

    i = searchsortedlast(x, v)
    if i == length(x) || x[i] == v
        return y[i]
    end

    weight = (v - x[i]) / (x[i+1] - x[i])
    return (1.0 - weight) * y[i] + weight * y[i+1]
end

#=
    trajectory_profiles(trajectory)

Return callable density and temperature profiles as a named tuple:

    profiles = trajectory_profiles(traj)
    rho_t = profiles.rho(t)
    T9_t = profiles.T9(t)
=#
function trajectory_profiles(trajectory::Trajectory)
    rho_profile(t) = _linear_profile(trajectory.time, trajectory.rho, t)
    T9_profile(t) = _linear_profile(trajectory.time, trajectory.T9, t)
    return (rho=rho_profile, T9=T9_profile)
end

function first_cooling_threshold_time(trajectory::Trajectory, threshold_T9::Real)
    threshold = Float64(threshold_T9)
    threshold > 0.0 || throw(ArgumentError("threshold_T9 must be positive"))

    peak_index = argmax(trajectory.T9)
    for i in peak_index:(length(trajectory.time)-1)
        T0 = trajectory.T9[i]
        T1 = trajectory.T9[i+1]
        t0 = trajectory.time[i]
        t1 = trajectory.time[i+1]

        T0 == threshold && return t0
        if T0 > threshold && T1 <= threshold
            fraction = (threshold - T0) / (T1 - T0)
            return t0 + fraction * (t1 - t0)
        end
    end

    return nothing
end

function read_starlib(path::AbstractString=_default_starlib_path(); warn_unsupported::Bool=false)
    tables = ReactionRateTable[]
    unsupported_counts = Dict{Int,Int}()

    open(path, "r") do io
        line_number = 0
        while !eof(io)
            header = strip(readline(io))
            line_number += 1
            isempty(header) && continue

            fields = split(header)
            length(fields) >= 4 || error("Malformed STARLIB header at line $line_number: $header")

            chapter = parse(Int, fields[1])
            source = fields[end-1]
            q_value = parse(Float64, fields[end])
            species = normalize_species_name.(fields[2:end-2])
            if !_supported_starlib_layout(chapter, length(species))
                unsupported_counts[chapter] = get(unsupported_counts, chapter, 0) + 1
            end
            reactants, products = _split_starlib_species(chapter, species)

            T9 = Float64[]
            rate = Float64[]
            factor_uncertainty = Float64[]
            sizehint!(T9, STARLIB_ROWS_PER_REACTION)
            sizehint!(rate, STARLIB_ROWS_PER_REACTION)
            sizehint!(factor_uncertainty, STARLIB_ROWS_PER_REACTION)

            for _ in 1:STARLIB_ROWS_PER_REACTION
                eof(io) && error("Unexpected end of file after STARLIB header at line $line_number")
                row = split(strip(readline(io)))
                line_number += 1
                length(row) >= 3 || error("Malformed STARLIB data row at line $line_number")

                push!(T9, parse(Float64, row[1]))
                push!(rate, parse(Float64, row[2]))
                push!(factor_uncertainty, parse(Float64, row[3]))
            end

            push!(tables, ReactionRateTable(chapter, reactants, products, source, q_value, T9, rate, factor_uncertainty))
        end
    end

    if warn_unsupported && !isempty(unsupported_counts)
        summary = join(["chapter $chapter: $count" for (chapter, count) in sort(collect(unsupported_counts))], ", ")
        @warn "Unsupported STARLIB chapter layouts were kept as raw reactants with empty products: $summary"
    end

    return tables
end

"""
    ReaclibSet

One REACLIB fit set: a single 7-coefficient term of a reaction rate. A REACLIB
rate is the sum of all sets sharing the same reaction and set label, e.g. a
non-resonant plus one or more resonant terms.

`resonance` is the REACLIB flag character (' ' or 'n' non-resonant, 'r'
resonant, 'w' weak). `reverse` marks detailed-balance reverse fits ('v' flag),
which REACLIB derives without partition functions.
"""
struct ReaclibSet
    chapter::Int
    reactants::Vector{String}
    products::Vector{String}
    label::String
    resonance::Char
    reverse::Bool
    q_value::Float64
    coefficients::NTuple{7,Float64}
end

# REACLIB2 chapter layouts (reactants, products). Chapters 1-8 match STARLIB;
# REACLIB2 additionally defines 9-11.
const _REACLIB_ARITIES = Dict(
    1 => (1, 1),
    2 => (1, 2),
    3 => (1, 3),
    4 => (2, 1),
    5 => (2, 2),
    6 => (2, 3),
    7 => (2, 4),
    8 => (3, 1),
    9 => (3, 2),
    10 => (4, 2),
    11 => (1, 4),
)

# The standard STARLIB 60-point temperature grid in GK. REACLIB analytic fits
# are evaluated onto this grid so REACLIB-derived tables interpolate exactly
# like STARLIB tables.
const STARLIB_T9_GRID = [
    1.0e-3, 2.0e-3, 3.0e-3, 4.0e-3, 5.0e-3, 6.0e-3, 7.0e-3, 8.0e-3, 9.0e-3, 1.0e-2,
    1.1e-2, 1.2e-2, 1.3e-2, 1.4e-2, 1.5e-2, 1.6e-2, 1.8e-2, 2.0e-2, 2.5e-2, 3.0e-2,
    4.0e-2, 5.0e-2, 6.0e-2, 7.0e-2, 8.0e-2, 9.0e-2, 1.0e-1, 1.1e-1, 1.2e-1, 1.3e-1,
    1.4e-1, 1.5e-1, 1.6e-1, 1.8e-1, 2.0e-1, 2.5e-1, 3.0e-1, 3.5e-1, 4.0e-1, 4.5e-1,
    5.0e-1, 6.0e-1, 7.0e-1, 8.0e-1, 9.0e-1, 1.0e0, 1.25e0, 1.5e0, 1.75e0, 2.0e0,
    2.5e0, 3.0e0, 3.5e0, 4.0e0, 5.0e0, 6.0e0, 7.0e0, 8.0e0, 9.0e0, 1.0e1,
]

# Cap on the REACLIB fit exponent so pathological extrapolation outside a fit's
# validity range yields a large finite rate instead of Inf.
const REACLIB_EXPONENT_CAP = 500.0

function _parse_reaclib_float(field::AbstractString)
    s = strip(field)
    isempty(s) && return 0.0
    return parse(Float64, s)
end

"""
    read_reaclib(path=DEFAULT_REACLIB_PATH)

Read a JINA REACLIB library (Reaclib1 or Reaclib2 layout) into `ReaclibSet`
entries. Every fit set is kept, including alternate set labels for the same
reaction and detailed-balance reverse fits; use `reaclib_rate_tables` or
`iliadis2002_rate_tables` to select and evaluate them.
"""
function read_reaclib(path::AbstractString=_default_reaclib_path())
    sets = ReaclibSet[]
    open(path, "r") do io
        chapter = 0
        line_number = 0
        while !eof(io)
            line = readline(io)
            line_number += 1
            stripped = strip(line)
            isempty(stripped) && continue

            if occursin(r"^\d+$", stripped)
                chapter = parse(Int, stripped)
                haskey(_REACLIB_ARITIES, chapter) || error("Unsupported REACLIB chapter $chapter at line $line_number of $path")
                continue
            end

            chapter >= 1 || error("REACLIB set header before any chapter marker at line $line_number of $path")
            header = rpad(line, 74)
            raw_species = [strip(header[(6 + 5 * i):(10 + 5 * i)]) for i in 0:5]
            label = strip(header[44:47])
            resonance = header[48]
            reverse = header[49] == 'v'
            q_value = _parse_reaclib_float(header[53:64])

            eof(io) && error("Unexpected end of file after REACLIB header at line $line_number of $path")
            first_coefficients = rpad(readline(io), 52)
            line_number += 1
            eof(io) && error("Unexpected end of file inside REACLIB set at line $line_number of $path")
            second_coefficients = rpad(readline(io), 39)
            line_number += 1

            coefficients = (
                _parse_reaclib_float(first_coefficients[1:13]),
                _parse_reaclib_float(first_coefficients[14:26]),
                _parse_reaclib_float(first_coefficients[27:39]),
                _parse_reaclib_float(first_coefficients[40:52]),
                _parse_reaclib_float(second_coefficients[1:13]),
                _parse_reaclib_float(second_coefficients[14:26]),
                _parse_reaclib_float(second_coefficients[27:39]),
            )

            nreactants, nproducts = _REACLIB_ARITIES[chapter]
            species = [normalize_species_name(s) for s in raw_species if !isempty(s)]
            length(species) == nreactants + nproducts ||
                error("REACLIB chapter $chapter set at line $(line_number - 2) of $path has $(length(species)) species, expected $(nreactants + nproducts)")

            push!(sets, ReaclibSet(
                chapter,
                species[1:nreactants],
                species[(nreactants + 1):end],
                label,
                resonance,
                reverse,
                q_value,
                coefficients,
            ))
        end
    end
    return sets
end

"""
    reaclib_rate(set, T9)
    reaclib_rate(sets, T9)

Evaluate the REACLIB analytic rate parameterization

    rate = exp(a1 + a2/T9 + a3*T9^(-1/3) + a4*T9^(1/3) + a5*T9 + a6*T9^(5/3) + a7*ln(T9))

for one fit set, or the sum over a vector of sets belonging to one reaction.
"""
function reaclib_rate(set::ReaclibSet, T9::Real)
    T = Float64(T9)
    T > 0.0 || throw(ArgumentError("T9 must be positive"))
    T13 = cbrt(T)
    a = set.coefficients
    exponent = a[1] + a[2] / T + a[3] / T13 + a[4] * T13 + a[5] * T + a[6] * T * T13 * T13 + a[7] * log(T)
    return exp(min(exponent, REACLIB_EXPONENT_CAP))
end

function reaclib_rate(sets::AbstractVector{ReaclibSet}, T9::Real)
    isempty(sets) && throw(ArgumentError("cannot evaluate an empty set of REACLIB fits"))
    return sum(reaclib_rate(set, T9) for set in sets)
end

_reaclib_group_key(set::ReaclibSet) = (set.chapter, Tuple(set.reactants), Tuple(set.products), lowercase(set.label), set.reverse)

#=
Build one `ReactionRateTable` from the fit sets of a single reaction/label
group. The source string follows the STARLIB convention of label plus flag
suffixes, e.g. `"nacr"`, `"wc12w"`, `"nacrv"`, so weak-rate and source
handling downstream behaves identically for STARLIB and REACLIB tables.
REACLIB carries no rate uncertainties, so the factor uncertainty is 1.
=#
function _reaclib_group_table(sets::AbstractVector{ReaclibSet}, T9_grid::AbstractVector{<:Real})
    representative = first(sets)
    grid = Float64.(collect(T9_grid))
    rate = Float64[max(reaclib_rate(sets, T9), LOG_INTERPOLATION_FLOOR) for T9 in grid]
    weak = any(set -> set.resonance == 'w', sets)
    source = representative.label * (weak ? "w" : "") * (representative.reverse ? "v" : "")

    return ReactionRateTable(
        representative.chapter,
        copy(representative.reactants),
        copy(representative.products),
        source,
        representative.q_value,
        grid,
        rate,
        ones(Float64, length(grid)),
    )
end

_reaclib_set_fingerprint(set::ReaclibSet) =
    (_reaclib_group_key(set), set.resonance, set.q_value, set.coefficients)

#=
Drop exact duplicate fit sets while keeping file order. Needed because the
label-specific REACLIB downloads (`data/reaclib_il01.dat`,
`data/reaclib_nacr.dat`) overlap with the sets a library snapshot already
contains; summing a duplicated set would double the rate.
=#
function _dedup_reaclib_sets(sets::AbstractVector{ReaclibSet})
    seen = Set{Any}()
    unique_sets = ReaclibSet[]
    for set in sets
        fingerprint = _reaclib_set_fingerprint(set)
        fingerprint in seen && continue
        push!(seen, fingerprint)
        push!(unique_sets, set)
    end
    return unique_sets
end

function _group_reaclib_sets(sets::AbstractVector{ReaclibSet}, keep::Function)
    order = Any[]
    groups = Dict{Any,Vector{ReaclibSet}}()
    for set in _dedup_reaclib_sets(sets)
        keep(set) || continue
        key = _reaclib_group_key(set)
        if !haskey(groups, key)
            groups[key] = ReaclibSet[]
            push!(order, key)
        end
        push!(groups[key], set)
    end
    return order, groups
end

"""
    reaclib_rate_tables(sets; T9_grid=STARLIB_T9_GRID, labels=nothing, include_reverse=false)

Evaluate REACLIB fit sets into `ReactionRateTable` entries on a temperature
grid. Sets are grouped by reaction and set label, and all sets of a group
(non-resonant plus resonant terms) are summed. `labels` restricts the output
to specific set labels, e.g. `["nacr", "il01"]`. Detailed-balance reverse
fits are excluded unless `include_reverse=true`.
"""
function reaclib_rate_tables(
    sets::AbstractVector{ReaclibSet};
    T9_grid::AbstractVector{<:Real}=STARLIB_T9_GRID,
    labels=nothing,
    include_reverse::Bool=false,
)
    wanted = labels === nothing ? nothing : Set(lowercase(strip(String(label))) for label in labels)
    keep(set) = (include_reverse || !set.reverse) && (wanted === nothing || lowercase(set.label) in wanted)
    order, groups = _group_reaclib_sets(sets, keep)
    return [_reaclib_group_table(groups[key], T9_grid) for key in order]
end

_reaction_display(reactants::AbstractVector{String}, products::AbstractVector{String}) =
    join(reactants, "+") * " -> " * join(products, "+")

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
        paper = ReactionRateTable[]
        isfile(DEFAULT_ILIADIS2001_PATH) && append!(paper, read_iliadis2001_rates(; partition_functions=pf))
        isfile(DEFAULT_NACRE_PATH) && append!(paper, read_nacre_rates(; partition_functions=pf))
        return iliadis2002_rate_tables(sets; paper_tables=paper, kwargs...)
    end
    return iliadis2002_rate_tables(sets; kwargs...)
end

"""
    find_rate(tables, label; source=nothing)

Find STARLIB rate tables matching a reaction label like `"18F(p,α)15O"`.
Returns all matching tables because STARLIB can contain multiple sources.
"""
function find_rate(tables::AbstractVector{ReactionRateTable}, label::AbstractString; source=nothing)
    reactants, products = parse_reaction_label(label)
    matches = filter(t -> t.reactants == reactants && t.products == products, tables)

    source === nothing && return matches
    wanted = lowercase(strip(string(source)))
    return filter(t -> lowercase(t.source) == wanted, matches)
end

"""
    find_reverse_rate(tables, label; source=nothing)

Find STARLIB tables whose parsed reactants/products are the exact reverse of a
reaction label. This detects explicit reverse rates already present in STARLIB;
it does not synthesize reciprocal-rule reverse rates.
"""
function find_reverse_rate(tables::AbstractVector{ReactionRateTable}, label::AbstractString; source=nothing)
    reactants, products = parse_reaction_label(label)
    matches = filter(t -> t.reactants == products && t.products == reactants, tables)

    source === nothing && return matches
    wanted = lowercase(strip(string(source)))
    return filter(t -> lowercase(t.source) == wanted, matches)
end

function _reaction_participant_key(table::ReactionRateTable)
    return (Tuple(table.reactants), Tuple(table.products))
end

function _unique_reaction_tables(tables::AbstractVector{ReactionRateTable})
    seen = Set{Any}()
    selected = ReactionRateTable[]
    for table in tables
        key = _reaction_participant_key(table)
        key in seen && continue
        push!(selected, table)
        push!(seen, key)
    end
    return selected
end

function _reverse_table_lookup(tables::AbstractVector{ReactionRateTable})
    lookup = Dict{Any,ReactionRateTable}()
    for table in tables
        _valid_nuclear_bookkeeping(table) || continue
        key = _reaction_participant_key(table)
        haskey(lookup, key) || (lookup[key] = table)
    end
    return lookup
end

function _mass_number_factor(names::AbstractVector{String})
    factor = 1.0
    for name in names
        factor *= species_from_name(name).A
    end
    return factor
end

function _can_generate_detailed_balance_reverse(table::ReactionRateTable)
    return length(table.reactants) == 2 &&
           length(table.products) == 1 &&
           table.q_value > 0.0
end

const DETAILED_BALANCE_CAPTURE_CONSTANT = 9.86851e9
const DETAILED_BALANCE_Q_FACTOR = 11.605
const GENERATED_REVERSE_RATE_FLOOR = 1.0e-300
const LOG_INTERPOLATION_FLOOR = 1.0e-300

"""
    PartitionFunctionTable

Ground-state statistical weights `g0 = 2J+1`, temperature-dependent
normalized partition functions `G(T9)`, and mass excesses (MeV) from a JINA
`winvne` file, keyed by normalized species name. `G` is tabulated on the
standard 24-point winvn temperature grid and is defined so that the total
partition function is `g0 * G(T9)`.
"""
struct PartitionFunctionTable
    T9::Vector{Float64}
    g0::Dict{String,Float64}
    G::Dict{String,Vector{Float64}}
    mass_excess::Dict{String,Float64}
end

const DEFAULT_WINVNE_PATH = joinpath(dirname(@__DIR__), "data", "winvne_v2.0.dat")

const _WINVNE_T9_GRID = [
    0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5,
    2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
]

"""
    read_winvne(path=DEFAULT_WINVNE_PATH)

Read a JINA `winvne` nuclear-data file into a `PartitionFunctionTable`.

The file layout is: a title line, a packed temperature-grid line, a species
directory, then one record per nuclide consisting of a header line
`name A Z N spin mass_excess label` followed by 24 partition-function values.
Records are recognized by their header shape, so the directory needs no
special handling. If a normalized name appears more than once (for example
`al-6` normalizing onto `al26`), the first record wins.
"""
function read_winvne(path::AbstractString=DEFAULT_WINVNE_PATH)
    isfile(path) || error("winvne file not found at $path; run data/download_rates.sh to fetch it")

    g0 = Dict{String,Float64}()
    G = Dict{String,Vector{Float64}}()
    mass_excess = Dict{String,Float64}()

    open(path, "r") do io
        line_number = 0
        while !eof(io)
            line = readline(io)
            line_number += 1
            fields = split(line)
            length(fields) >= 6 || continue

            spin = tryparse(Float64, fields[5])
            mass_number = tryparse(Float64, fields[2])
            excess = tryparse(Float64, fields[6])
            (spin === nothing || mass_number === nothing) && continue

            name = normalize_species_name(fields[1])
            values = Float64[]
            sizehint!(values, length(_WINVNE_T9_GRID))
            while length(values) < length(_WINVNE_T9_GRID)
                eof(io) && error("unexpected end of winvne file inside the record of '$name' at line $line_number of $path")
                row = split(readline(io))
                line_number += 1
                for token in row
                    push!(values, parse(Float64, token))
                end
            end
            length(values) == length(_WINVNE_T9_GRID) ||
                error("winvne record of '$name' at line $line_number of $path has $(length(values)) partition-function values, expected $(length(_WINVNE_T9_GRID))")

            if !haskey(g0, name)
                g0[name] = 2.0 * spin + 1.0
                G[name] = values
                excess === nothing || (mass_excess[name] = excess)
            end
        end
    end

    return PartitionFunctionTable(copy(_WINVNE_T9_GRID), g0, G, mass_excess)
end

"""
    reaction_q_value(pf, reactants, products)

Q-value in MeV from winvne mass excesses: `Q = sum(ME reactants) - sum(ME
products)`. Returns `nothing` when any participant is missing from the table.
"""
function reaction_q_value(pf::PartitionFunctionTable, reactants, products)
    total = 0.0
    for name in reactants
        excess = get(pf.mass_excess, normalize_species_name(name), nothing)
        excess === nothing && return nothing
        total += excess
    end
    for name in products
        excess = get(pf.mass_excess, normalize_species_name(name), nothing)
        excess === nothing && return nothing
        total -= excess
    end
    return total
end

"""
    read_tabulated_rates(path; source, partition_functions=nothing,
                         include_total_variants=false)

Read a ReacNetJl tabulated-rate file (as produced by the paper-table
extraction scripts) into `ReactionRateTable` entries. The format is blocks of

    reaction: p ne20 -> na21 ; label: 20Ne(p,g) ; table: 3 ; variant: standard
    <T9> <rate>
    ...
    end

Q-values are computed from winvne mass excesses when `partition_functions`
is given (0 otherwise, which disables generated detailed-balance reverses
for those tables). Rates with `variant: total` duplicate the sum of their
ground/isomer siblings and are skipped unless `include_total_variants=true`.

Rows may carry two columns (`T9 rate`) or four (`T9 rate lower upper`); the
limits become STARLIB-style factor uncertainties `sqrt(upper/lower)`.
Threshold reactions are tabulated only above their onset temperature, where
the papers' omitted rows mean "negligible": each table is padded down to the
STARLIB grid start with the rate floor (and flat up to 10 GK if needed) so
trajectory evaluation never leaves the interpolation domain.
"""
function read_tabulated_rates(
    path::AbstractString;
    source::AbstractString,
    partition_functions=nothing,
    include_total_variants::Bool=false,
)
    isfile(path) || error("tabulated rate file not found at $path")
    tables = ReactionRateTable[]

    reactants = String[]
    products = String[]
    variant = "standard"
    T9 = Float64[]
    rate = Float64[]
    factor_uncertainty = Float64[]

    header_pattern = r"^reaction:\s+(.+?)\s+->\s+(.+?)\s+;\s+label:.*;\s+variant:\s+(\w+)\s*$"

    for (line_number, raw_line) in enumerate(eachline(path))
        line = strip(raw_line)
        (isempty(line) || startswith(line, "#")) && continue

        if startswith(line, "reaction:")
            m = match(header_pattern, line)
            m === nothing && error("malformed reaction header at line $line_number of $path: $line")
            reactants = normalize_species_name.(split(m.captures[1]))
            products = normalize_species_name.(split(m.captures[2]))
            variant = m.captures[3]
            T9 = Float64[]
            rate = Float64[]
            factor_uncertainty = Float64[]
        elseif line == "end"
            isempty(reactants) && error("'end' without a reaction header at line $line_number of $path")
            if !(variant == "total" && !include_total_variants) && !isempty(T9)
                q = 0.0
                if partition_functions !== nothing
                    computed = reaction_q_value(partition_functions, reactants, products)
                    computed === nothing || (q = computed)
                end
                if first(T9) > first(STARLIB_T9_GRID)
                    pushfirst!(T9, first(STARLIB_T9_GRID))
                    pushfirst!(rate, GENERATED_REVERSE_RATE_FLOOR)
                    pushfirst!(factor_uncertainty, 1.0)
                end
                if last(T9) < last(STARLIB_T9_GRID)
                    push!(T9, last(STARLIB_T9_GRID))
                    push!(rate, last(rate))
                    push!(factor_uncertainty, last(factor_uncertainty))
                end
                chapter = length(reactants) == 2 ? (length(products) == 1 ? 4 : 5) : 1
                push!(tables, ReactionRateTable(
                    chapter,
                    copy(reactants),
                    copy(products),
                    String(source),
                    q,
                    copy(T9),
                    copy(rate),
                    copy(factor_uncertainty),
                ))
            end
            reactants = String[]
            products = String[]
        else
            fields = split(line)
            length(fields) in (2, 4) || error("malformed data row at line $line_number of $path: $line")
            push!(T9, parse(Float64, fields[1]))
            push!(rate, parse(Float64, fields[2]))
            if length(fields) == 4
                lower = parse(Float64, fields[3])
                upper = parse(Float64, fields[4])
                (lower > 0.0 && upper >= lower) || error("invalid rate limits at line $line_number of $path: $line")
                push!(factor_uncertainty, sqrt(upper / lower))
            else
                push!(factor_uncertainty, 1.0)
            end
        end
    end

    return tables
end

"""
    read_iliadis2001_rates(path=DEFAULT_ILIADIS2001_PATH; kwargs...)

Read the recommended reaction rates of Iliadis et al. (2001, ApJS 134, 151),
extracted from the paper's Tables 3-9, as `ReactionRateTable`s with source
`"il01tab"`. See `read_tabulated_rates` for the keyword arguments.
"""
read_iliadis2001_rates(path::AbstractString=DEFAULT_ILIADIS2001_PATH; kwargs...) =
    read_tabulated_rates(path; source="il01tab", kwargs...)

"""
    read_nacre_rates(path=DEFAULT_NACRE_PATH; kwargs...)

Read the adopted reaction rates of the NACRE compilation (Angulo et al.
1999), extracted from the paper's rate tables, as `ReactionRateTable`s with
source `"nacrtab"`. See `read_tabulated_rates` for the keyword arguments.
"""
read_nacre_rates(path::AbstractString=DEFAULT_NACRE_PATH; kwargs...) =
    read_tabulated_rates(path; source="nacrtab", kwargs...)

#=
Interpolate the normalized partition function of one species. Outside the
tabulated range the edge values are used; at nova temperatures (below the
0.1 GK grid start) G is 1 to excellent accuracy for all relevant nuclides.
Returns `nothing` when the species is not in the table.
=#
function _partition_function_at(pf::PartitionFunctionTable, name::AbstractString, T9::Real)
    key = normalize_species_name(name)
    haskey(pf.G, key) || return nothing
    values = pf.G[key]
    T = Float64(T9)
    T <= first(pf.T9) && return values[1]
    T >= last(pf.T9) && return values[end]
    return _interpolate_loglog(pf.T9, values, T; value_name="partition function")
end

_statistical_weight(pf::PartitionFunctionTable, name::AbstractString) =
    get(pf.g0, normalize_species_name(name), nothing)

"""
    generated_detailed_balance_reverse_table(table; partition_functions=nothing)

Generate a detailed-balance reverse table for a radiative-capture style
two-body forward reaction `a + b -> c`:

    lambda_rev = C * T9^(3/2) * (Aa*Ab/Ac)^(3/2) * (ga*gb/gc) * (Ga*Gb/Gc)(T9)
                 * exp(-11.605*Q/T9) * N_A<sigma v>_fwd

With `partition_functions` (a `PartitionFunctionTable` from `read_winvne`)
the spin factors `g = 2J+1` and normalized partition functions `G(T9)` are
included; participants missing from the table fall back to `g*G = 1`.
Without it only the mass-factor and Boltzmann terms are used, which is the
historical approximate behavior. Prefer explicit reverse tables from the
rate library when they exist.
"""
function generated_detailed_balance_reverse_table(
    table::ReactionRateTable;
    source::AbstractString="detail_balance",
    partition_functions=nothing,
)
    _can_generate_detailed_balance_reverse(table) || throw(ArgumentError("can only generate detailed-balance reverse rates for exothermic two-body to one-body reactions"))

    product_mass = _mass_number_factor(table.products)
    reactant_mass = _mass_number_factor(table.reactants)
    mass_factor = (reactant_mass / product_mass)^(3.0 / 2.0)

    spin_factor = 1.0
    if partition_functions !== nothing
        weights = [_statistical_weight(partition_functions, name) for name in vcat(table.reactants, table.products)]
        if !any(isnothing, weights)
            spin_factor = (weights[1] * weights[2]) / weights[3]
        end
    end

    reverse_rate = Float64[]
    sizehint!(reverse_rate, length(table.T9))
    for (T9, forward_rate) in zip(table.T9, table.rate)
        boltzmann = exp(-DETAILED_BALANCE_Q_FACTOR * table.q_value / T9)
        pf_ratio = 1.0
        if partition_functions !== nothing
            G_a = _partition_function_at(partition_functions, table.reactants[1], T9)
            G_b = _partition_function_at(partition_functions, table.reactants[2], T9)
            G_c = _partition_function_at(partition_functions, table.products[1], T9)
            if G_a !== nothing && G_b !== nothing && G_c !== nothing
                pf_ratio = (G_a * G_b) / G_c
            end
        end
        value = DETAILED_BALANCE_CAPTURE_CONSTANT * T9^(3.0 / 2.0) * mass_factor * spin_factor * pf_ratio * forward_rate * boltzmann
        push!(reverse_rate, max(value, GENERATED_REVERSE_RATE_FLOOR))
    end

    return ReactionRateTable(
        1,
        copy(table.products),
        copy(table.reactants),
        string(source),
        -table.q_value,
        copy(table.T9),
        reverse_rate,
        copy(table.factor_uncertainty),
    )
end

"""
    add_reverse_reaction_tables(all_tables, forward_tables;
                                generate_detailed_balance=true,
                                partition_functions=nothing)

Return a named tuple containing forward tables plus reverse tables. Exact
reverse entries from the rate library are preferred. If no explicit reverse
exists and the forward table is a supported radiative-capture-style `2 -> 1`
reaction, a detailed-balance reverse table is generated; pass
`partition_functions` from `read_winvne` to include spin factors and
partition-function ratios in the generated rates.
"""
function add_reverse_reaction_tables(
    all_tables::AbstractVector{ReactionRateTable},
    forward_tables::AbstractVector{ReactionRateTable};
    generate_detailed_balance::Bool=true,
    partition_functions=nothing,
)
    reverse_lookup = _reverse_table_lookup(all_tables)
    selected = _unique_reaction_tables(forward_tables)
    seen = Set{Any}(_reaction_participant_key(table) for table in selected)
    explicit = 0
    generated = 0
    missing = 0

    for table in copy(selected)
        reverse_key = (Tuple(table.products), Tuple(table.reactants))
        if haskey(reverse_lookup, reverse_key)
            reverse_table = reverse_lookup[reverse_key]
            key = _reaction_participant_key(reverse_table)
            if !(key in seen)
                push!(selected, reverse_table)
                push!(seen, key)
            end
            explicit += 1
        elseif generate_detailed_balance && _can_generate_detailed_balance_reverse(table)
            reverse_table = generated_detailed_balance_reverse_table(table; partition_functions=partition_functions)
            key = _reaction_participant_key(reverse_table)
            if !(key in seen)
                push!(selected, reverse_table)
                push!(seen, key)
            end
            generated += 1
        else
            missing += 1
        end
    end

    return (tables=selected, explicit=explicit, generated=generated, missing=missing)
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

function network_from_tables(tables::AbstractVector{ReactionRateTable}; species=nothing)
    reactions = [Reaction(table) for table in _unique_reaction_tables(tables)]
    network_species = species === nothing ? _infer_species_from_reactions(reactions) : collect(species)
    return ReactionNetwork(network_species, reactions)
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

function _decay_generator_matrix(network::ReactionNetwork, T9::Real)
    n = length(network.species)
    generator = zeros(Float64, n, n)
    for (reaction, compiled) in zip(network.reactions, network.compiled_reactions)
        length(compiled.reactant_species_indices) == 1 || continue
        source_index = only(compiled.reactant_species_indices)
        rate = interpolate_rate(reaction.rate_table, T9)
        for row in 1:n
            delta = compiled.stoichiometric_delta[row]
            delta == 0.0 && continue
            generator[row, source_index] += delta * rate
        end
    end
    return generator
end

function decay_mass_fractions(
    tables::AbstractVector{ReactionRateTable},
    X::AbstractDict,
    decay_time_s::Real;
    T9::Real=0.1,
)
    decay_time = Float64(decay_time_s)
    decay_time >= 0.0 || throw(ArgumentError("decay_time_s must be non-negative"))

    normalized_X = Dict(normalize_species_name(string(name)) => Float64(value) for (name, value) in X)
    decay_time == 0.0 && return (
        mass_fractions=copy(normalized_X),
        decay_tables=ReactionRateTable[],
        network=nothing,
        times=Float64[0.0],
    )

    decay_tables = select_decay_reaction_tables(tables, keys(normalized_X))
    isempty(decay_tables) && return (
        mass_fractions=copy(normalized_X),
        decay_tables=decay_tables,
        network=nothing,
        times=Float64[0.0, decay_time],
    )

    species = sort!(collect(union(Set(keys(normalized_X)), Set(_infer_species_from_reactions(Reaction.(decay_tables))))); by=name -> begin
        info = species_from_name(name)
        (info.Z, info.A, name)
    end)
    network = network_from_tables(decay_tables; species=species)
    X0 = Dict(name => get(normalized_X, name, 0.0) for name in network.species)
    Y0 = abundances_from_mass_fractions(network, X0)

    generator = _decay_generator_matrix(network, T9)
    Y_final = exp(decay_time * generator) * Y0
    for i in eachindex(Y_final)
        if Y_final[i] < 0.0 && abs(Y_final[i]) <= 1.0e-30
            Y_final[i] = 0.0
        end
    end
    decayed_X = mass_fractions_from_abundances(network, Y_final)
    for (name, value) in decayed_X
        normalized_X[name] = value
    end

    return (
        mass_fractions=normalized_X,
        decay_tables=decay_tables,
        network=network,
        times=Float64[0.0, decay_time],
    )
end

#=
    network_from_labels(tables, labels; species=nothing, source=nothing, on_multiple=:error)

Build a `ReactionNetwork` directly from STARLIB reaction labels. If `species` is
not supplied, the network species list is inferred from the selected reactions.
=#
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

#=
    interpolate_rate(table, T9)

Linearly interpolate a reaction rate at temperature `T9` in GK.

Interpolation is linear in `log(rate)` versus `log(T9)`, which is usually more
reasonable for tabulated thermonuclear rates than linear-linear interpolation.
=#
function _interpolate_loglog(grid::AbstractVector{<:Real}, values::AbstractVector{<:Real}, T9::Real; value_name::AbstractString="value")
    length(grid) == length(values) || throw(ArgumentError("grid and $value_name arrays must have the same length"))
    T = Float64(T9)

    first(grid) <= T <= last(grid) || throw(ArgumentError("T9=$T is outside the table range $(first(grid))–$(last(grid))"))

    i = searchsortedlast(grid, T)
    if i == length(grid) || grid[i] == T
        return Float64(values[i])
    end

    x0 = log(grid[i])
    x1 = log(grid[i+1])
    y0 = log(max(Float64(values[i]), LOG_INTERPOLATION_FLOOR))
    y1 = log(max(Float64(values[i+1]), LOG_INTERPOLATION_FLOOR))
    weight = (log(T) - x0) / (x1 - x0)
    return exp((1 - weight) * y0 + weight * y1)
end

function interpolate_rate(table::ReactionRateTable, T9::Real)
    return _interpolate_loglog(table.T9, table.rate, T9; value_name="rate")
end

"""
    interpolate_factor_uncertainty(table, T9)

Interpolate STARLIB's multiplicative factor uncertainty at temperature `T9`.
The interpolation is log-log, matching `interpolate_rate`.
"""
function interpolate_factor_uncertainty(table::ReactionRateTable, T9::Real)
    return _interpolate_loglog(table.T9, table.factor_uncertainty, T9; value_name="factor uncertainty")
end

"""
    sampled_interpolate_rate(table, T9, p)

Return a STARLIB lognormal sampled reaction rate:

    sampled_rate = recommended_rate * factor_uncertainty^p

where `p` is usually a standard normal deviate held fixed for a reaction during
one Monte Carlo network run.
"""
function sampled_interpolate_rate(table::ReactionRateTable, T9::Real, p::Real)
    rate = interpolate_rate(table, T9)
    factor_uncertainty = interpolate_factor_uncertainty(table, T9)
    return rate * factor_uncertainty^Float64(p)
end

function _factorial(n::Int)
    value = 1
    for i in 2:n
        value *= i
    end
    return value
end

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

function _screening_composition_factor(network::ReactionNetwork, Y::AbstractVector{<:Real})
    factor = 0.0
    for (i, species) in pairs(network.species_info)
        species.Z <= 0 && continue
        factor += (species.Z^2 + species.Z) * max(Float64(Y[i]), 0.0)
    end
    return factor
end

"""
    weak_screening_multiplier(network, reaction, Y, rho, T9)

Return an approximate weak-screening multiplier for charged-particle reactions.
This is a Salpeter-style diagnostic multiplier using the current abundance
composition. Reactions with fewer than two charged reactants return `1.0`.
"""
function weak_screening_multiplier(
    network::ReactionNetwork,
    reaction::Reaction,
    Y::AbstractVector{<:Real},
    rho::Real,
    T9::Real;
    max_exponent::Real=300.0,
)
    length(reaction.reactants) >= 2 || return 1.0
    T6 = 1000.0 * Float64(T9)
    T6 > 0.0 || throw(ArgumentError("T9 must be positive for screening"))
    rho_value = Float64(rho)
    rho_value > 0.0 || throw(ArgumentError("rho must be positive for screening"))

    zeta = _screening_composition_factor(network, Y)
    zeta > 0.0 || return 1.0

    exponent = 0.0
    reactant_info = [species_from_name(name) for name in reaction.reactants]
    for i in 1:(length(reactant_info)-1)
        Zi = reactant_info[i].Z
        Zi <= 0 && continue
        for j in (i+1):length(reactant_info)
            Zj = reactant_info[j].Z
            Zj <= 0 && continue
            exponent += 0.188 * Zi * Zj * sqrt(rho_value * zeta / T6^3)
        end
    end

    exponent <= 0.0 && return 1.0
    return exp(min(exponent, Float64(max_exponent)))
end

# CGS constants for the Chugunov screening evaluation (CODATA 2018).
const _ELEMENTARY_CHARGE_ESU = 4.80320471257e-10
const _BOLTZMANN_ERG_PER_K = 1.380649e-16
const _HBAR_ERG_S = 1.054571817e-27
const _ATOMIC_MASS_UNIT_G = 1.66053906660e-24

# Half-cosine transition between y=x and y=limit, starting at x=start.
# Ported from pynucastro's screening module.
function _smooth_clip(x::Float64, limit::Float64, start::Float64)
    lower, upper = limit < start ? (limit, x) : (x, limit)
    x < min(limit, start) && return lower
    x > max(limit, start) && return upper
    fraction = (1.0 - cos(pi * (x - min(limit, start)) / (start - limit))) / 2.0
    return (1.0 - fraction) * lower + fraction * upper
end

#=
Composition-dependent plasma quantities for Chugunov screening: the electron
number density and the temperature-independent part of the electron Coulomb
coupling parameter. `n_e = rho * sum(Z_i Y_i) / m_u`.
=#
function _ion_plasma_state(network::ReactionNetwork, Y::AbstractVector{<:Real}, rho::Real)
    charge_abundance = 0.0
    for (i, species) in pairs(network.species_info)
        species.Z <= 0 && continue
        charge_abundance += species.Z * max(Float64(Y[i]), 0.0)
    end
    charge_abundance > 0.0 || return nothing

    n_e = Float64(rho) * charge_abundance / _ATOMIC_MASS_UNIT_G
    gamma_e_fac = _ELEMENTARY_CHARGE_ESU^2 / _BOLTZMANN_ERG_PER_K * cbrt(4.0 * pi / 3.0) * cbrt(n_e)
    return (n_e=n_e, gamma_e_fac=gamma_e_fac)
end

#=
Screening exponent h = ln(f_screen) of one ion pair following Chugunov,
DeWitt & Yakovlev (2007), extended to multi-component plasmas as in Yakovlev
et al. (2006). Ported from pynucastro's `chugunov_2007`, which documents the
substitutions (Z^2 -> Z1*Z2, n_i -> n_e/ztilde^3, m_i -> 2*mu12*m_u). Valid
from the weak-screening limit through the strong-screening regime; the fit
caps at Gamma ~ 600 and T ~ 0.1 T_p.
=#
function _chugunov_pair_exponent(Z1::Float64, A1::Float64, Z2::Float64, A2::Float64, temperature_K::Float64, n_e::Float64, gamma_e_fac::Float64)
    ztilde = (cbrt(Z1) + cbrt(Z2)) / 2.0
    reduced_mass = A1 * A2 / (A1 + A2)
    n_i = n_e / ztilde^3
    m_i = 2.0 * reduced_mass * _ATOMIC_MASS_UNIT_G

    T_p = _HBAR_ERG_S / _BOLTZMANN_ERG_PER_K * _ELEMENTARY_CHARGE_ESU * sqrt(4.0 * pi * Z1 * Z2 * n_i / m_i)
    T_norm = _smooth_clip(temperature_K / T_p, 0.1, 0.2)
    Gamma = _smooth_clip(gamma_e_fac * Z1 * Z2 / (ztilde * T_norm * T_p), 600.0, 590.0)

    zeta = cbrt(4.0 / (3.0 * pi^2 * T_norm^2))
    poly = 1.0 + zeta * (0.022 + zeta * ((0.41 - 0.6 / Gamma) + (0.06 + 2.2 / Gamma) * zeta))
    gamtilde = Gamma / cbrt(poly)
    gamtilde2 = gamtilde^2

    A1_fit = 2.7822
    A2_fit = 98.34
    A3_fit = sqrt(3.0) - A1_fit / sqrt(A2_fit)
    B1_fit = -1.7476
    B2_fit = 66.07
    B3_fit = 1.12
    B4_fit = 65.0

    h = gamtilde^1.5 * (A1_fit / sqrt(A2_fit + gamtilde) + A3_fit / (1.0 + gamtilde)) +
        B1_fit * gamtilde2 / (B2_fit + gamtilde) +
        B3_fit * gamtilde2 / (B4_fit + gamtilde2)
    return max(h, 0.0)
end

function _chugunov_reaction_multiplier(compiled::CompiledReaction, context, T9::Real)
    (!context.plasma_active || isempty(compiled.screening_pairs)) && return 1.0
    temperature_K = 1.0e9 * Float64(T9)
    exponent = 0.0
    for (Z1, A1, Z2, A2) in compiled.screening_pairs
        exponent += _chugunov_pair_exponent(Z1, A1, Z2, A2, temperature_K, context.n_e, context.gamma_e_fac)
    end
    return exp(min(exponent, _SCREENING_MAX_EXPONENT))
end

"""
    chugunov_screening_multiplier(network, reaction, Y, rho, T9)

Chugunov, DeWitt & Yakovlev (2007) screening multiplier for one reaction at
the current composition. Available network-wide with `screening=:chugunov`.
"""
function chugunov_screening_multiplier(
    network::ReactionNetwork,
    reaction::Reaction,
    Y::AbstractVector{<:Real},
    rho::Real,
    T9::Real,
)
    T9 > 0.0 || throw(ArgumentError("T9 must be positive for screening"))
    rho > 0.0 || throw(ArgumentError("rho must be positive for screening"))
    length(reaction.reactants) >= 2 || return 1.0

    plasma = _ion_plasma_state(network, Y, rho)
    plasma === nothing && return 1.0

    temperature_K = 1.0e9 * Float64(T9)
    exponent = 0.0
    accumulated_Z = 0
    accumulated_A = 0
    for (k, name) in pairs(reaction.reactants)
        species = try
            species_from_name(name)
        catch
            Species(String(name), 0, 0)
        end
        if k > 1 && accumulated_Z > 0 && species.Z > 0
            exponent += _chugunov_pair_exponent(
                Float64(accumulated_Z), Float64(accumulated_A),
                Float64(species.Z), Float64(species.A),
                temperature_K, plasma.n_e, plasma.gamma_e_fac,
            )
        end
        accumulated_Z += species.Z
        accumulated_A += species.A
    end

    return exp(min(exponent, _SCREENING_MAX_EXPONENT))
end

function _screening_multiplier(screening, network::ReactionNetwork, reaction::Reaction, Y::AbstractVector{<:Real}, rho::Real, T9::Real)
    if screening === nothing || screening === false
        return 1.0
    elseif screening == :weak
        return weak_screening_multiplier(network, reaction, Y, rho, T9)
    elseif screening == :chugunov
        return chugunov_screening_multiplier(network, reaction, Y, rho, T9)
    elseif screening isa Function
        return Float64(screening(network, reaction, Y, rho, T9))
    end

    throw(ArgumentError("unsupported screening=$screening; use nothing, :weak, :chugunov, or a function `(network, reaction, Y, rho, T9) -> multiplier`"))
end

"""
    reaction_flux(reaction, Y, species_index, rho, T9; rate_multiplier=1.0)

Calculate the abundance flux for one reaction.

`Y` is the abundance vector. `species_index` maps species names to positions in
`Y`. `rho` is mass density in g cm^-3, and `T9` is temperature in GK.

For a one-body reaction, the flux is `rate * Yᵢ`. For a two-body reaction using
STARLIB's usual `N_A <σv>` rate, the flux is `rho * rate * Yᵢ * Yⱼ`, with the
standard symmetry correction for identical reactants.
"""
function reaction_flux(
    reaction::Reaction,
    Y::AbstractVector{<:Real},
    species_index::AbstractDict,
    rho::Real,
    T9::Real;
    rate_multiplier::Real=1.0,
    rate_p_value=nothing,
)
    nreactants = length(reaction.reactants)
    nreactants >= 1 || throw(ArgumentError("reaction must have at least one reactant"))

    base_rate = rate_p_value === nothing ? interpolate_rate(reaction.rate_table, T9) : sampled_interpolate_rate(reaction.rate_table, T9, rate_p_value)
    rate = rate_multiplier * base_rate
    abundance_product = 1.0
    for name in reaction.reactants
        abundance_product *= Y[_species_index(species_index, name)]
    end

    density_factor = Float64(rho)^(nreactants - 1)
    return density_factor * rate * abundance_product / _symmetry_factor(reaction.reactants)
end

function _reaction_flux(
    network::ReactionNetwork,
    reaction::Reaction,
    compiled::CompiledReaction,
    Y::AbstractVector{<:Real},
    rho::Real,
    T9::Real;
    rate_multiplier::Real=1.0,
    rate_p_value=nothing,
    screening=nothing,
)
    compiled.nreactants >= 1 || throw(ArgumentError("reaction must have at least one reactant"))

    base_rate = rate_p_value === nothing ? interpolate_rate(reaction.rate_table, T9) : sampled_interpolate_rate(reaction.rate_table, T9, rate_p_value)
    rate = rate_multiplier * _screening_multiplier(screening, network, reaction, Y, rho, T9) * base_rate
    abundance_product = 1.0
    for index in compiled.reactant_indices
        abundance_product *= Y[index]
    end

    density_factor = Float64(rho)^(compiled.nreactants - 1)
    return density_factor * rate * abundance_product / compiled.symmetry_factor
end

#=
    network_rhs(Y, reactions, species_index, rho, T9; rate_multipliers=nothing)

Calculate `dY/dt` for a single-zone reaction network at fixed density and
temperature.

`rate_multipliers` can be supplied as a vector the same length as `reactions`.
This gives us a simple hook for later STARLIB uncertainty studies.
=#
function network_rhs(
    Y::AbstractVector{<:Real},
    reactions::AbstractVector{Reaction},
    species_index::AbstractDict,
    rho::Real,
    T9::Real;
    rate_multipliers=nothing,
    rate_p_values=nothing,
)
    dYdt = zeros(Float64, length(Y))

    if rate_multipliers !== nothing && length(rate_multipliers) != length(reactions)
        throw(ArgumentError("rate_multipliers must have the same length as reactions"))
    end
    if rate_p_values !== nothing && length(rate_p_values) != length(reactions)
        throw(ArgumentError("rate_p_values must have the same length as reactions"))
    end

    for (reaction_number, reaction) in pairs(reactions)
        multiplier = rate_multipliers === nothing ? 1.0 : rate_multipliers[reaction_number]
        p_value = rate_p_values === nothing ? nothing : rate_p_values[reaction_number]
        flux = reaction_flux(reaction, Y, species_index, rho, T9; rate_multiplier=multiplier, rate_p_value=p_value)

        for name in reaction.reactants
            dYdt[_species_index(species_index, name)] -= flux
        end

        for name in reaction.products
            dYdt[_species_index(species_index, name)] += flux
        end
    end

    return dYdt
end

function network_rhs(
    Y::AbstractVector{<:Real},
    network::ReactionNetwork,
    rho::Real,
    T9::Real;
    rate_multipliers=nothing,
    rate_p_values=nothing,
    screening=nothing,
)
    dYdt = zeros(Float64, length(Y))

    length(Y) == length(network.species) || throw(ArgumentError("Y length must match the number of network species"))
    if rate_multipliers !== nothing && length(rate_multipliers) != length(network.reactions)
        throw(ArgumentError("rate_multipliers must have the same length as reactions"))
    end
    if rate_p_values !== nothing && length(rate_p_values) != length(network.reactions)
        throw(ArgumentError("rate_p_values must have the same length as reactions"))
    end

    for (reaction_number, reaction) in pairs(network.reactions)
        compiled = network.compiled_reactions[reaction_number]
        multiplier = rate_multipliers === nothing ? 1.0 : rate_multipliers[reaction_number]
        p_value = rate_p_values === nothing ? nothing : rate_p_values[reaction_number]
        flux = _reaction_flux(network, reaction, compiled, Y, rho, T9; rate_multiplier=multiplier, rate_p_value=p_value, screening=screening)

        for (index, count) in zip(compiled.reactant_species_indices, compiled.reactant_species_counts)
            dYdt[index] -= count * flux
        end
        for (index, count) in zip(compiled.product_species_indices, compiled.product_species_counts)
            dYdt[index] += count * flux
        end
    end

    return dYdt
end

#=
Per-step solver cache.

Within one solver step the temperature and density are fixed, so the
interpolated rate, the density power, the symmetry factor, and any rate
multipliers are constant. `NetworkStepCache` hoists all of that out of the
abundance loops: an RHS or Jacobian evaluation then only multiplies cached
prefactors by abundance products. This changes no mathematics — only where
the temperature-dependent work is done.

Custom screening functions cannot be decomposed this way, so
`_build_step_cache` returns `nothing` for them and callers fall back to the
uncached path.
=#
struct NetworkStepCache
    rho::Float64
    T9::Float64
    prefactors::Vector{Float64}
    screening::Symbol
    screening_scale::Float64
end

function _build_step_cache(
    network::ReactionNetwork,
    rho_t::Real,
    T9_t::Real;
    rate_multipliers=nothing,
    rate_p_values=nothing,
    screening=nothing,
)
    screening isa Function && return nothing
    mode = screening === nothing || screening === false ? :none :
           screening == :weak ? :weak :
           screening == :chugunov ? :chugunov :
           throw(ArgumentError("unsupported screening=$screening; use nothing, :weak, :chugunov, or a function `(network, reaction, Y, rho, T9) -> multiplier`"))

    nreactions = length(network.reactions)
    if rate_multipliers !== nothing && length(rate_multipliers) != nreactions
        throw(ArgumentError("rate_multipliers must have the same length as reactions"))
    end
    if rate_p_values !== nothing && length(rate_p_values) != nreactions
        throw(ArgumentError("rate_p_values must have the same length as reactions"))
    end

    rho_value = Float64(rho_t)
    T9_value = Float64(T9_t)
    prefactors = Vector{Float64}(undef, nreactions)
    for r in 1:nreactions
        table = network.reactions[r].rate_table
        compiled = network.compiled_reactions[r]
        p_value = rate_p_values === nothing ? nothing : rate_p_values[r]
        base_rate = p_value === nothing ? interpolate_rate(table, T9_value) : sampled_interpolate_rate(table, T9_value, p_value)
        multiplier = rate_multipliers === nothing ? 1.0 : Float64(rate_multipliers[r])
        prefactors[r] = multiplier * base_rate * rho_value^(compiled.nreactants - 1) / compiled.symmetry_factor
    end

    screening_scale = 0.0
    if mode != :none
        T9_value > 0.0 || throw(ArgumentError("T9 must be positive for screening"))
        rho_value > 0.0 || throw(ArgumentError("rho must be positive for screening"))
        if mode == :weak
            T6 = 1000.0 * T9_value
            screening_scale = 0.188 * sqrt(rho_value / T6^3)
        end
    end

    return NetworkStepCache(rho_value, T9_value, prefactors, mode, screening_scale)
end

const _SCREENING_MAX_EXPONENT = 300.0

function _screening_zeta_scale(cache::NetworkStepCache, network::ReactionNetwork, Y::AbstractVector{Float64})
    cache.screening == :weak || return 0.0
    zeta = _screening_composition_factor(network, Y)
    zeta > 0.0 || return 0.0
    return cache.screening_scale * sqrt(zeta)
end

@inline function _cached_screening_multiplier(compiled::CompiledReaction, zeta_scale::Float64)
    (zeta_scale > 0.0 && compiled.charge_pair_sum > 0.0) || return 1.0
    return exp(min(compiled.charge_pair_sum * zeta_scale, _SCREENING_MAX_EXPONENT))
end

# Per-evaluation screening context: the composition-dependent scalars shared
# by every reaction of one RHS/Jacobian evaluation. A concrete struct keeps
# the per-reaction screening call dispatch-free and allocation-free.
struct ScreeningContext
    zeta_scale::Float64
    n_e::Float64
    gamma_e_fac::Float64
    plasma_active::Bool
end

function _screening_context(cache::NetworkStepCache, network::ReactionNetwork, Y::AbstractVector{Float64})
    if cache.screening == :weak
        return ScreeningContext(_screening_zeta_scale(cache, network, Y), 0.0, 0.0, false)
    elseif cache.screening == :chugunov
        plasma = _ion_plasma_state(network, Y, cache.rho)
        plasma === nothing && return ScreeningContext(0.0, 0.0, 0.0, false)
        return ScreeningContext(0.0, plasma.n_e, plasma.gamma_e_fac, true)
    end
    return ScreeningContext(0.0, 0.0, 0.0, false)
end

@inline function _cached_reaction_screening(cache::NetworkStepCache, compiled::CompiledReaction, context::ScreeningContext)
    if cache.screening == :weak
        return _cached_screening_multiplier(compiled, context.zeta_scale)
    elseif cache.screening == :chugunov
        return _chugunov_reaction_multiplier(compiled, context, cache.T9)
    end
    return 1.0
end

# In-place dY/dt with all temperature-dependent factors taken from the cache.
function _cached_network_rhs!(dYdt::Vector{Float64}, network::ReactionNetwork, cache::NetworkStepCache, Y::AbstractVector{Float64})
    fill!(dYdt, 0.0)
    screening_context = _screening_context(cache, network, Y)

    @inbounds for r in eachindex(network.compiled_reactions)
        compiled = network.compiled_reactions[r]
        flux = cache.prefactors[r] * _cached_reaction_screening(cache, compiled, screening_context)
        for index in compiled.reactant_indices
            flux *= Y[index]
        end
        flux == 0.0 && continue

        for (index, count) in zip(compiled.reactant_species_indices, compiled.reactant_species_counts)
            dYdt[index] -= count * flux
        end
        for (index, count) in zip(compiled.product_species_indices, compiled.product_species_counts)
            dYdt[index] += count * flux
        end
    end

    return dYdt
end

#=
Analytic Jacobian d(dY/dt)/dY. The flux of each reaction is a polynomial in
the reactant abundances, so the derivative follows from the product rule over
the distinct reactant species. The weak-screening multiplier is treated as
constant with respect to Y; Newton's converged answer is fixed by the exact
residual alone, so this only shapes the iteration path, not the solution.
=#
function _cached_network_jacobian!(J::Matrix{Float64}, network::ReactionNetwork, cache::NetworkStepCache, Y::AbstractVector{Float64})
    fill!(J, 0.0)
    screening_context = _screening_context(cache, network, Y)

    @inbounds for r in eachindex(network.compiled_reactions)
        compiled = network.compiled_reactions[r]
        base = cache.prefactors[r] * _cached_reaction_screening(cache, compiled, screening_context)
        base == 0.0 && continue

        reactant_indices = compiled.reactant_species_indices
        reactant_counts = compiled.reactant_species_counts
        for jpos in eachindex(reactant_indices)
            jindex = reactant_indices[jpos]
            count_j = reactant_counts[jpos]
            derivative = Float64(count_j)
            count_j > 1 && (derivative *= Y[jindex]^(count_j - 1))
            for kpos in eachindex(reactant_indices)
                kpos == jpos && continue
                derivative *= Y[reactant_indices[kpos]]^reactant_counts[kpos]
            end
            dflux = base * derivative
            dflux == 0.0 && continue

            for (index, count) in zip(reactant_indices, reactant_counts)
                J[index, jindex] -= count * dflux
            end
            for (index, count) in zip(compiled.product_species_indices, compiled.product_species_counts)
                J[index, jindex] += count * dflux
            end
        end
    end

    return J
end

#=
    reaction_string(reaction)

Return a compact, normalized reaction label for display or diagnostics.

For common two-body STARLIB reactions, this returns labels like
`"f18(p,he4)o15"` or `"f18(p,γ)ne19"`.
=#
function _display_species_name(name::AbstractString)
    normalized = normalize_species_name(name)
    normalized == "p" && return "p"
    normalized == "n" && return "n"
    normalized == "d" && return "d"
    normalized == "t" && return "t"
    normalized == "he4" && return "α"
    normalized == "he3" && return "3He"

    m = match(r"^([a-z]+)\*(\d+)$", normalized)
    if m !== nothing
        species = species_from_name(normalized)
        return string(species.A, uppercasefirst(m.captures[1]), "*")
    end

    species = species_from_name(normalized)
    symbol_match = match(r"^([a-z]+)\d+$", normalized)
    symbol = symbol_match === nothing ? normalized : symbol_match.captures[1]
    return string(species.A, uppercasefirst(symbol))
end

function _format_species_group(names::AbstractVector{String})
    isempty(names) && return "γ"

    parts = String[]
    seen = Set{String}()
    for name in names
        normalized = normalize_species_name(name)
        normalized in seen && continue
        multiplicity = Base.count(==(normalized), normalize_species_name.(names))
        label = _display_species_name(normalized)
        push!(parts, multiplicity == 1 ? label : string(multiplicity, label))
        push!(seen, normalized)
    end
    return join(parts, "+")
end

function _weak_symbol(delta_Z::Int)
    delta_Z == -1 && return "β+"
    delta_Z == 1 && return "β-"
    return "weak"
end

function reaction_string(reaction::Reaction)
    conservation = reaction_conservation(reaction)

    if conservation.is_weak_decay
        weak = _weak_symbol(conservation.delta_Z)
        if length(reaction.reactants) == 1
            return string(_format_species_group(reaction.reactants), "(", weak, ")", _format_species_group(reaction.products))
        elseif length(reaction.reactants) == 2
            projectile, target = reaction.reactants
            return string(
                _display_species_name(target),
                "(",
                _display_species_name(projectile),
                ",eν)",
                _format_species_group(reaction.products),
            )
        end
        return string(_format_species_group(reaction.reactants), "(", weak, ")->", _format_species_group(reaction.products))
    end

    if length(reaction.reactants) == 2
        projectile, target = reaction.reactants
        if length(reaction.products) == 1
            return string(_display_species_name(target), "(", _display_species_name(projectile), ",γ)", _format_species_group(reaction.products))
        elseif length(reaction.products) >= 2
            ejectiles = reaction.products[1:end-1]
            product = reaction.products[end:end]
            return string(
                _display_species_name(target),
                "(",
                _display_species_name(projectile),
                ",",
                _format_species_group(ejectiles),
                ")",
                _format_species_group(product),
            )
        end
    end

    return _format_species_group(reaction.reactants) * "->" * _format_species_group(reaction.products)
end

#=
    reaction_fluxes(network, Y, rho, T9; rate_multipliers=nothing)

Calculate instantaneous reaction fluxes for every reaction in a network.
The returned vector is ordered like `network.reactions`.
=#
function reaction_fluxes(
    network::ReactionNetwork,
    Y::AbstractVector{<:Real},
    rho::Real,
    T9::Real;
    rate_multipliers=nothing,
    rate_p_values=nothing,
    screening=nothing,
)
    if rate_multipliers !== nothing && length(rate_multipliers) != length(network.reactions)
        throw(ArgumentError("rate_multipliers must have the same length as network.reactions"))
    end
    if rate_p_values !== nothing && length(rate_p_values) != length(network.reactions)
        throw(ArgumentError("rate_p_values must have the same length as network.reactions"))
    end

    fluxes = zeros(Float64, length(network.reactions))
    for (i, reaction) in pairs(network.reactions)
        multiplier = rate_multipliers === nothing ? 1.0 : rate_multipliers[i]
        p_value = rate_p_values === nothing ? nothing : rate_p_values[i]
        fluxes[i] = _reaction_flux(network, reaction, network.compiled_reactions[i], Y, rho, T9; rate_multiplier=multiplier, rate_p_value=p_value, screening=screening)
    end
    return fluxes
end

#=
    reaction_flux_history(network, history, times, rho, T9; rate_multipliers=nothing)

Calculate reaction fluxes at every saved timestep from `solve_network` output.
`rho` and `T9` may be constants or functions of time.

Returns a matrix where `flux_history[n, r]` is the flux of reaction `r` at
`times[n]`.
=#
function reaction_flux_history(
    network::ReactionNetwork,
    history::AbstractMatrix{<:Real},
    times::AbstractVector{<:Real},
    rho,
    T9;
    rate_multipliers=nothing,
    rate_p_values=nothing,
    screening=nothing,
)
    length(times) == size(history, 1) || throw(ArgumentError("times length must match the number of history rows"))
    size(history, 2) == length(network.species) || throw(ArgumentError("history column count must match the number of network species"))

    flux_history = Matrix{Float64}(undef, length(times), length(network.reactions))
    for (n, t) in pairs(times)
        rho_t = _profile_value(rho, t)
        T9_t = _profile_value(T9, t)
        flux_history[n, :] .= reaction_fluxes(network, view(history, n, :), rho_t, T9_t; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)
    end
    return flux_history
end

#=
    integrated_fluxes(times, flux_history)

Integrate reaction flux histories in time using the trapezoid rule.
Returns one integrated flux per reaction.
=#
function integrated_fluxes(times::AbstractVector{<:Real}, flux_history::AbstractMatrix{<:Real})
    length(times) == size(flux_history, 1) || throw(ArgumentError("times length must match the number of flux-history rows"))
    length(times) >= 2 || throw(ArgumentError("at least two time points are required"))

    totals = zeros(Float64, size(flux_history, 2))
    for n in 1:(length(times)-1)
        dt = Float64(times[n+1] - times[n])
        dt >= 0.0 || throw(ArgumentError("times must be monotonically increasing"))
        totals .+= 0.5 * dt .* (view(flux_history, n, :) .+ view(flux_history, n + 1, :))
    end
    return totals
end

#=
    energy_generation_rate(network, Y, rho, T9; ...)

Return diagnostic nuclear energy generation in erg g^-1 s^-1 from reaction
Q-values and instantaneous reaction fluxes. This does not feed back into the
temperature trajectory.
=#
function energy_generation_rate(
    network::ReactionNetwork,
    Y::AbstractVector{<:Real},
    rho::Real,
    T9::Real;
    rate_multipliers=nothing,
    rate_p_values=nothing,
    screening=nothing,
)
    fluxes = reaction_fluxes(network, Y, rho, T9; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)
    epsilon = 0.0
    for (i, reaction) in pairs(network.reactions)
        epsilon += reaction.rate_table.q_value * fluxes[i]
    end
    return epsilon * AVOGADRO * MEV_TO_ERG
end

#=
    energy_generation_history(network, history, times, rho, T9; ...)

Return diagnostic nuclear energy generation in erg g^-1 s^-1 at every saved
history row.
=#
function energy_generation_history(
    network::ReactionNetwork,
    history::AbstractMatrix{<:Real},
    times::AbstractVector{<:Real},
    rho,
    T9;
    rate_multipliers=nothing,
    rate_p_values=nothing,
    screening=nothing,
)
    length(times) == size(history, 1) || throw(ArgumentError("times length must match the number of history rows"))
    size(history, 2) == length(network.species) || throw(ArgumentError("history column count must match the number of network species"))

    epsilon = Vector{Float64}(undef, length(times))
    for (n, t) in pairs(times)
        rho_t = _profile_value(rho, t)
        T9_t = _profile_value(T9, t)
        epsilon[n] = energy_generation_rate(network, view(history, n, :), rho_t, T9_t; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)
    end
    return epsilon
end

#=
    integrated_energy_generation(times, epsilon_history)

Integrate diagnostic energy generation in time using the trapezoid rule.
Returns specific energy release in erg g^-1.
=#
function integrated_energy_generation(times::AbstractVector{<:Real}, epsilon_history::AbstractVector{<:Real})
    length(times) == length(epsilon_history) || throw(ArgumentError("times length must match energy-history length"))
    length(times) >= 2 || throw(ArgumentError("at least two time points are required"))

    total = 0.0
    for n in 1:(length(times)-1)
        dt = Float64(times[n+1] - times[n])
        dt >= 0.0 || throw(ArgumentError("times must be monotonically increasing"))
        total += 0.5 * dt * (epsilon_history[n] + epsilon_history[n+1])
    end
    return total
end

"""
    species_flux_balance(network, Y, rho, T9; rate_multipliers=nothing)

Calculate instantaneous production, destruction, and net `dY/dt` contributions
for every species.

Returns `(production=..., destruction=..., net=...)`, with vectors ordered like
`network.species`.
"""
function species_flux_balance(
    network::ReactionNetwork,
    Y::AbstractVector{<:Real},
    rho::Real,
    T9::Real;
    rate_multipliers=nothing,
    rate_p_values=nothing,
    screening=nothing,
)
    fluxes = reaction_fluxes(network, Y, rho, T9; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)
    production = zeros(Float64, length(network.species))
    destruction = zeros(Float64, length(network.species))

    for (reaction_index, reaction) in pairs(network.reactions)
        compiled = network.compiled_reactions[reaction_index]
        flux = fluxes[reaction_index]

        for (index, count) in zip(compiled.product_species_indices, compiled.product_species_counts)
            production[index] += count * flux
        end

        for (index, count) in zip(compiled.reactant_species_indices, compiled.reactant_species_counts)
            destruction[index] += count * flux
        end
    end

    return (production=production, destruction=destruction, net=production .- destruction)
end

#=
    reaction_edges(network)

Return a flattened graph-like edge list for diagnostics or external plotting.
Each edge connects one reactant species to one product species for a reaction.
This is a graph approximation of the reaction hypergraph.
=#
function reaction_edges(network::ReactionNetwork)
    edges = NamedTuple[]
    for (reaction_index, reaction) in pairs(network.reactions)
        label = reaction_string(reaction)
        for reactant in reaction.reactants
            for product in reaction.products
                push!(edges, (reaction_index=reaction_index, reaction=label, from=reactant, to=product))
            end
        end
    end
    return edges
end

function _species_property_sum(names::AbstractVector{String}, property::Symbol)
    total = 0
    for name in names
        species = species_from_name(name)
        if property == :A
            total += species.A
        elseif property == :Z
            total += species.Z
        else
            throw(ArgumentError("unsupported species property $property"))
        end
    end
    return total
end

"""
    reaction_conservation(reaction)

Check baryon-number (`A`) and charge (`Z`) conservation for one reaction.
Returns a named tuple with reactant/product totals and boolean conservation flags.

Photons, if present, are treated as `A=0`, `Z=0`.
"""
function reaction_conservation(reaction::Reaction)
    reactant_A = _species_property_sum(reaction.reactants, :A)
    product_A = _species_property_sum(reaction.products, :A)
    reactant_Z = _species_property_sum(reaction.reactants, :Z)
    product_Z = _species_property_sum(reaction.products, :Z)
    delta_A = product_A - reactant_A
    delta_Z = product_Z - reactant_Z
    is_weak_decay = delta_A == 0 && abs(delta_Z) == 1

    return (
        reaction=join(reaction.reactants, "+") * "->" * join(reaction.products, "+"),
        reactant_A=reactant_A,
        product_A=product_A,
        delta_A=delta_A,
        conserves_A=reactant_A == product_A,
        reactant_Z=reactant_Z,
        product_Z=product_Z,
        delta_Z=delta_Z,
        conserves_Z=reactant_Z == product_Z,
        is_weak_decay=is_weak_decay,
        valid_nuclear_bookkeeping=(reactant_A == product_A) && ((reactant_Z == product_Z) || is_weak_decay),
    )
end

function _network_species_issues(network::ReactionNetwork)
    issues = String[]

    if length(unique(network.species)) != length(network.species)
        push!(issues, "network species list contains duplicates")
    end

    if length(network.species_info) != length(network.species)
        push!(issues, "network species_info length does not match species length")
    end

    for (i, name) in pairs(network.species)
        haskey(network.species_index, name) || push!(issues, "species '$name' is missing from species_index")
        get(network.species_index, name, 0) == i || push!(issues, "species_index for '$name' does not match species order")

        if i <= length(network.species_info)
            info = network.species_info[i]
            info.name == name || push!(issues, "species_info[$i] name '$(info.name)' does not match species '$name'")
        end
    end

    return issues
end

function _reaction_validation_report(network::ReactionNetwork, reaction_index::Int, reaction::Reaction)
    issues = String[]

    for name in Iterators.flatten((reaction.reactants, reaction.products))
        haskey(network.species_index, name) || push!(issues, "species '$name' is not present in the network")
    end

    conservation = try
        reaction_conservation(reaction)
    catch err
        push!(issues, "could not check conservation: $(err)")
        nothing
    end

    if conservation !== nothing
        conservation.conserves_A || push!(issues, "baryon number is not conserved: ΔA=$(conservation.delta_A)")
        if !conservation.conserves_Z && !conservation.is_weak_decay
            push!(issues, "charge is not conserved: ΔZ=$(conservation.delta_Z)")
        end
    end

    return (
        reaction_index=reaction_index,
        reaction=reaction_string(reaction),
        reactants=copy(reaction.reactants),
        products=copy(reaction.products),
        conservation=conservation,
        valid=isempty(issues),
        issues=issues,
    )
end

"""
    network_validation_report(network; throw_on_error=false)

Validate basic reaction-network bookkeeping and physics constraints.

Checks currently include:
- consistency of `species`, `species_info`, and `species_index`
- all reaction species are present in the network
- baryon-number conservation for every reaction
- charge conservation for every reaction

Returns a named tuple containing summary information, per-reaction reports, and
all issues. If `throw_on_error=true`, invalid networks raise an `ArgumentError`.
"""
function network_validation_report(network::ReactionNetwork; throw_on_error::Bool=false)
    issues = _network_species_issues(network)
    reaction_reports = NamedTuple[]

    for (reaction_index, reaction) in pairs(network.reactions)
        report = _reaction_validation_report(network, reaction_index, reaction)
        push!(reaction_reports, report)
        for issue in report.issues
            push!(issues, "reaction $(reaction_index) ($(reaction_string(reaction))): $issue")
        end
    end

    report = (
        valid=isempty(issues),
        num_species=length(network.species),
        num_reactions=length(network.reactions),
        species=copy(network.species),
        reaction_reports=reaction_reports,
        issues=issues,
    )

    if throw_on_error && !report.valid
        throw(ArgumentError("network validation failed:\n" * join(report.issues, "\n")))
    end

    return report
end

_profile_value(value::Real, t::Real) = Float64(value)
_profile_value(value, t::Real) = Float64(value(t))

function _validate_time_inputs(t_start::Float64, t_end::Float64, dt::Float64)
    t_end > t_start || throw(ArgumentError("tspan must have t_end > t_start"))
    dt > 0 || throw(ArgumentError("dt must be positive"))
end

function _time_grid(tspan::Tuple{<:Real,<:Real}, dt::Real)
    t_start = Float64(tspan[1])
    t_end = Float64(tspan[2])
    step = Float64(dt)
    _validate_time_inputs(t_start, t_end, step)

    times = collect(t_start:step:t_end)
    if isempty(times) || times[end] < t_end
        push!(times, t_end)
    elseif times[end] > t_end
        times[end] = t_end
    end

    return times
end

function _checked_initial_abundances(Y0::AbstractVector{<:Real}, network::ReactionNetwork)
    length(Y0) == length(network.species) || throw(ArgumentError("Y0 length must match the number of network species"))
    return Float64.(Y0)
end

function _rhs_at(network::ReactionNetwork, Y::AbstractVector{<:Real}, rho, T9, t::Real; rate_multipliers=nothing, rate_p_values=nothing, screening=nothing)
    rho_t = _profile_value(rho, t)
    T9_t = _profile_value(T9, t)
    return network_rhs(Y, network, rho_t, T9_t; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)
end

function _step_cache_at(network::ReactionNetwork, rho, T9, t::Float64; rate_multipliers=nothing, rate_p_values=nothing, screening=nothing)
    return _build_step_cache(
        network,
        _profile_value(rho, t),
        _profile_value(T9, t);
        rate_multipliers=rate_multipliers,
        rate_p_values=rate_p_values,
        screening=screening,
    )
end

function _euler_step(network::ReactionNetwork, Y::Vector{Float64}, t::Float64, dt::Float64, rho, T9; rate_multipliers=nothing, rate_p_values=nothing, screening=nothing)
    cache = _step_cache_at(network, rho, T9, t; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)
    if cache === nothing
        k1 = _rhs_at(network, Y, rho, T9, t; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)
    else
        k1 = _cached_network_rhs!(similar(Y), network, cache, Y)
    end
    return Y .+ dt .* k1
end

function _rk4_step(network::ReactionNetwork, Y::Vector{Float64}, t::Float64, dt::Float64, rho, T9; rate_multipliers=nothing, rate_p_values=nothing, screening=nothing)
    cache_start = _step_cache_at(network, rho, T9, t; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)

    if cache_start === nothing
        k1 = _rhs_at(network, Y, rho, T9, t; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)
        k2 = _rhs_at(network, Y .+ 0.5 * dt .* k1, rho, T9, t + 0.5 * dt; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)
        k3 = _rhs_at(network, Y .+ 0.5 * dt .* k2, rho, T9, t + 0.5 * dt; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)
        k4 = _rhs_at(network, Y .+ dt .* k3, rho, T9, t + dt; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)
        return Y .+ (dt / 6.0) .* (k1 .+ 2.0 .* k2 .+ 2.0 .* k3 .+ k4)
    end

    cache_mid = _step_cache_at(network, rho, T9, t + 0.5 * dt; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)
    cache_end = _step_cache_at(network, rho, T9, t + dt; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)
    k1 = _cached_network_rhs!(similar(Y), network, cache_start, Y)
    k2 = _cached_network_rhs!(similar(Y), network, cache_mid, Y .+ 0.5 * dt .* k1)
    k3 = _cached_network_rhs!(similar(Y), network, cache_mid, Y .+ 0.5 * dt .* k2)
    k4 = _cached_network_rhs!(similar(Y), network, cache_end, Y .+ dt .* k3)
    return Y .+ (dt / 6.0) .* (k1 .+ 2.0 .* k2 .+ 2.0 .* k3 .+ k4)
end

function _max_abs(values::AbstractVector{<:Real})
    maximum(abs, values; init=0.0)
end

function _all_finite(values::AbstractVector{<:Real})
    all(isfinite, values)
end

function _positivity_limited_alpha(Y::Vector{Float64}, correction::Vector{Float64})
    alpha = 1.0
    for i in eachindex(Y)
        correction[i] < 0.0 || continue
        Y[i] > 0.0 || continue
        alpha = min(alpha, 0.9 * Y[i] / -correction[i])
    end
    return clamp(alpha, 0.0, 1.0)
end

function _clamp_tiny_negative_trials!(Y::Vector{Float64}, floor::Float64)
    for i in eachindex(Y)
        if Y[i] < 0.0 && abs(Y[i]) <= floor
            Y[i] = 0.0
        end
    end
    return Y
end

struct NewtonConvergenceError <: Exception
    iterations::Int
    residual_norm::Float64
end

function Base.showerror(io::IO, err::NewtonConvergenceError)
    print(io, "backward Euler Newton iteration failed to converge after $(err.iterations) iterations; residual norm = $(err.residual_norm)")
end

function _backward_euler_residual(network::ReactionNetwork, Y_next::Vector{Float64}, Y::Vector{Float64}, t_next::Float64, dt::Float64, rho, T9; rate_multipliers=nothing, rate_p_values=nothing, screening=nothing)
    return Y_next .- Y .- dt .* _rhs_at(network, Y_next, rho, T9, t_next; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)
end

function _backward_euler_jacobian(network::ReactionNetwork, Y_next::Vector{Float64}, residual::Vector{Float64}, Y::Vector{Float64}, t_next::Float64, dt::Float64, rho, T9, finite_difference_epsilon::Float64; rate_multipliers=nothing, rate_p_values=nothing, screening=nothing)
    n = length(Y_next)
    jacobian = Matrix{Float64}(undef, n, n)

    Base.Threads.@threads for j in 1:n
        perturbed = copy(Y_next)
        saved = Y_next[j]
        step = finite_difference_epsilon * max(abs(saved), 1.0)
        perturbed[j] = saved + step
        perturbed_residual = _backward_euler_residual(network, perturbed, Y, t_next, dt, rho, T9; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)
        jacobian[:, j] .= (perturbed_residual .- residual) ./ step
    end

    return jacobian
end

function _backward_euler_step(
    network::ReactionNetwork,
    Y::Vector{Float64},
    t::Float64,
    dt::Float64,
    rho,
    T9;
    rate_multipliers=nothing,
    rate_p_values=nothing,
    screening=nothing,
    newton_tolerance::Real=1.0e-10,
    max_newton_iterations::Integer=20,
    finite_difference_epsilon::Real=sqrt(eps(Float64)),
    jacobian::Symbol=:analytic,
    newton_iterations=nothing,
)
    max_newton_iterations > 0 || throw(ArgumentError("max_newton_iterations must be positive"))
    newton_tolerance > 0.0 || throw(ArgumentError("newton_tolerance must be positive"))
    finite_difference_epsilon > 0.0 || throw(ArgumentError("finite_difference_epsilon must be positive"))
    jacobian in (:analytic, :finite_difference) || throw(ArgumentError("unsupported jacobian=$jacobian; use :analytic or :finite_difference"))

    t_next = t + dt
    Y_next = _euler_step(network, Y, t, dt, rho, T9; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)

    # All Newton residuals of this step share the fixed (rho, T9) at t_next.
    cache = _step_cache_at(network, rho, T9, t_next; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)
    n = length(Y)
    rhs_buffer = cache === nothing ? nothing : Vector{Float64}(undef, n)

    residual_at = function (Y_trial::Vector{Float64})
        if cache === nothing
            return _backward_euler_residual(network, Y_trial, Y, t_next, dt, rho, T9; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)
        end
        _cached_network_rhs!(rhs_buffer, network, cache, Y_trial)
        return Y_trial .- Y .- dt .* rhs_buffer
    end

    residual = residual_at(Y_next)
    tolerance = Float64(newton_tolerance) * max(_max_abs(Y_next), 1.0)
    if _max_abs(residual) <= tolerance
        newton_iterations !== nothing && push!(newton_iterations, 0)
        return Y_next
    end

    use_analytic_jacobian = cache !== nothing && jacobian == :analytic
    jacobian_buffer = use_analytic_jacobian ? Matrix{Float64}(undef, n, n) : nothing

    for iteration in 1:max_newton_iterations
        residual_norm = _max_abs(residual)
        if use_analytic_jacobian
            # Newton matrix of the residual: I - dt * d(dY/dt)/dY.
            _cached_network_jacobian!(jacobian_buffer, network, cache, Y_next)
            @inbounds for j in 1:n
                for i in 1:n
                    jacobian_buffer[i, j] = -dt * jacobian_buffer[i, j]
                end
                jacobian_buffer[j, j] += 1.0
            end
            jacobian_matrix = jacobian_buffer
        else
            jacobian_matrix = _backward_euler_jacobian(
                network,
                Y_next,
                residual,
                Y,
                t_next,
                dt,
                rho,
                T9,
                Float64(finite_difference_epsilon);
                rate_multipliers=rate_multipliers,
                rate_p_values=rate_p_values,
                screening=screening,
            )
        end
        # A singular Newton matrix (e.g. dt so large that I - dt*J loses the
        # identity part to rounding along a conserved direction) is a step
        # failure, not a fatal error: report non-convergence so adaptive
        # drivers shrink dt and retry.
        correction = try
            jacobian_matrix \ (-residual)
        catch err
            if err isa LinearAlgebra.SingularException || err isa LinearAlgebra.LAPACKException
                throw(NewtonConvergenceError(iteration, residual_norm))
            end
            rethrow()
        end

        _all_finite(correction) || throw(NewtonConvergenceError(iteration, residual_norm))

        accepted_trial = false
        best_trial = copy(Y_next)
        best_residual = copy(residual)
        best_residual_norm = residual_norm
        alpha = _positivity_limited_alpha(Y_next, correction)
        alpha = alpha > 0.0 ? alpha : 1.0

        for _ in 1:12
            trial = Y_next .+ alpha .* correction
            _clamp_tiny_negative_trials!(trial, 1.0e-30)
            if _all_finite(trial) && minimum(trial; init=0.0) >= 0.0
                trial_residual = residual_at(trial)
                trial_residual_norm = _max_abs(trial_residual)
                if trial_residual_norm < best_residual_norm
                    best_trial = trial
                    best_residual = trial_residual
                    best_residual_norm = trial_residual_norm
                end
                if trial_residual_norm < residual_norm || trial_residual_norm <= tolerance
                    Y_next = trial
                    residual = trial_residual
                    accepted_trial = true
                    break
                end
            end
            alpha *= 0.5
        end

        if !accepted_trial
            if best_residual_norm < residual_norm
                Y_next = best_trial
                residual = best_residual
            else
                Y_next .+= correction
                residual = residual_at(Y_next)
            end
        end

        tolerance = Float64(newton_tolerance) * max(_max_abs(Y_next), 1.0)
        if _max_abs(residual) <= tolerance || _max_abs(correction) <= tolerance
            newton_iterations !== nothing && push!(newton_iterations, iteration)
            return Y_next
        end
    end

    throw(NewtonConvergenceError(Int(max_newton_iterations), _max_abs(residual)))
end

function _newton_iteration_summary(iterations::Vector{Int}, failed_steps::Integer)
    if isempty(iterations)
        return (
            newton_iterations=Int[],
            mean_newton_iterations=NaN,
            max_newton_iterations=0,
            newton_failed_steps=Int(failed_steps),
        )
    end

    return (
        newton_iterations=copy(iterations),
        mean_newton_iterations=sum(iterations) / length(iterations),
        max_newton_iterations=maximum(iterations),
        newton_failed_steps=Int(failed_steps),
    )
end

function _fixed_solver_stats(times::AbstractVector{<:Real}, newton_iterations::Vector{Int}, newton_failed_steps::Integer)
    step_sizes = diff(Float64.(times))
    return (
        accepted_steps=length(times) - 1,
        rejected_steps=0,
        min_dt=minimum(step_sizes),
        max_dt=maximum(step_sizes),
        final_dt=last(step_sizes),
        max_fractional_change=NaN,
        max_absolute_change=NaN,
        reached_dt_min=false,
        _newton_iteration_summary(newton_iterations, newton_failed_steps)...,
    )
end

#=
    solve_network(network, Y0, tspan, dt, rho, T9; method=:rk4, rate_multipliers=nothing, clamp_negative=true)

Evolve a single-zone reaction network in time.

This solves the ordinary differential equation system `dY/dt = f(Y, rho, T9)`.
For a single-zone post-processing network there are no spatial derivatives, so
this is an ODE system rather than a PDE.

Arguments:
- `network`: a `ReactionNetwork`.
- `Y0`: initial abundance vector ordered like `network.species`.
- `tspan`: `(t_start, t_end)` in seconds.
- `dt`: fixed timestep in seconds.
- `rho`: density in g cm^-3, or a function `rho(t)`.
- `T9`: temperature in GK, or a function `T9(t)`.

Supported methods are `:euler`, `:rk4`, and `:backward_euler`. RK4 is usually
more accurate for the same timestep, while backward Euler is a dependency-free
implicit option for stiffer exploratory runs.

Returns `(times, Y_history)`, where `Y_history[i, j]` is the abundance of species
`network.species[j]` at `times[i]`.
=#
function solve_network(
    network::ReactionNetwork,
    Y0::AbstractVector{<:Real},
    tspan::Tuple{<:Real,<:Real},
    dt::Real,
    rho,
    T9;
    method::Symbol=:rk4,
    rate_multipliers=nothing,
    rate_p_values=nothing,
    clamp_negative::Bool=true,
    screening=nothing,
    newton_tolerance::Real=1.0e-10,
    max_newton_iterations::Integer=20,
    finite_difference_epsilon::Real=sqrt(eps(Float64)),
    jacobian::Symbol=:analytic,
    return_stats::Bool=false,
)
    times = _time_grid(tspan, dt)
    Y = _checked_initial_abundances(Y0, network)
    Y_history = Matrix{Float64}(undef, length(times), length(Y))
    Y_history[1, :] .= Y
    newton_iterations = Int[]
    newton_failed_steps = 0

    for step_index in 1:(length(times)-1)
        t = times[step_index]
        step_dt = times[step_index+1] - t

        if method == :euler
            Y = _euler_step(network, Y, t, step_dt, rho, T9; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)
        elseif method == :rk4
            Y = _rk4_step(network, Y, t, step_dt, rho, T9; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)
        elseif method == :backward_euler
            Y = _backward_euler_step(
                network,
                Y,
                t,
                step_dt,
                rho,
                T9;
                rate_multipliers=rate_multipliers,
                rate_p_values=rate_p_values,
                screening=screening,
                newton_tolerance=newton_tolerance,
                max_newton_iterations=max_newton_iterations,
                finite_difference_epsilon=finite_difference_epsilon,
                jacobian=jacobian,
                newton_iterations=newton_iterations,
            )
        else
            throw(ArgumentError("unsupported method $method; use :euler, :rk4, or :backward_euler"))
        end

        if clamp_negative
            for i in eachindex(Y)
                Y[i] < 0.0 && (Y[i] = 0.0)
            end
        end

        Y_history[step_index+1, :] .= Y
    end

    stats = _fixed_solver_stats(times, newton_iterations, newton_failed_steps)
    return return_stats ? (times, Y_history, stats) : (times, Y_history)
end

function _single_step(
    network::ReactionNetwork,
    Y::Vector{Float64},
    t::Float64,
    dt::Float64,
    rho,
    T9,
    method::Symbol;
    rate_multipliers=nothing,
    rate_p_values=nothing,
    screening=nothing,
    newton_tolerance::Real=1.0e-10,
    max_newton_iterations::Integer=20,
    finite_difference_epsilon::Real=sqrt(eps(Float64)),
    jacobian::Symbol=:analytic,
    newton_iterations=nothing,
)
    if method == :euler
        return _euler_step(network, Y, t, dt, rho, T9; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)
    elseif method == :rk4
        return _rk4_step(network, Y, t, dt, rho, T9; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)
    elseif method == :backward_euler
        return _backward_euler_step(
            network,
            Y,
            t,
            dt,
            rho,
            T9;
            rate_multipliers=rate_multipliers,
            rate_p_values=rate_p_values,
            screening=screening,
            newton_tolerance=newton_tolerance,
            max_newton_iterations=max_newton_iterations,
            finite_difference_epsilon=finite_difference_epsilon,
            jacobian=jacobian,
            newton_iterations=newton_iterations,
        )
    end
    throw(ArgumentError("unsupported method $method; use :euler, :rk4, or :backward_euler"))
end

function _max_fractional_change(Y::Vector{Float64}, Y_next::Vector{Float64}, abundance_floor::Float64)
    max_change = 0.0
    for i in eachindex(Y)
        scale = max(abs(Y[i]), abundance_floor)
        max_change = max(max_change, abs(Y_next[i] - Y[i]) / scale)
    end
    return max_change
end

function _max_absolute_change(Y::Vector{Float64}, Y_next::Vector{Float64})
    max_change = 0.0
    for i in eachindex(Y)
        max_change = max(max_change, abs(Y_next[i] - Y[i]))
    end
    return max_change
end

function _adaptive_factor(
    fractional_change::Float64,
    max_fractional_change::Float64,
    absolute_change::Float64,
    max_absolute_change::Float64,
    safety::Float64,
    shrink_factor::Float64,
    growth_factor::Float64,
)
    factor = growth_factor
    if fractional_change > 0.0
        factor = min(factor, safety * max_fractional_change / fractional_change)
    end
    if isfinite(max_absolute_change) && absolute_change > 0.0
        factor = min(factor, safety * max_absolute_change / absolute_change)
    end
    return clamp(factor, shrink_factor, growth_factor)
end

#=
    solve_network_adaptive(network, Y0, tspan, dt_initial, rho, T9; ...)

Evolve a network with simple adaptive explicit timestepping. A proposed step is
accepted when the maximum fractional abundance change is below
`max_fractional_change`, using `abundance_floor` to avoid division by zero for
trace species. If `max_absolute_change` is finite, the step must also satisfy
that absolute abundance-change limit.

This is still an explicit method, not a stiff implicit solver, but it is safer
than a fixed timestep for exploratory post-processing.
=#
function solve_network_adaptive(
    network::ReactionNetwork,
    Y0::AbstractVector{<:Real},
    tspan::Tuple{<:Real,<:Real},
    dt_initial::Real,
    rho,
    T9;
    method::Symbol=:rk4,
    max_fractional_change::Real=0.05,
    max_absolute_change::Real=Inf,
    abundance_floor::Real=1.0e-30,
    dt_min::Real=1.0e-12,
    dt_max::Real=Inf,
    safety::Real=0.8,
    growth_factor::Real=2.0,
    shrink_factor::Real=0.5,
    max_steps::Integer=1_000_000,
    rate_multipliers=nothing,
    rate_p_values=nothing,
    clamp_negative::Bool=true,
    screening=nothing,
    return_stats::Bool=false,
    newton_tolerance::Real=1.0e-10,
    max_newton_iterations::Integer=20,
    finite_difference_epsilon::Real=sqrt(eps(Float64)),
    jacobian::Symbol=:analytic,
)
    t_start = Float64(tspan[1])
    t_end = Float64(tspan[2])
    dt = min(Float64(dt_initial), Float64(dt_max))
    _validate_time_inputs(t_start, t_end, dt)
    max_fractional_change > 0.0 || throw(ArgumentError("max_fractional_change must be positive"))
    max_absolute_change > 0.0 || throw(ArgumentError("max_absolute_change must be positive"))
    abundance_floor > 0.0 || throw(ArgumentError("abundance_floor must be positive"))
    dt_min > 0.0 || throw(ArgumentError("dt_min must be positive"))
    dt_max > 0.0 || throw(ArgumentError("dt_max must be positive"))
    safety > 0.0 || throw(ArgumentError("safety must be positive"))
    growth_factor > 1.0 || throw(ArgumentError("growth_factor must be greater than 1"))
    0.0 < shrink_factor < 1.0 || throw(ArgumentError("shrink_factor must be between 0 and 1"))
    max_steps > 0 || throw(ArgumentError("max_steps must be positive"))

    Y = _checked_initial_abundances(Y0, network)
    times = Float64[t_start]
    history_rows = Vector{Float64}[copy(Y)]
    t = t_start
    accepted_steps = 0
    rejected_steps = 0
    min_dt_used = Inf
    max_dt_used = 0.0
    max_fractional_change_seen = 0.0
    max_absolute_change_seen = 0.0
    reached_dt_min = false
    newton_iterations = Int[]
    newton_failed_steps = 0

    while t < t_end && accepted_steps < max_steps
        step_dt = min(dt, t_end - t)
        proposed_newton_iterations = Int[]
        Y_next = try
            _single_step(
                network,
                Y,
                t,
                step_dt,
                rho,
                T9,
                method;
                rate_multipliers=rate_multipliers,
                rate_p_values=rate_p_values,
                screening=screening,
                newton_tolerance=newton_tolerance,
                max_newton_iterations=max_newton_iterations,
                finite_difference_epsilon=finite_difference_epsilon,
                jacobian=jacobian,
                newton_iterations=proposed_newton_iterations,
            )
        catch err
            if err isa NewtonConvergenceError
                rejected_steps += 1
                newton_failed_steps += 1
                at_dt_min = step_dt <= dt_min
                at_dt_min && throw(ArgumentError("backward Euler Newton iteration failed at dt_min=$dt_min: $err"))
                dt = max(step_dt * shrink_factor, Float64(dt_min))
                continue
            end
            rethrow()
        end

        if clamp_negative
            for i in eachindex(Y_next)
                Y_next[i] < 0.0 && (Y_next[i] = 0.0)
            end
        end

        fractional_change = _max_fractional_change(Y, Y_next, Float64(abundance_floor))
        absolute_change = _max_absolute_change(Y, Y_next)
        max_fractional_change_seen = max(max_fractional_change_seen, fractional_change)
        max_absolute_change_seen = max(max_absolute_change_seen, absolute_change)

        fractional_ok = fractional_change <= max_fractional_change
        absolute_ok = absolute_change <= max_absolute_change
        at_dt_min = step_dt <= dt_min
        if (fractional_ok && absolute_ok) || at_dt_min
            t += step_dt
            Y = Y_next
            push!(times, t)
            push!(history_rows, copy(Y))
            append!(newton_iterations, proposed_newton_iterations)
            accepted_steps += 1
            min_dt_used = min(min_dt_used, step_dt)
            max_dt_used = max(max_dt_used, step_dt)
            reached_dt_min |= at_dt_min

            if fractional_change == 0.0 && absolute_change == 0.0
                dt = min(step_dt * growth_factor, Float64(dt_max))
            else
                factor = _adaptive_factor(
                    fractional_change,
                    Float64(max_fractional_change),
                    absolute_change,
                    Float64(max_absolute_change),
                    Float64(safety),
                    Float64(shrink_factor),
                    Float64(growth_factor),
                )
                dt = min(max(step_dt * factor, Float64(dt_min)), Float64(dt_max))
            end
        else
            rejected_steps += 1
            dt = max(step_dt * shrink_factor, Float64(dt_min))
            if at_dt_min
                throw(ArgumentError("adaptive timestep reached dt_min=$dt_min but proposed fractional change=$fractional_change and absolute change=$absolute_change exceed limits"))
            end
        end
    end

    accepted_steps < max_steps || throw(ArgumentError("adaptive solver exceeded max_steps=$max_steps"))

    history = Matrix{Float64}(undef, length(history_rows), length(Y))
    for (i, row) in pairs(history_rows)
        history[i, :] .= row
    end

    stats = (
        accepted_steps=accepted_steps,
        rejected_steps=rejected_steps,
        min_dt=min_dt_used,
        max_dt=max_dt_used,
        final_dt=dt,
        max_fractional_change=max_fractional_change_seen,
        max_absolute_change=max_absolute_change_seen,
        reached_dt_min=reached_dt_min,
        _newton_iteration_summary(newton_iterations, newton_failed_steps)...,
    )

    return return_stats ? (times, history, stats) : (times, history)
end

"""
    run_monte_carlo(network, Y0, tspan, dt, rho, T9; nruns, seed=nothing, method=:rk4, store_histories=false)

Run repeated single-zone network calculations with STARLIB lognormal rate
sampling. For each run and each reaction, a random `p ~ Normal(0, 1)` is drawn
and held fixed for that reaction throughout the run:

    sampled_rate(T) = recommended_rate(T) * factor_uncertainty(T)^p

Returns a named tuple containing final abundances, sampled `p` values, and
optionally all abundance histories.
"""
function run_monte_carlo(
    network::ReactionNetwork,
    Y0::AbstractVector{<:Real},
    tspan::Tuple{<:Real,<:Real},
    dt::Real,
    rho,
    T9;
    nruns::Integer,
    seed=nothing,
    method::Symbol=:rk4,
    rate_multipliers=nothing,
    clamp_negative::Bool=true,
    store_histories::Bool=false,
    screening=nothing,
)
    nruns > 0 || throw(ArgumentError("nruns must be positive"))

    rng = seed === nothing ? Random.default_rng() : Random.MersenneTwister(seed)
    nreactions = length(network.reactions)
    nspecies = length(network.species)
    sampled_p_values = randn(rng, nruns, nreactions)
    final_abundances = Matrix{Float64}(undef, nruns, nspecies)
    histories = Matrix{Float64}[]
    saved_times = Float64[]

    for run_index in 1:nruns
        p_values = view(sampled_p_values, run_index, :)
        times, history = solve_network(
            network,
            Y0,
            tspan,
            dt,
            rho,
            T9;
            method=method,
            rate_multipliers=rate_multipliers,
            rate_p_values=p_values,
            clamp_negative=clamp_negative,
            screening=screening,
        )

        if run_index == 1
            saved_times = collect(times)
        end

        final_abundances[run_index, :] .= history[end, :]
        store_histories && push!(histories, history)
    end

    return (
        times=saved_times,
        final_abundances=final_abundances,
        rate_p_values=sampled_p_values,
        histories=histories,
        species=copy(network.species),
        reactions=[reaction_string(reaction) for reaction in network.reactions],
    )
end

"""
    solve_single_zone(tables, labels, X0, tspan, dt, rho, T9; adaptive=true, ...)

Build and run a single-zone post-processing network from user-facing inputs.

This is the convenience workflow for interactive use:
- select STARLIB reactions by label
- infer and validate the network
- convert initial mass fractions `X0` to abundances
- run fixed-step or adaptive integration
- return mass-fraction diagnostics with the raw abundance history
"""
function solve_single_zone(
    tables::AbstractVector{ReactionRateTable},
    labels::AbstractVector{<:AbstractString},
    X0::AbstractDict,
    tspan::Tuple{<:Real,<:Real},
    dt::Real,
    rho,
    T9;
    species=nothing,
    source=nothing,
    on_multiple::Symbol=:error,
    adaptive::Bool=true,
    method::Symbol=:rk4,
    normalize_mass_fractions::Bool=false,
    check_mass_fraction_sum::Bool=false,
    mass_fraction_atol::Real=1.0e-8,
    validate::Bool=true,
    max_fractional_change::Real=0.05,
    max_absolute_change::Real=Inf,
    abundance_floor::Real=1.0e-30,
    dt_min::Real=1.0e-12,
    dt_max::Real=Inf,
    safety::Real=0.8,
    growth_factor::Real=2.0,
    shrink_factor::Real=0.5,
    max_steps::Integer=1_000_000,
    rate_multipliers=nothing,
    rate_p_values=nothing,
    clamp_negative::Bool=true,
    screening=nothing,
    newton_tolerance::Real=1.0e-10,
    max_newton_iterations::Integer=20,
    finite_difference_epsilon::Real=sqrt(eps(Float64)),
    jacobian::Symbol=:analytic,
)
    network = network_from_labels(tables, labels; species=species, source=source, on_multiple=on_multiple)
    validation = validate ? network_validation_report(network; throw_on_error=true) : network_validation_report(network)
    Y0 = abundances_from_mass_fractions(
        network,
        X0;
        normalize=normalize_mass_fractions,
        check_sum=check_mass_fraction_sum,
        atol=mass_fraction_atol,
    )

    times, history, solver_stats = if adaptive
        solve_network_adaptive(
            network,
            Y0,
            tspan,
            dt,
            rho,
            T9;
            method=method,
            max_fractional_change=max_fractional_change,
            max_absolute_change=max_absolute_change,
            abundance_floor=abundance_floor,
            dt_min=dt_min,
            dt_max=dt_max,
            safety=safety,
            growth_factor=growth_factor,
            shrink_factor=shrink_factor,
            max_steps=max_steps,
            rate_multipliers=rate_multipliers,
            rate_p_values=rate_p_values,
            clamp_negative=clamp_negative,
            screening=screening,
            return_stats=true,
            newton_tolerance=newton_tolerance,
            max_newton_iterations=max_newton_iterations,
            finite_difference_epsilon=finite_difference_epsilon,
            jacobian=jacobian,
        )
    else
        solve_network(
            network,
            Y0,
            tspan,
            dt,
            rho,
            T9;
            method=method,
            rate_multipliers=rate_multipliers,
            rate_p_values=rate_p_values,
            clamp_negative=clamp_negative,
            screening=screening,
            newton_tolerance=newton_tolerance,
            max_newton_iterations=max_newton_iterations,
            finite_difference_epsilon=finite_difference_epsilon,
            jacobian=jacobian,
            return_stats=true,
        )
    end

    flux_history = reaction_flux_history(
        network,
        history,
        times,
        rho,
        T9;
        rate_multipliers=rate_multipliers,
        rate_p_values=rate_p_values,
        screening=screening,
    )
    epsilon_history = energy_generation_history(
        network,
        history,
        times,
        rho,
        T9;
        rate_multipliers=rate_multipliers,
        rate_p_values=rate_p_values,
        screening=screening,
    )

    return (
        network=network,
        validation=validation,
        times=times,
        abundances=history,
        mass_fraction_history=mass_fraction_history(network, history),
        mass_fraction_drift=mass_fraction_drift(network, history),
        abundance_diagnostics=abundance_diagnostics(network, history),
        initial_mass_fractions=mass_fractions_from_abundances(network, view(history, 1, :)),
        final_mass_fractions=mass_fractions_from_abundances(network, view(history, size(history, 1), :)),
        reaction_fluxes=flux_history,
        integrated_fluxes=integrated_fluxes(times, flux_history),
        energy_generation=epsilon_history,
        integrated_energy_generation=integrated_energy_generation(times, epsilon_history),
        solver_stats=solver_stats,
        adaptive=adaptive,
    )
end

"""
    solve_network_fbdf(network, Y0, tspan, rho, T9; kwargs...)

Solve the network with the variable-order stiff `FBDF` integrator from
OrdinaryDiffEqBDF, using ReacNetJl's cached RHS and analytic Jacobian.
Higher-order than backward Euler, so it takes far fewer steps at the same
accuracy.

This function is provided by a package extension: install and load the solver
package first with `import Pkg; Pkg.add("OrdinaryDiffEqBDF")` and
`using OrdinaryDiffEqBDF`. Returns `(times, history, stats)` like
`solve_network_adaptive`. Also available through `run_ppn(...; method=:fbdf)`.
"""
function solve_network_fbdf(network::ReactionNetwork, Y0, tspan, rho, T9; kwargs...)
    error("solve_network_fbdf requires the OrdinaryDiffEqBDF extension; run `import Pkg; Pkg.add(\"OrdinaryDiffEqBDF\")` and add `using OrdinaryDiffEqBDF` before calling it")
end

"""
    run_ppn(trajectory_path, abundance_path;
            rates=:starlib, screening=:weak, neutron_captures=true, kwargs...)

Run a complete single-zone post-processing nucleosynthesis calculation from a
trajectory file and an initial-abundance file, in one call.

The pipeline is: read the trajectory and abundances, build the H-Ca nova
network from the chosen rate library, add reverse rates, validate the
network, and integrate the abundances over the full trajectory with the
adaptive backward-Euler solver and analytic Jacobian.

# Arguments
- `trajectory_path`: trajectory file with `AGEUNIT`/`TUNIT`/`RHOUNIT` metadata.
- `abundance_path`: initial abundance table (`Z name A X` rows).

# Keywords
- `rates`: `:starlib` (default) or `:iliadis2002` for the NACRE (A < 20) plus
  Iliadis et al. 2001 (A = 20-40) REACLIB baseline of the 2002 nova
  sensitivity study.
- `tables`: pass a prebuilt `Vector{ReactionRateTable}` to skip library
  loading (overrides `rates`).
- `screening`: `nothing`, `:weak`, or `:chugunov`.
- `neutron_captures`: include neutron-induced reactions in the network
  selection (default `true`).
- `partition_functions`: `:auto` (default; uses `data/winvne_v2.0.dat` when
  present), `nothing`, or a `PartitionFunctionTable` — applied to generated
  detailed-balance reverse rates.
- `jacobian`: `:analytic` (default) or `:finite_difference`.
- `method`: `:backward_euler` (default), `:euler`, `:rk4`, or `:fbdf` for the
  high-order stiff integrator (requires `using OrdinaryDiffEqBDF`).
- `dt_initial`, `dt_min`, `dt_max`: solver step controls; by default they are
  chosen from the trajectory duration like the nova example driver.
- `max_fractional_change`, `max_absolute_change`, `abundance_floor`,
  `max_newton_iterations`, `max_steps`: adaptive controller settings,
  defaulting to the values validated by the nova example driver.

# Returns
A named tuple with the `network`, `trajectory`, solution `times` and
abundance `history`, `initial_mass_fractions` and `final_mass_fractions`
dictionaries, `inert_mass_fractions` for species outside the network,
`solver_stats`, `reverse_summary`, the `rate_policy_report` (for
`rates=:iliadis2002`), and the network `validation` report.
"""
function run_ppn(
    trajectory_path::AbstractString,
    abundance_path::AbstractString;
    rates::Symbol=:starlib,
    tables=nothing,
    screening=:weak,
    neutron_captures::Bool=true,
    generate_detailed_balance::Bool=true,
    partition_functions=:auto,
    jacobian::Symbol=:analytic,
    method::Symbol=:backward_euler,
    max_fractional_change::Real=0.50,
    max_absolute_change::Real=1.0e-4,
    abundance_floor::Real=1.0e-8,
    max_newton_iterations::Integer=80,
    dt_initial=nothing,
    dt_min=nothing,
    dt_max=nothing,
    max_steps::Integer=1_000_000,
)
    trajectory = read_trajectory(trajectory_path)
    profiles = trajectory_profiles(trajectory)
    X_raw = read_initial_abundances(abundance_path)
    X_normalized = read_initial_abundances(abundance_path; normalize=true)

    rate_policy_report = nothing
    if tables === nothing
        if rates == :iliadis2002
            selection = iliadis2002_rate_tables(; include_reverse=true)
            tables = selection.tables
            rate_policy_report = selection.report
        elseif rates == :starlib
            tables = read_starlib()
        else
            throw(ArgumentError("unsupported rates=$rates; use :starlib or :iliadis2002"))
        end
    end

    projectiles = neutron_captures ? ("p", "he4", "he3", "d", "n") : ("p", "he4", "he3", "d")
    forward_tables = select_h_ca_reaction_tables(tables, keys(X_raw); projectiles=projectiles)

    pf = partition_functions
    if pf === :auto
        pf = isfile(DEFAULT_WINVNE_PATH) ? read_winvne() : nothing
    end
    reverse_summary = add_reverse_reaction_tables(
        tables,
        forward_tables;
        generate_detailed_balance=generate_detailed_balance,
        partition_functions=pf,
    )
    network = network_from_tables(reverse_summary.tables)
    validation = network_validation_report(network; throw_on_error=true)

    X0 = Dict(name => value for (name, value) in X_normalized if haskey(network.species_index, name))
    inert_mass_fractions = Dict(name => value for (name, value) in X_normalized if !haskey(network.species_index, name))
    Y0 = abundances_from_mass_fractions(network, X0)

    t_start = first(trajectory.time)
    t_end = last(trajectory.time)
    duration = t_end - t_start
    step_initial = dt_initial === nothing ? (duration > 100.0 ? 1.0 : 0.02) : Float64(dt_initial)
    step_min = dt_min === nothing ? (duration > 100.0 ? 1.0e-8 : 1.0e-10) : Float64(dt_min)
    step_max = dt_max === nothing ? (duration > 100.0 ? 20.0 : 0.05) : Float64(dt_max)

    times, history, solver_stats = if method == :fbdf
        solve_network_fbdf(
            network,
            Y0,
            (t_start, t_end),
            profiles.rho,
            profiles.T9;
            screening=screening,
        )
    else
        solve_network_adaptive(
            network,
            Y0,
            (t_start, t_end),
            step_initial,
            profiles.rho,
            profiles.T9;
            method=method,
            screening=screening,
            jacobian=jacobian,
            max_fractional_change=max_fractional_change,
            max_absolute_change=max_absolute_change,
            abundance_floor=abundance_floor,
            max_newton_iterations=max_newton_iterations,
            dt_min=step_min,
            dt_max=step_max,
            max_steps=max_steps,
            return_stats=true,
        )
    end

    return (
        network=network,
        trajectory=trajectory,
        times=times,
        history=history,
        initial_mass_fractions=mass_fractions_from_abundances(network, view(history, 1, :)),
        final_mass_fractions=mass_fractions_from_abundances(network, view(history, size(history, 1), :)),
        inert_mass_fractions=inert_mass_fractions,
        mass_fraction_drift=mass_fraction_drift(network, history),
        solver_stats=solver_stats,
        reverse_summary=(explicit=reverse_summary.explicit, generated=reverse_summary.generated, missing=reverse_summary.missing),
        rate_policy_report=rate_policy_report,
        validation=validation,
    )
end

# Precompile the solver hot path on a miniature network so first use in a
# fresh session skips most of the compilation latency.
@setup_workload begin
    _pc_grid = STARLIB_T9_GRID
    _pc_unit = ones(length(_pc_grid))
    _pc_tables = [
        ReactionRateTable(4, ["p", "o17"], ["f18"], "pc", 5.607, _pc_grid, fill(2.0e4, length(_pc_grid)), _pc_unit),
        ReactionRateTable(5, ["p", "f18"], ["he4", "o15"], "pc", 2.882, _pc_grid, fill(6.0e4, length(_pc_grid)), _pc_unit),
        ReactionRateTable(1, ["o15"], ["n15"], "pcw", 2.754, _pc_grid, fill(4.5e-3, length(_pc_grid)), _pc_unit),
    ]

    @compile_workload begin
        _pc_network = network_from_tables(_pc_tables)
        _pc_X0 = Dict("p" => 0.7, "o17" => 0.1, "f18" => 1.0e-4, "he4" => 0.2, "o15" => 0.0, "n15" => 0.0)
        _pc_Y0 = abundances_from_mass_fractions(_pc_network, _pc_X0; normalize=true)
        for _pc_screening in (nothing, :weak, :chugunov)
            network_rhs(_pc_Y0, _pc_network, 500.0, 0.25; screening=_pc_screening)
            solve_network_adaptive(
                _pc_network, _pc_Y0, (0.0, 1.0e-6), 1.0e-7, 500.0, 0.25;
                method=:backward_euler, screening=_pc_screening, return_stats=true,
                abundance_floor=1.0e-8, max_newton_iterations=80,
            )
        end
        solve_network(_pc_network, _pc_Y0, (0.0, 1.0e-7), 1.0e-8, 500.0, 0.25; method=:rk4, screening=:weak)
        mass_fractions_from_abundances(_pc_network, _pc_Y0)
    end
end

end # module
