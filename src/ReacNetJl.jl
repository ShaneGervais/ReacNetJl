module ReacNetJl

using LinearAlgebra
using Random

export Species,
    Trajectory,
    ReactionRateTable,
    Reaction,
    ReactionNetwork,
    species_from_name,
    abundance_from_mass_fraction,
    mass_fraction_from_abundance,
    normalize_species_name,
    parse_reaction_label,
    read_starlib,
    read_trajectory,
    trajectory_profiles,
    starlib_chapter_report,
    find_rate,
    find_reverse_rate,
    reaction_from_label,
    network_from_labels,
    weak_screening_multiplier,
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
    species_flux_balance,
    reaction_edges,
    reaction_conservation,
    network_validation_report,
    network_rhs,
    solve_network,
    solve_network_adaptive,
    run_monte_carlo,
    solve_single_zone

#=
    Species

A nuclear species identified by charge number `Z`, mass number `A`, and a display
`name`, such as `Species("f18", 9, 18)`.
=#
struct Species
    name::String
    Z::Int
    A::Int
end

#=
    Trajectory

A single-zone thermodynamic trajectory with time in seconds, temperature `T9` in
GK, and density `rho` in g cm^-3.
=#
struct Trajectory
    time::Vector{Float64}
    T9::Vector{Float64}
    rho::Vector{Float64}
end

#=
    ReactionRateTable

A temperature-dependent thermonuclear reaction rate from STARLIB.

`T9` is the temperature grid in GK. `rate` is the recommended STARLIB rate.
`factor_uncertainty` is the STARLIB multiplicative factor uncertainty.
=#
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

#=
    Reaction

A network reaction built from a `ReactionRateTable`.

The first version stores only the reaction participants and a STARLIB rate table.
Later we can add screening, reverse-rate handling, and uncertainty sampling here
without changing the network RHS interface.
=#
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
end

#=
    ReactionNetwork

A single-zone nuclear reaction network.

`species` fixes the order of the abundance vector `Y`. `species_index` maps each
species name to its position in `Y`, and `reactions` stores the reactions that
contribute to `dY/dt`.
=#
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

const DEFAULT_STARLIB_PATH = joinpath(dirname(@__DIR__), "starlib_v610_120222.dat/starlib.dat")
const STARLIB_ROWS_PER_REACTION = 60

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

#=
    normalize_species_name(name)

Normalize common isotope notation to STARLIB-style lowercase names.

Examples:
- `"18F" -> "f18"`
- `"4He" -> "he4"`
- `"α" -> "he4"`
- `"γ" -> "gamma"`
=#
function normalize_species_name(name::AbstractString)
    s = lowercase(strip(name))
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

#=
    species_from_name(name)

Create a `Species` object from a normalized or common isotope name.

Examples:
- `species_from_name("18F") == Species("f18", 9, 18)`
- `species_from_name("he4") == Species("he4", 2, 4)`
- `species_from_name("p") == Species("p", 1, 1)`
=#
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
=#
function _is_beta_token(token::AbstractString)
    s = lowercase(strip(token))
    s = replace(s, " " => "", "β" => "beta")
    return s in ("beta+", "b+", "β+", "betaplus", "ec", "electroncapture", "beta-", "b-", "β-", "betaminus")
end

function parse_reaction_label(label::AbstractString)
    decay_match = match(r"^\s*([^\(]+)\(([^,\)]+)\)(.+?)\s*$", label)
    if decay_match !== nothing && _is_beta_token(decay_match.captures[2])
        parent = normalize_species_name(decay_match.captures[1])
        daughter = normalize_species_name(decay_match.captures[3])
        return [parent], [daughter]
    end

    m = match(r"^\s*([^\(]+)\(([^,]+),([^\)]+)\)(.+?)\s*$", label)
    m === nothing && throw(ArgumentError("reaction label must look like target(projectile,ejectile)product or parent(β+)daughter"))

    target = normalize_species_name(m.captures[1])
    projectile = normalize_species_name(m.captures[2])
    ejectile = normalize_species_name(m.captures[3])
    product = normalize_species_name(m.captures[4])

    reactants = [projectile, target]
    products = ejectile == "gamma" ? [product] : [ejectile, product]
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
    if chapter == 1 && length(species) == 2
        return [species[1]], [species[2]]
    elseif chapter in (2, 4) && length(species) == 3
        return species[1:2], [species[3]]
    elseif chapter in (4, 5) && length(species) == 4
        return species[1:2], species[3:4]
    end

    # Conservative fallback: keep the raw STARLIB order if this chapter is not
    # supported yet. We will expand this as the network grows.
    return species, String[]
end

function _supported_starlib_layout(chapter::Int, nspecies::Int)
    return (chapter == 1 && nspecies == 2) ||
           (chapter in (2, 4) && nspecies == 3) ||
           (chapter in (4, 5) && nspecies == 4)
end

function _supported_rate_table(table::ReactionRateTable)
    return !isempty(table.reactants) && !isempty(table.products)
end

#=
    starlib_chapter_report(tables)

Summarize which STARLIB chapter layouts were parsed into supported reactant and
product bookkeeping. Unsupported rows are kept by `read_starlib`, but they are
not suitable for network construction until their chapter layout is implemented.
=#
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

function read_starlib(path::AbstractString=DEFAULT_STARLIB_PATH; warn_unsupported::Bool=false)
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

#=
    find_rate(tables, label; source=nothing)

Find STARLIB rate tables matching a reaction label like `"18F(p,α)15O"`.
Returns all matching tables because STARLIB can contain multiple sources.
=#
function find_rate(tables::AbstractVector{ReactionRateTable}, label::AbstractString; source=nothing)
    reactants, products = parse_reaction_label(label)
    matches = filter(t -> t.reactants == reactants && t.products == products, tables)

    source === nothing && return matches
    wanted = lowercase(strip(string(source)))
    return filter(t -> lowercase(t.source) == wanted, matches)
end

#=
    find_reverse_rate(tables, label; source=nothing)

Find STARLIB tables whose parsed reactants/products are the exact reverse of a
reaction label. This detects explicit reverse rates already present in STARLIB;
it does not synthesize reciprocal-rule reverse rates.
=#
function find_reverse_rate(tables::AbstractVector{ReactionRateTable}, label::AbstractString; source=nothing)
    reactants, products = parse_reaction_label(label)
    matches = filter(t -> t.reactants == products && t.products == reactants, tables)

    source === nothing && return matches
    wanted = lowercase(strip(string(source)))
    return filter(t -> lowercase(t.source) == wanted, matches)
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

#=
    reaction_from_label(tables, label; source=nothing, on_multiple=:error)

Find a STARLIB rate table by reaction label and wrap it as a `Reaction`.
If multiple entries match, choose a `source`, or set `on_multiple=:first`.
=#
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

    values[i] > 0.0 || throw(ArgumentError("cannot log-interpolate non-positive $value_name at lower grid point"))
    values[i+1] > 0.0 || throw(ArgumentError("cannot log-interpolate non-positive $value_name at upper grid point"))

    x0 = log(grid[i])
    x1 = log(grid[i+1])
    y0 = log(values[i])
    y1 = log(values[i+1])
    weight = (log(T) - x0) / (x1 - x0)
    return exp((1 - weight) * y0 + weight * y1)
end

function interpolate_rate(table::ReactionRateTable, T9::Real)
    return _interpolate_loglog(table.T9, table.rate, T9; value_name="rate")
end

#=
    interpolate_factor_uncertainty(table, T9)

Interpolate STARLIB's multiplicative factor uncertainty at temperature `T9`.
The interpolation is log-log, matching `interpolate_rate`.
=#
function interpolate_factor_uncertainty(table::ReactionRateTable, T9::Real)
    return _interpolate_loglog(table.T9, table.factor_uncertainty, T9; value_name="factor uncertainty")
end

#=
    sampled_interpolate_rate(table, T9, p)

Return a STARLIB lognormal sampled reaction rate:

    sampled_rate = recommended_rate * factor_uncertainty^p

where `p` is usually a standard normal deviate held fixed for a reaction during
one Monte Carlo network run.
=#
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

#=
    weak_screening_multiplier(network, reaction, Y, rho, T9)

Return an approximate weak-screening multiplier for charged-particle reactions.
This is a Salpeter-style diagnostic multiplier using the current abundance
composition. Reactions with fewer than two charged reactants return `1.0`.
=#
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

function _screening_multiplier(screening, network::ReactionNetwork, reaction::Reaction, Y::AbstractVector{<:Real}, rho::Real, T9::Real)
    if screening === nothing || screening === false
        return 1.0
    elseif screening == :weak
        return weak_screening_multiplier(network, reaction, Y, rho, T9)
    elseif screening isa Function
        return Float64(screening(network, reaction, Y, rho, T9))
    end

    throw(ArgumentError("unsupported screening=$screening; use nothing, :weak, or a function `(network, reaction, Y, rho, T9) -> multiplier`"))
end

#=
    reaction_flux(reaction, Y, species_index, rho, T9; rate_multiplier=1.0)

Calculate the abundance flux for one reaction.

`Y` is the abundance vector. `species_index` maps species names to positions in
`Y`. `rho` is mass density in g cm^-3, and `T9` is temperature in GK.

For a one-body reaction, the flux is `rate * Yᵢ`. For a two-body reaction using
STARLIB's usual `N_A <σv>` rate, the flux is `rho * rate * Yᵢ * Yⱼ`, with the
standard symmetry correction for identical reactants.
=#
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

        for i in eachindex(dYdt)
            dYdt[i] += compiled.stoichiometric_delta[i] * flux
        end
    end

    return dYdt
end

#=
    reaction_string(reaction)

Return a compact, normalized reaction label for display or diagnostics.

For common two-body STARLIB reactions, this returns labels like
`"f18(p,he4)o15"` or `"f18(p,γ)ne19"`.
=#
function reaction_string(reaction::Reaction)
    if length(reaction.reactants) == 2 && length(reaction.products) == 1
        projectile, target = reaction.reactants
        product = only(reaction.products)
        return "$(target)($(projectile),γ)$(product)"
    elseif length(reaction.reactants) == 2 && length(reaction.products) == 2
        projectile, target = reaction.reactants
        ejectile, product = reaction.products
        return "$(target)($(projectile),$(ejectile))$(product)"
    elseif length(reaction.reactants) == 1 && length(reaction.products) == 1
        parent = species_from_name(only(reaction.reactants))
        daughter = species_from_name(only(reaction.products))
        if parent.A == daughter.A && parent.Z == daughter.Z + 1
            return "$(parent.name)(β+)$(daughter.name)"
        elseif parent.A == daughter.A && parent.Z + 1 == daughter.Z
            return "$(parent.name)(β-)$(daughter.name)"
        end
        return "$(only(reaction.reactants))->$(only(reaction.products))"
    end

    return join(reaction.reactants, "+") * "->" * join(reaction.products, "+")
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
    species_flux_balance(network, Y, rho, T9; rate_multipliers=nothing)

Calculate instantaneous production, destruction, and net `dY/dt` contributions
for every species.

Returns `(production=..., destruction=..., net=...)`, with vectors ordered like
`network.species`.
=#
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

#=
    reaction_conservation(reaction)

Check baryon-number (`A`) and charge (`Z`) conservation for one reaction.
Returns a named tuple with reactant/product totals and boolean conservation flags.

Photons, if present, are treated as `A=0`, `Z=0`.
=#
function reaction_conservation(reaction::Reaction)
    reactant_A = _species_property_sum(reaction.reactants, :A)
    product_A = _species_property_sum(reaction.products, :A)
    reactant_Z = _species_property_sum(reaction.reactants, :Z)
    product_Z = _species_property_sum(reaction.products, :Z)
    delta_A = product_A - reactant_A
    delta_Z = product_Z - reactant_Z
    is_weak_decay = length(reaction.reactants) == 1 && length(reaction.products) == 1 && delta_A == 0 && abs(delta_Z) == 1

    return (
        reaction=reaction_string(reaction),
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

#=
    network_validation_report(network; throw_on_error=false)

Validate basic reaction-network bookkeeping and physics constraints.

Checks currently include:
- consistency of `species`, `species_info`, and `species_index`
- all reaction species are present in the network
- baryon-number conservation for every reaction
- charge conservation for every reaction

Returns a named tuple containing summary information, per-reaction reports, and
all issues. If `throw_on_error=true`, invalid networks raise an `ArgumentError`.
=#
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

function _euler_step(network::ReactionNetwork, Y::Vector{Float64}, t::Float64, dt::Float64, rho, T9; rate_multipliers=nothing, rate_p_values=nothing, screening=nothing)
    k1 = _rhs_at(network, Y, rho, T9, t; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)
    return Y .+ dt .* k1
end

function _rk4_step(network::ReactionNetwork, Y::Vector{Float64}, t::Float64, dt::Float64, rho, T9; rate_multipliers=nothing, rate_p_values=nothing, screening=nothing)
    k1 = _rhs_at(network, Y, rho, T9, t; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)
    k2 = _rhs_at(network, Y .+ 0.5 * dt .* k1, rho, T9, t + 0.5 * dt; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)
    k3 = _rhs_at(network, Y .+ 0.5 * dt .* k2, rho, T9, t + 0.5 * dt; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)
    k4 = _rhs_at(network, Y .+ dt .* k3, rho, T9, t + dt; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)
    return Y .+ (dt / 6.0) .* (k1 .+ 2.0 .* k2 .+ 2.0 .* k3 .+ k4)
end

function _max_abs(values::AbstractVector{<:Real})
    maximum(abs, values; init=0.0)
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
    perturbed = copy(Y_next)

    for j in 1:n
        saved = perturbed[j]
        step = finite_difference_epsilon * max(abs(saved), 1.0)
        perturbed[j] = saved + step
        perturbed_residual = _backward_euler_residual(network, perturbed, Y, t_next, dt, rho, T9; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)
        jacobian[:, j] .= (perturbed_residual .- residual) ./ step
        perturbed[j] = saved
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
    newton_iterations=nothing,
)
    max_newton_iterations > 0 || throw(ArgumentError("max_newton_iterations must be positive"))
    newton_tolerance > 0.0 || throw(ArgumentError("newton_tolerance must be positive"))
    finite_difference_epsilon > 0.0 || throw(ArgumentError("finite_difference_epsilon must be positive"))

    t_next = t + dt
    Y_next = _euler_step(network, Y, t, dt, rho, T9; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)

    residual = _backward_euler_residual(network, Y_next, Y, t_next, dt, rho, T9; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)
    tolerance = Float64(newton_tolerance) * max(_max_abs(Y_next), 1.0)
    if _max_abs(residual) <= tolerance
        newton_iterations !== nothing && push!(newton_iterations, 0)
        return Y_next
    end

    for iteration in 1:max_newton_iterations
        jacobian = _backward_euler_jacobian(
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
        correction = jacobian \ (-residual)
        Y_next .+= correction

        residual = _backward_euler_residual(network, Y_next, Y, t_next, dt, rho, T9; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)
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

#=
    run_monte_carlo(network, Y0, tspan, dt, rho, T9; nruns, seed=nothing, method=:rk4, store_histories=false)

Run repeated single-zone network calculations with STARLIB lognormal rate
sampling. For each run and each reaction, a random `p ~ Normal(0, 1)` is drawn
and held fixed for that reaction throughout the run:

    sampled_rate(T) = recommended_rate(T) * factor_uncertainty(T)^p

Returns a named tuple containing final abundances, sampled `p` values, and
optionally all abundance histories.
=#
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

#=
    solve_single_zone(tables, labels, X0, tspan, dt, rho, T9; adaptive=true, ...)

Build and run a single-zone post-processing network from user-facing inputs.

This is the convenience workflow for interactive use:
- select STARLIB reactions by label
- infer and validate the network
- convert initial mass fractions `X0` to abundances
- run fixed-step or adaptive integration
- return mass-fraction diagnostics with the raw abundance history
=#
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
        solver_stats=solver_stats,
        adaptive=adaptive,
    )
end

end # module
