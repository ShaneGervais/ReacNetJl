module ReacNetJl

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
    find_rate,
    reaction_from_label,
    network_from_labels,
    abundances_from_mass_fractions,
    mass_fractions_from_abundances,
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
    run_monte_carlo

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

    return ReactionNetwork(normalized_species, species_info, species_index, normalized_reactions)
end

const DEFAULT_STARLIB_PATH = joinpath(dirname(@__DIR__), "starlib.dat")
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
    m === nothing && throw(ArgumentError("could not infer species information from '$name'"))

    symbol = m.captures[1]
    A = parse(Int, m.captures[2])
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

Lines beginning with `#` and blank lines are ignored.
=#
function read_trajectory(path::AbstractString)
    time = Float64[]
    T9 = Float64[]
    rho = Float64[]

    open(path, "r") do io
        for raw_line in eachline(io)
            line = strip(split(raw_line, '#'; limit=2)[1])
            isempty(line) && continue
            fields = split(replace(line, ',' => ' '))
            length(fields) >= 3 || throw(ArgumentError("trajectory row must contain at least three columns: $raw_line"))

            push!(time, parse(Float64, fields[1]))
            push!(T9, parse(Float64, fields[2]))
            push!(rho, parse(Float64, fields[3]))
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

function read_starlib(path::AbstractString=DEFAULT_STARLIB_PATH)
    tables = ReactionRateTable[]

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

function _species_counts(species::Vector{String})
    counts = Dict{String,Int}()
    for name in species
        counts[name] = get(counts, name, 0) + 1
    end
    return counts
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
    for count in values(_species_counts(species))
        factor *= _factorial(count)
    end
    return factor
end

function _species_index(index::AbstractDict, name::String)
    haskey(index, name) && return index[name]
    throw(ArgumentError("species '$name' is missing from the species index"))
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

        for (name, count) in _species_counts(reaction.reactants)
            dYdt[_species_index(species_index, name)] -= count * flux
        end

        for (name, count) in _species_counts(reaction.products)
            dYdt[_species_index(species_index, name)] += count * flux
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
)
    return network_rhs(Y, network.reactions, network.species_index, rho, T9; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values)
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
        fluxes[i] = reaction_flux(reaction, Y, network.species_index, rho, T9; rate_multiplier=multiplier, rate_p_value=p_value)
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
)
    length(times) == size(history, 1) || throw(ArgumentError("times length must match the number of history rows"))
    size(history, 2) == length(network.species) || throw(ArgumentError("history column count must match the number of network species"))

    flux_history = Matrix{Float64}(undef, length(times), length(network.reactions))
    for (n, t) in pairs(times)
        rho_t = _profile_value(rho, t)
        T9_t = _profile_value(T9, t)
        flux_history[n, :] .= reaction_fluxes(network, view(history, n, :), rho_t, T9_t; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values)
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
)
    fluxes = reaction_fluxes(network, Y, rho, T9; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values)
    production = zeros(Float64, length(network.species))
    destruction = zeros(Float64, length(network.species))

    for (reaction_index, reaction) in pairs(network.reactions)
        flux = fluxes[reaction_index]

        for (name, count) in _species_counts(reaction.products)
            production[_species_index(network.species_index, name)] += count * flux
        end

        for (name, count) in _species_counts(reaction.reactants)
            destruction[_species_index(network.species_index, name)] += count * flux
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

function _rhs_at(network::ReactionNetwork, Y::AbstractVector{<:Real}, rho, T9, t::Real; rate_multipliers=nothing, rate_p_values=nothing)
    rho_t = _profile_value(rho, t)
    T9_t = _profile_value(T9, t)
    return network_rhs(Y, network, rho_t, T9_t; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values)
end

function _euler_step(network::ReactionNetwork, Y::Vector{Float64}, t::Float64, dt::Float64, rho, T9; rate_multipliers=nothing, rate_p_values=nothing)
    k1 = _rhs_at(network, Y, rho, T9, t; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values)
    return Y .+ dt .* k1
end

function _rk4_step(network::ReactionNetwork, Y::Vector{Float64}, t::Float64, dt::Float64, rho, T9; rate_multipliers=nothing, rate_p_values=nothing)
    k1 = _rhs_at(network, Y, rho, T9, t; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values)
    k2 = _rhs_at(network, Y .+ 0.5 * dt .* k1, rho, T9, t + 0.5 * dt; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values)
    k3 = _rhs_at(network, Y .+ 0.5 * dt .* k2, rho, T9, t + 0.5 * dt; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values)
    k4 = _rhs_at(network, Y .+ dt .* k3, rho, T9, t + dt; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values)
    return Y .+ (dt / 6.0) .* (k1 .+ 2.0 .* k2 .+ 2.0 .* k3 .+ k4)
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

Supported methods are `:euler` and `:rk4`. RK4 is usually more accurate for the
same timestep, while Euler is useful for simple debugging.

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
)
    times = _time_grid(tspan, dt)
    Y = _checked_initial_abundances(Y0, network)
    Y_history = Matrix{Float64}(undef, length(times), length(Y))
    Y_history[1, :] .= Y

    for step_index in 1:(length(times)-1)
        t = times[step_index]
        step_dt = times[step_index+1] - t

        if method == :euler
            Y = _euler_step(network, Y, t, step_dt, rho, T9; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values)
        elseif method == :rk4
            Y = _rk4_step(network, Y, t, step_dt, rho, T9; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values)
        else
            throw(ArgumentError("unsupported method $method; use :euler or :rk4"))
        end

        if clamp_negative
            for i in eachindex(Y)
                Y[i] < 0.0 && (Y[i] = 0.0)
            end
        end

        Y_history[step_index+1, :] .= Y
    end

    return times, Y_history
end

function _single_step(network::ReactionNetwork, Y::Vector{Float64}, t::Float64, dt::Float64, rho, T9, method::Symbol; rate_multipliers=nothing, rate_p_values=nothing)
    if method == :euler
        return _euler_step(network, Y, t, dt, rho, T9; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values)
    elseif method == :rk4
        return _rk4_step(network, Y, t, dt, rho, T9; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values)
    end
    throw(ArgumentError("unsupported method $method; use :euler or :rk4"))
end

function _max_fractional_change(Y::Vector{Float64}, Y_next::Vector{Float64}, abundance_floor::Float64)
    max_change = 0.0
    for i in eachindex(Y)
        scale = max(abs(Y[i]), abundance_floor)
        max_change = max(max_change, abs(Y_next[i] - Y[i]) / scale)
    end
    return max_change
end

#=
    solve_network_adaptive(network, Y0, tspan, dt_initial, rho, T9; ...)

Evolve a network with simple adaptive explicit timestepping. A proposed step is
accepted when the maximum fractional abundance change is below
`max_fractional_change`, using `abundance_floor` to avoid division by zero for
trace species.

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
)
    t_start = Float64(tspan[1])
    t_end = Float64(tspan[2])
    dt = min(Float64(dt_initial), Float64(dt_max))
    _validate_time_inputs(t_start, t_end, dt)
    max_fractional_change > 0.0 || throw(ArgumentError("max_fractional_change must be positive"))
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

    while t < t_end && accepted_steps < max_steps
        step_dt = min(dt, t_end - t)
        Y_next = _single_step(network, Y, t, step_dt, rho, T9, method; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values)

        if clamp_negative
            for i in eachindex(Y_next)
                Y_next[i] < 0.0 && (Y_next[i] = 0.0)
            end
        end

        change = _max_fractional_change(Y, Y_next, Float64(abundance_floor))
        if change <= max_fractional_change || step_dt <= dt_min
            t += step_dt
            Y = Y_next
            push!(times, t)
            push!(history_rows, copy(Y))
            accepted_steps += 1

            if change == 0.0
                dt = min(step_dt * growth_factor, Float64(dt_max))
            else
                factor = clamp(safety * max_fractional_change / change, shrink_factor, growth_factor)
                dt = min(max(step_dt * factor, Float64(dt_min)), Float64(dt_max))
            end
        else
            dt = max(step_dt * shrink_factor, Float64(dt_min))
            if step_dt <= dt_min
                throw(ArgumentError("adaptive timestep reached dt_min=$dt_min but fractional change=$change exceeds limit=$max_fractional_change"))
            end
        end
    end

    accepted_steps < max_steps || throw(ArgumentError("adaptive solver exceeded max_steps=$max_steps"))

    history = Matrix{Float64}(undef, length(history_rows), length(Y))
    for (i, row) in pairs(history_rows)
        history[i, :] .= row
    end

    return times, history
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

end # module
