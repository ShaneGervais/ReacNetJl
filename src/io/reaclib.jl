# JINA REACLIB analytic-fit reader and evaluator.

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

Evaluate the standard REACLIB analytic rate parameterization,

```math
R(T_9) = \\exp\\left(a_1 + \\frac{a_2}{T_9} + \\frac{a_3}{T_9^{1/3}} + a_4 T_9^{1/3} + a_5 T_9 + a_6 T_9^{5/3} + a_7 \\ln T_9\\right)
```

for one 7-coefficient fit `set`, or the sum over a vector of sets belonging to
one reaction (a REACLIB rate is often the sum of a non-resonant term plus one
or more narrow-resonance terms, each its own 7-coefficient set with the same
reaction/label but different `resonance` flags). The exponent is capped at
`REACLIB_EXPONENT_CAP` so evaluating a fit far outside its validated `T9`
range gives a large finite number instead of `Inf`/`NaN`.
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

# Build one `ReactionRateTable` from the fit sets of a single reaction/label
# group, by evaluating `reaclib_rate` (summed over the group) at every point
# of `T9_grid` and flooring at `LOG_INTERPOLATION_FLOOR` (rates are
# interpolated in log-space downstream, so a literal zero would break that).
# The source string follows the STARLIB convention of label plus flag
# suffixes, e.g. `"nacr"`, `"wc12w"`, `"nacrv"`, so weak-rate and source
# handling downstream behaves identically for STARLIB and REACLIB tables.
# REACLIB carries no rate uncertainties, so the factor uncertainty is 1.
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

