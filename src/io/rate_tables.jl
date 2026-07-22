# Generic ReactionRateTable lookup and interpolation, independent of data source.

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

# Hashable key identifying a reaction by its reactant/product species alone
# (ignoring rate source, Q-value, etc.) -- used throughout to detect "is this
# the same reaction" when merging/de-duplicating tables from different
# libraries (`_unique_reaction_tables`, `add_reverse_reaction_tables`,
# `iliadis2002_rate_tables`'s paper-table override matching).
#
# Sorted so reactant/product *order* doesn't matter: different rate
# libraries (or REACLIB label sets vs. hand-tabulated paper overrides) don't
# consistently list a reaction's participants in the same order (e.g.
# `["p", "na22"]` vs `["na22", "p"]`), and an order-sensitive key would treat
# those as two different reactions -- silently double-counting the same
# physical process once per ordering instead of overriding/deduplicating it.
function _reaction_participant_key(table::ReactionRateTable)
    return (Tuple(sort(table.reactants)), Tuple(sort(table.products)))
end

# Keep only the first table for each distinct reaction (by
# `_reaction_participant_key`), preserving input order. Used when building a
# network from a table list that may contain more than one source/label for
# the same reaction; `network_from_tables` wants exactly one rate per
# reaction, so this is the arbitrary-but-deterministic tie-break.
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

# Index every table by its `_reaction_participant_key`, for O(1) reverse-rate
# lookup in `add_reverse_reaction_tables` (is there an explicit table whose
# reactants/products are this forward reaction's products/reactants,
# reversed?). Tables that fail baryon/charge bookkeeping are excluded so a
# malformed row can't be mistaken for a real reverse rate.
function _reverse_table_lookup(tables::AbstractVector{ReactionRateTable})
    lookup = Dict{Any,ReactionRateTable}()
    for table in tables
        _valid_nuclear_bookkeeping(table) || continue
        key = _reaction_participant_key(table)
        haskey(lookup, key) || (lookup[key] = table)
    end
    return lookup
end

# Shared log-log interpolation for both `interpolate_rate` and
# `interpolate_factor_uncertainty`: linear in log(value) versus log(T9)
# ( ln y = (1-w) ln y0 + w ln y1, with w = (ln T9 - ln T0)/(ln T1 - ln T0) ),
# the standard choice for thermonuclear rates since both T9 grids and rate
# magnitudes span many orders of magnitude (a linear-linear interpolant would
# badly misrepresent the curvature between grid points). Values are floored
# at `LOG_INTERPOLATION_FLOOR` before taking the log so a tabulated zero
# (astrophysically negligible, but not loggable) doesn't error. Throws if
# `T9` falls outside the table's grid -- this is intentional: extrapolating a
# tabulated rate silently would be worse than failing loudly.
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

"""
    interpolate_rate(table, T9)

Interpolate a reaction rate at temperature `T9` (GK), linearly in
`log(rate)` versus `log(T9)`:

```math
\\ln R(T_9) = (1-w)\\ln R_0 + w\\ln R_1, \\qquad w = \\frac{\\ln T_9 - \\ln T_{9,0}}{\\ln T_{9,1} - \\ln T_{9,0}}
```

between the two grid points `T_{9,0} \\le T_9 \\le T_{9,1}` bracketing the
query. This is usually more physically reasonable than linear-linear
interpolation, since thermonuclear rates and the STARLIB/REACLIB temperature
grid both vary over many orders of magnitude. Throws if `T9` is outside the
table's tabulated range.
"""
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
const GENERATED_REVERSE_RATE_FLOOR = 1.0e-300
const LOG_INTERPOLATION_FLOOR = 1.0e-300

