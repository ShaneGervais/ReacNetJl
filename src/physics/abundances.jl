# Mass fraction <-> abundance conversion and bulk diagnostics.

function _normalize_mass_fraction_keys(X::AbstractDict)
    normalized = Dict{String,Float64}()
    for (name, value) in X
        normalized_name = normalize_species_name(string(name))
        haskey(normalized, normalized_name) && throw(ArgumentError("duplicate mass fraction entry for species '$normalized_name'"))
        normalized[normalized_name] = Float64(value)
    end
    return normalized
end

"""
    abundances_from_mass_fractions(network, X; normalize=false, check_sum=false, atol=1e-8)

Convert a dictionary of mass fractions `X_i` into an abundance vector `Y_i`
ordered like `network.species`, using

```math
Y_i = \\frac{X_i}{A_i}
```

per species (via `abundance_from_mass_fraction`). Abundance `Y_i` (mol g^-1)
is the quantity the ODE solvers actually evolve, since it turns nuclear
reaction counting into simple linear stoichiometry (`dY_i/dt` sums directly
over reactions, whereas mass fractions would need extra `A_i` factors
everywhere); mass fractions are recovered for reporting via
`mass_fractions_from_abundances`.

Missing network species are assigned zero mass fraction. Extra species not in
the network raise an error. If `normalize=true`, mass fractions are divided by
their total before conversion; if `check_sum=true` (and `normalize=false`), a
total that deviates from 1 by more than `atol` raises an error instead.
"""
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

"""
    mass_fractions_from_abundances(network, Y)

Convert an abundance vector `Y_i` ordered like `network.species` into a
dictionary of mass fractions keyed by normalized species name, using the
inverse of `abundances_from_mass_fractions`:

```math
X_i = A_i Y_i
```
"""
function mass_fractions_from_abundances(network::ReactionNetwork, Y::AbstractVector{<:Real})
    length(Y) == length(network.species) || throw(ArgumentError("Y length must match the number of network species"))

    X = Dict{String,Float64}()
    for (i, species) in pairs(network.species_info)
        species.A > 0 || throw(ArgumentError("cannot convert abundance for species '$(species.name)' with A=$(species.A)"))
        X[species.name] = mass_fraction_from_abundance(Y[i], species.A)
    end
    return X
end

"""
    mass_fraction_history(network, history)

Convert an abundance-history matrix (`history[n, i] = Y_i(t_n)`, one row per
saved solver time) into a mass-fraction-history matrix with the same shape,
applying `X_i = A_i Y_i` row by row. Columns remain ordered like
`network.species`.
"""
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

"""
    first_mass_fraction_threshold_crossing(network, times, history, species, threshold; direction=:down)

Find the first time a species' mass fraction `X(t)` crosses `threshold`,
linearly interpolating between the two bracketing saved states rather than
just returning the nearest grid point. For `direction=:down`, this locates
where `X` first falls to or below `threshold` (e.g. the hydrogen-exhaustion
time, `--stop-hydrogen` in the example driver); `direction=:up` finds the
first rise to or above it instead.

Between two adjacent saved states with mass fractions `x_0` at `t_0` and `x_1`
at `t_1`, the crossing fraction along the segment is

```math
f = \\frac{\\mathrm{threshold} - x_0}{x_1 - x_0}, \\qquad t_{\\mathrm{cross}} = t_0 + f\\,(t_1 - t_0)
```

and the full abundance state is linearly interpolated the same way
(`state = (1-f)\\,Y(t_0) + f\\,Y(t_1)`) so callers can resume/truncate a run
from the crossing point. Returns `nothing` if the threshold is never crossed;
returns immediately with `fraction=0.0` if the very first saved state already
satisfies the crossing condition.
"""
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

"""
    total_mass_fraction(network, Y)

Return the total mass fraction of the active network,

```math
\\sum_i X_i = \\sum_i A_i Y_i
```

for one abundance state. Should stay at (or very near) 1 throughout a run
since baryon number is conserved by every reaction (`network_validation_report`
checks this per-reaction); systematic drift here signals a solver-accuracy or
bookkeeping problem, not real physics.
"""
function total_mass_fraction(network::ReactionNetwork, Y::AbstractVector{<:Real})
    length(Y) == length(network.species) || throw(ArgumentError("Y length must match the number of network species"))

    total = 0.0
    for (i, species) in pairs(network.species_info)
        species.A > 0 || throw(ArgumentError("cannot convert abundance for species '$(species.name)' with A=$(species.A)"))
        total += mass_fraction_from_abundance(Y[i], species.A)
    end
    return total
end

"""
    total_mass_fraction_history(network, history)

Return `total_mass_fraction(network, Y(t_n))` at every saved abundance-history
row -- the raw data behind `mass_fraction_drift`'s conservation summary.
"""
function total_mass_fraction_history(network::ReactionNetwork, history::AbstractMatrix{<:Real})
    size(history, 2) == length(network.species) || throw(ArgumentError("history column count must match the number of network species"))

    totals = Vector{Float64}(undef, size(history, 1))
    for i in axes(history, 1)
        totals[i] = total_mass_fraction(network, view(history, i, :))
    end
    return totals
end

"""
    mass_fraction_drift(network, history)

Summarize total-mass-fraction (`\\sum_i X_i`, see `total_mass_fraction`) drift
over a history: the initial and final totals, net `drift = final - initial`,
the largest absolute deviation from the initial total at any saved time
(`max_abs_drift`), and the min/max totals reached. This is the network's
baryon-conservation health check -- a well-behaved solve should show `drift`
and `max_abs_drift` many orders of magnitude below 1.
"""
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

"""
    abundance_diagnostics(network, history)

Report positivity-oriented diagnostics for an abundance history: the minimum
abundance `Y_i` and minimum mass fraction `X_i` reached anywhere in the run,
each with the species and time-index where it occurred, plus boolean flags
for whether either dipped negative. Explicit and adaptive backward-Euler
solving can produce small negative abundances near a true zero from
finite-precision Newton iteration; `has_negative_abundance` /
`has_negative_mass_fraction` distinguish "solver noise at the floor" from a
real instability without requiring a full history scan by hand.
"""
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

