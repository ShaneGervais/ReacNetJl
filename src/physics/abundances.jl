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

