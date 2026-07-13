# Explicit Euler and RK4 steppers.


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
