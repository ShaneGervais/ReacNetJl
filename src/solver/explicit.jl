# Explicit Euler and RK4 steppers.

# Forward (explicit) Euler: Y(t+dt) = Y(t) + dt * f(Y(t), t), with
# f = network_rhs. Cheapest and least accurate of the fixed-step methods
# (local error O(dt^2), global error O(dt)); useful mainly as a baseline or
# for very small/non-stiff test networks, since nuclear reaction networks are
# usually too stiff for explicit methods to take practical timesteps.
# `_step_cache_at`/`_cached_network_rhs!` avoid recomputing rates when
# density/temperature/rate multipliers haven't changed since the last
# evaluation at this (t, rho, T9); see step_cache.jl.
function _euler_step(network::ReactionNetwork, Y::Vector{Float64}, t::Float64, dt::Float64, rho, T9; rate_multipliers=nothing, rate_p_values=nothing, screening=nothing)
    cache = _step_cache_at(network, rho, T9, t; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)
    if cache === nothing
        k1 = _rhs_at(network, Y, rho, T9, t; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)
    else
        k1 = _cached_network_rhs!(similar(Y), network, cache, Y)
    end
    return Y .+ dt .* k1
end

# Classical 4th-order Runge-Kutta:
#
#   k1 = f(Y, t)                  k2 = f(Y + dt/2 k1, t + dt/2)
#   k3 = f(Y + dt/2 k2, t + dt/2) k4 = f(Y + dt k3,   t + dt)
#   Y(t+dt) = Y(t) + dt/6 (k1 + 2 k2 + 2 k3 + k4)
#
# Local error O(dt^5), global error O(dt^4) -- much better accuracy per step
# than Euler, at 4x the RHS evaluations. Still an explicit method, so it
# shares Euler's stability limitations on stiff networks; `solve_network`'s
# default `:backward_euler` is the production choice for real nova networks,
# with `:rk4`/`:euler` mainly for cross-checking or mild test problems.
# Density/temperature are re-cached at each of the three distinct evaluation
# times (start, midpoint x2, end) since T9(t)/rho(t) generally vary across
# the substeps for a trajectory-driven run.
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

# Largest-magnitude entry of a vector, used by the adaptive step controller
# to measure how large a proposed step's change was; `init=0.0` makes this
# well-defined (0) for an empty vector instead of erroring.
function _max_abs(values::AbstractVector{<:Real})
    maximum(abs, values; init=0.0)
end

# Did a step produce a finite result? Used to reject steps that overflowed
# to Inf/NaN (e.g. from too large a `dt` at high reaction rates) before they
# corrupt the abundance history.
function _all_finite(values::AbstractVector{<:Real})
    all(isfinite, values)
end
