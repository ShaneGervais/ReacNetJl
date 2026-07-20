# Fixed-step time integration driver.


# Summarize backward-Euler Newton-iteration counts across a run (mean/max
# iterations per step, and how many steps failed to converge) for the
# `solver_stats` diagnostic returned alongside the abundance history.
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

# Assemble the `solver_stats` named tuple for a fixed-step run. Since every
# step uses the same `dt` by construction, `rejected_steps`/
# `max_fractional_change`/`reached_dt_min` are trivial placeholders that
# only carry real information for the adaptive solver (see adaptive.jl).
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


function _validate_time_inputs(t_start::Float64, t_end::Float64, dt::Float64)
    t_end > t_start || throw(ArgumentError("tspan must have t_end > t_start"))
    dt > 0 || throw(ArgumentError("dt must be positive"))
end

# Build the fixed-step time grid from t_start to t_end in steps of dt. The
# final grid point is snapped exactly to t_end (extending with a final short
# step if dt doesn't evenly divide the interval, or trimming an overshoot
# back to t_end) so the run always finishes exactly at the requested time
# rather than overshooting or falling short by a fractional step.
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
"""
    solve_network(network, Y0, tspan, dt, rho, T9; method=:rk4, rate_multipliers=nothing, clamp_negative=true)

Evolve a single-zone reaction network in time with a *fixed* timestep `dt`
(see `solve_network_adaptive` for automatic step-size control, the choice
`run_ppn` actually uses).

This solves the ordinary differential equation system `dY/dt = f(Y, rho, T9)`
(`f` = `network_rhs`). For a single-zone post-processing network there are no
spatial derivatives, so this is an ODE system rather than a PDE.

Arguments:
- `network`: a `ReactionNetwork`.
- `Y0`: initial abundance vector ordered like `network.species`.
- `tspan`: `(t_start, t_end)` in seconds.
- `dt`: fixed timestep in seconds.
- `rho`: density in g cm^-3, or a function `rho(t)` (e.g. from `trajectory_profiles`).
- `T9`: temperature in GK, or a function `T9(t)`.

Supported methods are `:euler`, `:rk4`, and `:backward_euler` (see
explicit.jl/backward_euler.jl for each method's update formula). RK4 is
usually more accurate for the same timestep among the explicit choices, while
backward Euler is the implicit, unconditionally-stable option needed for
stiff nova networks. With `clamp_negative=true` (default), any abundance that
goes slightly negative after a step is clamped to zero.

Returns `(times, Y_history)`, where `Y_history[i, j]` is the abundance of
species `network.species[j]` at `times[i]`; pass `return_stats=true` to also
get a `solver_stats` named tuple (Newton-iteration counts, if
`method=:backward_euler`).
"""
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

# Dispatch one timestep to the requested method's stepper; shared by
# `solve_network`'s fixed-step loop and `solve_network_adaptive`'s
# accept/reject loop (adaptive.jl), so both drivers stay in sync about which
# methods exist and how each is invoked.
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
