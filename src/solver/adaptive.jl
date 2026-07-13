# Adaptive step-size controller.


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

