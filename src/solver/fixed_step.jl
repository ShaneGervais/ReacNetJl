# Fixed-step time integration driver.


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
