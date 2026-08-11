# Implicit backward-Euler stepper with Newton iteration.
#
# Backward Euler solves the *implicit* update
#
#   Y(t+dt) = Y(t) + dt * f(Y(t+dt), t+dt)
#
# (contrast forward Euler's explicit Y(t)+dt*f(Y(t),t), see explicit.jl) by
# finding the root of the residual
#
#   R(Y_next) = Y_next - Y - dt * f(Y_next, t+dt) = 0
#
# via Newton's method: at each iteration, solve the linear system
# `J * correction = -R` for the Newton matrix `J = I - dt * df/dY` (the
# identity from differentiating Y_next w.r.t. itself, minus dt times the
# network's own Jacobian df/dY), then update `Y_next += correction`. This is
# unconditionally stable for stiff systems (unlike explicit methods), which
# is why it's the default integrator for real nova reaction networks: rate
# multiplier differences of many orders of magnitude between fast
# (e.g. proton captures) and slow (e.g. some weak decays) channels make the
# network numerically stiff.

# Newton corrections can occasionally overshoot into negative abundances
# (unphysical); rather than reject the whole step, shrink the step along the
# correction direction (a "damped Newton" line search) just enough that no
# abundance goes more than 10% negative, matching how far a physical
# solution could plausibly move in one Newton iteration.
function _positivity_limited_alpha(Y::Vector{Float64}, correction::Vector{Float64})
    alpha = 1.0
    for i in eachindex(Y)
        correction[i] < 0.0 || continue
        Y[i] > 0.0 || continue
        alpha = min(alpha, 0.9 * Y[i] / -correction[i])
    end
    return clamp(alpha, 0.0, 1.0)
end

# Clamp abundances that went negative only by a tiny amount (within `floor`
# of zero) up to exactly zero. This is finite-precision Newton-iteration
# noise around a true zero abundance, not a real instability; leaving it
# slightly negative would otherwise propagate into log-scale diagnostics
# (mass fractions, flux ratios) as spurious NaN/Inf.
function _clamp_tiny_negative_trials!(Y::Vector{Float64}, floor::Float64)
    for i in eachindex(Y)
        if Y[i] < 0.0 && abs(Y[i]) <= floor
            Y[i] = 0.0
        end
    end
    return Y
end

struct NewtonConvergenceError <: Exception
    iterations::Int
    residual_norm::Float64
end

function Base.showerror(io::IO, err::NewtonConvergenceError)
    print(io, "backward Euler Newton iteration failed to converge after $(err.iterations) iterations; residual norm = $(err.residual_norm)")
end

# R(Y_next) = Y_next - Y - dt*f(Y_next, t_next), the backward-Euler residual
# whose root Newton's method seeks (see the file-level comment above).
function _backward_euler_residual(network::ReactionNetwork, Y_next::Vector{Float64}, Y::Vector{Float64}, t_next::Float64, dt::Float64, rho, T9; rate_multipliers=nothing, rate_p_values=nothing, screening=nothing)
    return Y_next .- Y .- dt .* _rhs_at(network, Y_next, rho, T9, t_next; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)
end

# Finite-difference Jacobian dR/dY_next of the backward-Euler residual, one
# perturbed column at a time (in parallel across Julia threads): column j is
# `(R(Y_next + eps*e_j) - R(Y_next)) / eps`. Retained mainly for validating
# `jacobian=:analytic` (the default and much faster path, an ~360x speedup
# on the 144-species nova network per `examples/benchmark_solver.jl`, since
# it avoids one extra full network RHS evaluation per species).
function _backward_euler_jacobian(network::ReactionNetwork, Y_next::Vector{Float64}, residual::Vector{Float64}, Y::Vector{Float64}, t_next::Float64, dt::Float64, rho, T9, finite_difference_epsilon::Float64; rate_multipliers=nothing, rate_p_values=nothing, screening=nothing)
    n = length(Y_next)
    jacobian = Matrix{Float64}(undef, n, n)

    Base.Threads.@threads for j in 1:n
        perturbed = copy(Y_next)
        saved = Y_next[j]
        step = finite_difference_epsilon * max(abs(saved), 1.0)
        perturbed[j] = saved + step
        perturbed_residual = _backward_euler_residual(network, perturbed, Y, t_next, dt, rho, T9; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)
        jacobian[:, j] .= (perturbed_residual .- residual) ./ step
    end

    return jacobian
end

# One implicit backward-Euler timestep from (Y, t) to (Y_next, t+dt), solved
# by damped Newton iteration (see the file-level comment above). Starts the
# Newton iteration from an explicit-Euler predictor (`_euler_step`) rather
# than `Y` itself, giving the iteration a head start toward the true
# solution. `jacobian=:analytic` (default) uses the network's precomputed
# analytic Jacobian (`_cached_network_jacobian!`) forming `I - dt*J`
# directly; `:finite_difference` falls back to `_backward_euler_jacobian`.
# Each iteration's correction is damped by `_positivity_limited_alpha` and
# backtracked (halving `alpha` up to 12 times) until the residual norm
# actually decreases, so a single bad Newton step can't blow up the
# iteration. Throws `NewtonConvergenceError` if `max_newton_iterations` is
# exhausted without meeting `newton_tolerance` (relative to `max(|Y_next|,
# 1)`) -- the adaptive solver (`solve_network_adaptive`) catches this and
# retries with a smaller `dt` rather than treating it as fatal.
function _backward_euler_step(
    network::ReactionNetwork,
    Y::Vector{Float64},
    t::Float64,
    dt::Float64,
    rho,
    T9;
    rate_multipliers=nothing,
    rate_p_values=nothing,
    screening=nothing,
    newton_tolerance::Real=1.0e-10,
    max_newton_iterations::Integer=20,
    finite_difference_epsilon::Real=sqrt(eps(Float64)),
    jacobian::Symbol=:analytic,
    newton_iterations=nothing,
)
    max_newton_iterations > 0 || throw(ArgumentError("max_newton_iterations must be positive"))
    newton_tolerance > 0.0 || throw(ArgumentError("newton_tolerance must be positive"))
    finite_difference_epsilon > 0.0 || throw(ArgumentError("finite_difference_epsilon must be positive"))
    jacobian in (:analytic, :finite_difference, :sparse) || throw(ArgumentError("unsupported jacobian=$jacobian; use :analytic, :finite_difference, or :sparse"))

    t_next = t + dt
    Y_next = _euler_step(network, Y, t, dt, rho, T9; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)

    # All Newton residuals of this step share the fixed (rho, T9) at t_next.
    cache = _step_cache_at(network, rho, T9, t_next; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)
    n = length(Y)
    rhs_buffer = cache === nothing ? nothing : Vector{Float64}(undef, n)

    residual_at = function (Y_trial::Vector{Float64})
        if cache === nothing
            return _backward_euler_residual(network, Y_trial, Y, t_next, dt, rho, T9; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)
        end
        _cached_network_rhs!(rhs_buffer, network, cache, Y_trial)
        return Y_trial .- Y .- dt .* rhs_buffer
    end

    residual = residual_at(Y_next)
    tolerance = Float64(newton_tolerance) * max(_max_abs(Y_next), 1.0)
    if _max_abs(residual) <= tolerance
        newton_iterations !== nothing && push!(newton_iterations, 0)
        return Y_next
    end

    use_analytic_jacobian = cache !== nothing && jacobian == :analytic
    use_sparse_jacobian = cache !== nothing && jacobian == :sparse
    jacobian_buffer = use_analytic_jacobian ? Matrix{Float64}(undef, n, n) : nothing
    sparse_jacobian_buffer = use_sparse_jacobian ? sparse_jacobian_prototype(network) : nothing

    for iteration in 1:max_newton_iterations
        residual_norm = _max_abs(residual)
        if use_analytic_jacobian
            # Newton matrix of the residual: I - dt * d(dY/dt)/dY.
            _cached_network_jacobian!(jacobian_buffer, network, cache, Y_next)
            @inbounds for j in 1:n
                for i in 1:n
                    jacobian_buffer[i, j] = -dt * jacobian_buffer[i, j]
                end
                jacobian_buffer[j, j] += 1.0
            end
            jacobian_matrix = jacobian_buffer
        elseif use_sparse_jacobian
            # Same Newton matrix as the dense analytic path, assembled at the
            # network's precomputed sparsity structure instead (see
            # `_cached_network_jacobian_sparse!`); `colptr`/`rowval` never
            # change here, only `nzval`.
            _cached_network_jacobian_sparse!(sparse_jacobian_buffer, network, cache, Y_next)
            sparse_jacobian_buffer.nzval .*= -dt
            @inbounds for j in 1:n
                sparse_jacobian_buffer.nzval[_sparse_nzval_index(sparse_jacobian_buffer, j, j)] += 1.0
            end
            jacobian_matrix = sparse_jacobian_buffer
        else
            jacobian_matrix = _backward_euler_jacobian(
                network,
                Y_next,
                residual,
                Y,
                t_next,
                dt,
                rho,
                T9,
                Float64(finite_difference_epsilon);
                rate_multipliers=rate_multipliers,
                rate_p_values=rate_p_values,
                screening=screening,
            )
        end
        # A singular Newton matrix (e.g. dt so large that I - dt*J loses the
        # identity part to rounding along a conserved direction) is a step
        # failure, not a fatal error: report non-convergence so adaptive
        # drivers shrink dt and retry. The sparse path routes through
        # `_sparse_newton_solve` (KLU when the extension is loaded) rather
        # than SparseArrays' own default `\` (UMFPACK), so which sparse
        # backend actually solved the system stays explicit rather than an
        # implicit fallback.
        correction = try
            use_sparse_jacobian ? _sparse_newton_solve(jacobian_matrix, -residual) : jacobian_matrix \ (-residual)
        catch err
            if err isa LinearAlgebra.SingularException || err isa LinearAlgebra.LAPACKException
                throw(NewtonConvergenceError(iteration, residual_norm))
            end
            rethrow()
        end

        _all_finite(correction) || throw(NewtonConvergenceError(iteration, residual_norm))

        accepted_trial = false
        best_trial = copy(Y_next)
        best_residual = copy(residual)
        best_residual_norm = residual_norm
        alpha = _positivity_limited_alpha(Y_next, correction)
        alpha = alpha > 0.0 ? alpha : 1.0

        for _ in 1:12
            trial = Y_next .+ alpha .* correction
            _clamp_tiny_negative_trials!(trial, 1.0e-30)
            if _all_finite(trial) && minimum(trial; init=0.0) >= 0.0
                trial_residual = residual_at(trial)
                trial_residual_norm = _max_abs(trial_residual)
                if trial_residual_norm < best_residual_norm
                    best_trial = trial
                    best_residual = trial_residual
                    best_residual_norm = trial_residual_norm
                end
                if trial_residual_norm < residual_norm || trial_residual_norm <= tolerance
                    Y_next = trial
                    residual = trial_residual
                    accepted_trial = true
                    break
                end
            end
            alpha *= 0.5
        end

        if !accepted_trial
            if best_residual_norm < residual_norm
                Y_next = best_trial
                residual = best_residual
            else
                Y_next .+= correction
                residual = residual_at(Y_next)
            end
        end

        tolerance = Float64(newton_tolerance) * max(_max_abs(Y_next), 1.0)
        if _max_abs(residual) <= tolerance || _max_abs(correction) <= tolerance
            newton_iterations !== nothing && push!(newton_iterations, iteration)
            return Y_next
        end
    end

    throw(NewtonConvergenceError(Int(max_newton_iterations), _max_abs(residual)))
end
