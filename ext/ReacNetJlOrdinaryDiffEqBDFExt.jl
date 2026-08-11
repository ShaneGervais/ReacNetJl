# Package extension implementing `solve_network_fbdf` (see src/solver/fbdf.jl)
# with OrdinaryDiffEqBDF's variable-order stiff FBDF integrator, reusing
# ReacNetJl's own cached RHS and analytic Jacobian (step_cache.jl) rather than
# a finite-difference or AD-detected Jacobian.
module ReacNetJlOrdinaryDiffEqBDFExt

using ReacNetJl
using OrdinaryDiffEqBDF

import ReacNetJl: ReactionNetwork, _profile_value, _step_cache_at, _cached_network_rhs!, _cached_network_jacobian!, _checked_initial_abundances

# Bundles everything the RHS/Jacobian callbacks need as the ODEProblem's `p`,
# since SciML calls back as `f!(du, u, p, t)`/`jac!(J, u, p, t)` with no other
# hook for closed-over state. `rho`/`T9` may be a constant `Real` or a
# callable profile (`_profile_value` resolves either at a given `t`, matching
# every other solver in this package).
struct FBDFParams{N,R,Tm,RM,RP}
    network::N
    rho::R
    T9::Tm
    rate_multipliers::RM
    rate_p_values::RP
    screening::Union{Nothing,Symbol}
end

function _fbdf_rhs!(du, u, p::FBDFParams, t)
    rho_t = _profile_value(p.rho, t)
    T9_t = _profile_value(p.T9, t)
    cache = _step_cache_at(p.network, rho_t, T9_t, t; rate_multipliers=p.rate_multipliers, rate_p_values=p.rate_p_values, screening=p.screening)
    _cached_network_rhs!(du, p.network, cache, u)
    return nothing
end

function _fbdf_jac!(J, u, p::FBDFParams, t)
    rho_t = _profile_value(p.rho, t)
    T9_t = _profile_value(p.T9, t)
    cache = _step_cache_at(p.network, rho_t, T9_t, t; rate_multipliers=p.rate_multipliers, rate_p_values=p.rate_p_values, screening=p.screening)
    _cached_network_jacobian!(J, p.network, cache, u)
    return nothing
end

function ReacNetJl.solve_network_fbdf(
    network::ReactionNetwork,
    Y0,
    tspan::Tuple{<:Real,<:Real},
    rho,
    T9;
    rate_multipliers=nothing,
    rate_p_values=nothing,
    screening=nothing,
    reltol::Real=1.0e-8,
    abstol::Real=1.0e-10,
)
    # Custom screening functions can't be decomposed into the per-step cache
    # that both the RHS and the analytic Jacobian rely on here (see
    # `_build_step_cache`'s `screening isa Function` short-circuit); every
    # other solver in this package silently falls back to a slower uncached
    # path for that case, but FBDF's `jac_prototype` is fixed at problem
    # construction, so there is no equivalent fallback -- reject it plainly
    # instead of quietly ignoring the custom screening or guessing a
    # Jacobian for it.
    screening === nothing || screening isa Symbol && screening in (:weak, :chugunov) ||
        throw(ArgumentError("solve_network_fbdf only supports screening=nothing, :weak, or :chugunov (a custom screening function cannot be cached for the SciML Jacobian interface); got $(repr(screening))"))

    t_start = Float64(tspan[1])
    t_end = Float64(tspan[2])
    t_end > t_start || throw(ArgumentError("tspan must have t_end > t_start"))

    Y_checked = _checked_initial_abundances(Y0, network)
    n = length(Y_checked)
    params = FBDFParams(network, rho, T9, rate_multipliers, rate_p_values, screening)

    ode_function = ODEFunction(_fbdf_rhs!; jac=_fbdf_jac!, jac_prototype=Matrix{Float64}(undef, n, n))
    problem = ODEProblem(ode_function, Y_checked, (t_start, t_end), params)
    sol = solve(problem, FBDF(); reltol=reltol, abstol=abstol)

    sol.retcode == ReturnCode.Success || error("solve_network_fbdf: integration failed with retcode=$(sol.retcode)")

    times = sol.t
    history = Matrix{Float64}(undef, length(times), n)
    for (i, u) in pairs(sol.u)
        history[i, :] .= u
    end

    step_sizes = diff(times)
    stats = (
        accepted_steps=sol.stats.naccept,
        rejected_steps=sol.stats.nreject,
        min_dt=isempty(step_sizes) ? NaN : minimum(step_sizes),
        max_dt=isempty(step_sizes) ? NaN : maximum(step_sizes),
        final_dt=isempty(step_sizes) ? NaN : last(step_sizes),
        # not meaningful for this solver (no fixed fractional/absolute change
        # controller, no dt_min floor); kept for shape parity with the other
        # solvers' solver_stats rather than silently omitted.
        max_fractional_change=NaN,
        max_absolute_change=NaN,
        reached_dt_min=false,
        function_evaluations=sol.stats.nf,
        jacobian_evaluations=sol.stats.njacs,
        retcode=sol.retcode,
    )

    return (times, history, stats)
end

end # module
