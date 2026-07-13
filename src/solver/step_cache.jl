# Per-step rate/screening cache and the cached RHS/Jacobian the solvers share.

#=
Per-step solver cache.

Within one solver step the temperature and density are fixed, so the
interpolated rate, the density power, the symmetry factor, and any rate
multipliers are constant. `NetworkStepCache` hoists all of that out of the
abundance loops: an RHS or Jacobian evaluation then only multiplies cached
prefactors by abundance products. This changes no mathematics — only where
the temperature-dependent work is done.

Custom screening functions cannot be decomposed this way, so
`_build_step_cache` returns `nothing` for them and callers fall back to the
uncached path.
=#
struct NetworkStepCache
    rho::Float64
    T9::Float64
    prefactors::Vector{Float64}
    screening::Symbol
    screening_scale::Float64
end

function _build_step_cache(
    network::ReactionNetwork,
    rho_t::Real,
    T9_t::Real;
    rate_multipliers=nothing,
    rate_p_values=nothing,
    screening=nothing,
)
    screening isa Function && return nothing
    mode = screening === nothing || screening === false ? :none :
           screening == :weak ? :weak :
           screening == :chugunov ? :chugunov :
           throw(ArgumentError("unsupported screening=$screening; use nothing, :weak, :chugunov, or a function `(network, reaction, Y, rho, T9) -> multiplier`"))

    nreactions = length(network.reactions)
    if rate_multipliers !== nothing && length(rate_multipliers) != nreactions
        throw(ArgumentError("rate_multipliers must have the same length as reactions"))
    end
    if rate_p_values !== nothing && length(rate_p_values) != nreactions
        throw(ArgumentError("rate_p_values must have the same length as reactions"))
    end

    rho_value = Float64(rho_t)
    T9_value = Float64(T9_t)
    prefactors = Vector{Float64}(undef, nreactions)
    for r in 1:nreactions
        table = network.reactions[r].rate_table
        compiled = network.compiled_reactions[r]
        p_value = rate_p_values === nothing ? nothing : rate_p_values[r]
        base_rate = p_value === nothing ? interpolate_rate(table, T9_value) : sampled_interpolate_rate(table, T9_value, p_value)
        multiplier = rate_multipliers === nothing ? 1.0 : Float64(rate_multipliers[r])
        prefactors[r] = multiplier * base_rate * rho_value^(compiled.nreactants - 1) / compiled.symmetry_factor
    end

    screening_scale = 0.0
    if mode != :none
        T9_value > 0.0 || throw(ArgumentError("T9 must be positive for screening"))
        rho_value > 0.0 || throw(ArgumentError("rho must be positive for screening"))
        if mode == :weak
            T6 = 1000.0 * T9_value
            screening_scale = 0.188 * sqrt(rho_value / T6^3)
        end
    end

    return NetworkStepCache(rho_value, T9_value, prefactors, mode, screening_scale)
end

function _screening_zeta_scale(cache::NetworkStepCache, network::ReactionNetwork, Y::AbstractVector{Float64})
    cache.screening == :weak || return 0.0
    zeta = _screening_composition_factor(network, Y)
    zeta > 0.0 || return 0.0
    return cache.screening_scale * sqrt(zeta)
end

@inline function _cached_screening_multiplier(compiled::CompiledReaction, zeta_scale::Float64)
    (zeta_scale > 0.0 && compiled.charge_pair_sum > 0.0) || return 1.0
    return exp(min(compiled.charge_pair_sum * zeta_scale, _SCREENING_MAX_EXPONENT))
end

# Per-evaluation screening context: the composition-dependent scalars shared
# by every reaction of one RHS/Jacobian evaluation. A concrete struct keeps
# the per-reaction screening call dispatch-free and allocation-free.
struct ScreeningContext
    zeta_scale::Float64
    n_e::Float64
    gamma_e_fac::Float64
    plasma_active::Bool
end

function _screening_context(cache::NetworkStepCache, network::ReactionNetwork, Y::AbstractVector{Float64})
    if cache.screening == :weak
        return ScreeningContext(_screening_zeta_scale(cache, network, Y), 0.0, 0.0, false)
    elseif cache.screening == :chugunov
        plasma = _ion_plasma_state(network, Y, cache.rho)
        plasma === nothing && return ScreeningContext(0.0, 0.0, 0.0, false)
        return ScreeningContext(0.0, plasma.n_e, plasma.gamma_e_fac, true)
    end
    return ScreeningContext(0.0, 0.0, 0.0, false)
end

@inline function _cached_reaction_screening(cache::NetworkStepCache, compiled::CompiledReaction, context::ScreeningContext)
    if cache.screening == :weak
        return _cached_screening_multiplier(compiled, context.zeta_scale)
    elseif cache.screening == :chugunov
        return _chugunov_reaction_multiplier(compiled, context, cache.T9)
    end
    return 1.0
end

# In-place dY/dt with all temperature-dependent factors taken from the cache.
function _cached_network_rhs!(dYdt::Vector{Float64}, network::ReactionNetwork, cache::NetworkStepCache, Y::AbstractVector{Float64})
    fill!(dYdt, 0.0)
    screening_context = _screening_context(cache, network, Y)

    @inbounds for r in eachindex(network.compiled_reactions)
        compiled = network.compiled_reactions[r]
        flux = cache.prefactors[r] * _cached_reaction_screening(cache, compiled, screening_context)
        for index in compiled.reactant_indices
            flux *= Y[index]
        end
        flux == 0.0 && continue

        for (index, count) in zip(compiled.reactant_species_indices, compiled.reactant_species_counts)
            dYdt[index] -= count * flux
        end
        for (index, count) in zip(compiled.product_species_indices, compiled.product_species_counts)
            dYdt[index] += count * flux
        end
    end

    return dYdt
end

#=
Analytic Jacobian d(dY/dt)/dY. The flux of each reaction is a polynomial in
the reactant abundances, so the derivative follows from the product rule over
the distinct reactant species. The weak-screening multiplier is treated as
constant with respect to Y; Newton's converged answer is fixed by the exact
residual alone, so this only shapes the iteration path, not the solution.
=#
function _cached_network_jacobian!(J::Matrix{Float64}, network::ReactionNetwork, cache::NetworkStepCache, Y::AbstractVector{Float64})
    fill!(J, 0.0)
    screening_context = _screening_context(cache, network, Y)

    @inbounds for r in eachindex(network.compiled_reactions)
        compiled = network.compiled_reactions[r]
        base = cache.prefactors[r] * _cached_reaction_screening(cache, compiled, screening_context)
        base == 0.0 && continue

        reactant_indices = compiled.reactant_species_indices
        reactant_counts = compiled.reactant_species_counts
        for jpos in eachindex(reactant_indices)
            jindex = reactant_indices[jpos]
            count_j = reactant_counts[jpos]
            derivative = Float64(count_j)
            count_j > 1 && (derivative *= Y[jindex]^(count_j - 1))
            for kpos in eachindex(reactant_indices)
                kpos == jpos && continue
                derivative *= Y[reactant_indices[kpos]]^reactant_counts[kpos]
            end
            dflux = base * derivative
            dflux == 0.0 && continue

            for (index, count) in zip(reactant_indices, reactant_counts)
                J[index, jindex] -= count * dflux
            end
            for (index, count) in zip(compiled.product_species_indices, compiled.product_species_counts)
                J[index, jindex] += count * dflux
            end
        end
    end

    return J
end


_profile_value(value::Real, t::Real) = Float64(value)
_profile_value(value, t::Real) = Float64(value(t))

function _rhs_at(network::ReactionNetwork, Y::AbstractVector{<:Real}, rho, T9, t::Real; rate_multipliers=nothing, rate_p_values=nothing, screening=nothing)
    rho_t = _profile_value(rho, t)
    T9_t = _profile_value(T9, t)
    return network_rhs(Y, network, rho_t, T9_t; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)
end

function _step_cache_at(network::ReactionNetwork, rho, T9, t::Float64; rate_multipliers=nothing, rate_p_values=nothing, screening=nothing)
    return _build_step_cache(
        network,
        _profile_value(rho, t),
        _profile_value(T9, t);
        rate_multipliers=rate_multipliers,
        rate_p_values=rate_p_values,
        screening=screening,
    )
end
