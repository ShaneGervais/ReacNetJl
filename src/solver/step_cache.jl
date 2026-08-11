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

`prefactors` is parametrized over `T` (not hardcoded `Float64`) so it can hold
AD duals when `rate_multipliers`/`rate_p_values` are themselves being
differentiated (see `_rhs_eltype` in flux.jl and Tier 3 of
`reacnetjl_feature_spec.md`); `rho`/`T9`/`screening_scale` stay plain `Float64`
since they are fixed evaluation points here, not differentiation variables.
=#
struct NetworkStepCache{T}
    rho::Float64
    T9::Float64
    prefactors::Vector{T}
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
    prefactor_type = promote_type(
        rate_multipliers === nothing ? Float64 : eltype(rate_multipliers),
        rate_p_values === nothing ? Float64 : eltype(rate_p_values),
    )
    prefactors = Vector{prefactor_type}(undef, nreactions)
    for r in 1:nreactions
        table = network.reactions[r].rate_table
        compiled = network.compiled_reactions[r]
        p_value = rate_p_values === nothing ? nothing : rate_p_values[r]
        base_rate = p_value === nothing ? interpolate_rate(table, T9_value) : sampled_interpolate_rate(table, T9_value, p_value)
        multiplier = rate_multipliers === nothing ? 1.0 : rate_multipliers[r]
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

# Composition-dependent part of the weak-screening exponent that's shared by
# every reaction in one evaluation: `0.188*sqrt(rho*zeta/T6^3)` from
# `weak_screening_multiplier`, split into its rho/T9-only prefactor
# (`cache.screening_scale`, precomputed once per step) times sqrt(zeta)
# (composition-dependent, recomputed per Y). `_cached_screening_multiplier`
# then only needs `charge_pair_sum * zeta_scale` per reaction.
function _screening_zeta_scale(cache::NetworkStepCache, network::ReactionNetwork, Y::AbstractVector{<:Real})
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

function _screening_context(cache::NetworkStepCache, network::ReactionNetwork, Y::AbstractVector{<:Real})
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
# `dYdt`/`Y` are generic (not hardcoded `Vector{Float64}`) so this same
# function differentiates directly under ForwardDiff when `cache.prefactors`
# (and hence `dYdt`) carries duals -- see `NetworkStepCache`.
function _cached_network_rhs!(dYdt::AbstractVector{<:Real}, network::ReactionNetwork, cache::NetworkStepCache, Y::AbstractVector{<:Real})
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
function _cached_network_jacobian!(J::AbstractMatrix{<:Real}, network::ReactionNetwork, cache::NetworkStepCache, Y::AbstractVector{<:Real})
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

"""
    sparse_jacobian_prototype(network)

A fresh `SparseMatrixCSC{Float64,Int}` sharing `network.jacobian_sparsity`'s
structure (`colptr`/`rowval`), values zeroed. Pass this as the buffer to
`_cached_network_jacobian_sparse!`; its structure never changes across calls,
only `nzval` is overwritten each evaluation (feature spec Tier 0 #2).
"""
function sparse_jacobian_prototype(network::ReactionNetwork)
    pattern = network.jacobian_sparsity
    return SparseMatrixCSC(pattern.m, pattern.n, copy(pattern.colptr), copy(pattern.rowval), zeros(Float64, nnz(pattern)))
end

# Binary search for row `row`'s position within column `col`'s `nzrange` --
# valid because SparseMatrixCSC stores each column's row indices in strictly
# ascending order. Throws (rather than silently no-op'ing) if `(row, col)`
# falls outside the precomputed sparsity pattern: that would mean
# `_jacobian_sparsity_pattern` under-counted a real Jacobian entry, a bug in
# the pattern, not a normal runtime condition.
function _sparse_nzval_index(J::SparseMatrixCSC, row::Int, col::Int)
    range = nzrange(J, col)
    position = searchsortedfirst(view(J.rowval, range), row)
    (position <= length(range) && J.rowval[range[position]] == row) ||
        throw(ArgumentError("(row=$row, col=$col) is not in the Jacobian sparsity pattern"))
    return range[position]
end

# Sparse counterpart of `_cached_network_jacobian!`: identical math (same
# per-reaction product-rule derivative over distinct reactant species), only
# the write target differs -- accumulate into `J.nzval` at the precomputed
# sparsity structure's positions instead of a dense `Matrix` cell. `J` must
# have been built from `network.jacobian_sparsity` (e.g. via
# `sparse_jacobian_prototype`); the structure (`colptr`/`rowval`) is never
# touched here, only values.
function _cached_network_jacobian_sparse!(J::SparseMatrixCSC{Float64,Int}, network::ReactionNetwork, cache::NetworkStepCache, Y::AbstractVector{<:Real})
    fill!(J.nzval, 0.0)
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
                J.nzval[_sparse_nzval_index(J, index, jindex)] -= count * dflux
            end
            for (index, count) in zip(compiled.product_species_indices, compiled.product_species_counts)
                J.nzval[_sparse_nzval_index(J, index, jindex)] += count * dflux
            end
        end
    end

    return J
end


# Evaluate a trajectory quantity at time t, whether it's a constant Real
# (single-zone experiments with fixed rho/T9) or a callable profile (e.g.
# from `trajectory_profiles`). This dispatch is what lets every solver
# accept `rho`/`T9` as either form without duplicating each stepper.
_profile_value(value::Real, t::Real) = Float64(value)
_profile_value(value, t::Real) = Float64(value(t))

# Uncached RHS evaluation at time t: resolve rho(t)/T9(t), then call
# `network_rhs`. Used when no step cache applies (a custom screening
# function, see `_build_step_cache`).
function _rhs_at(network::ReactionNetwork, Y::AbstractVector{<:Real}, rho, T9, t::Real; rate_multipliers=nothing, rate_p_values=nothing, screening=nothing)
    rho_t = _profile_value(rho, t)
    T9_t = _profile_value(T9, t)
    return network_rhs(Y, network, rho_t, T9_t; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)
end

# Resolve rho(t)/T9(t) and build the `NetworkStepCache` for that instant, or
# `nothing` if this screening choice can't be cached (a custom function) --
# the entry point every stepper (`_euler_step`, `_rk4_step`,
# `_backward_euler_step`) uses at the start of a timestep.
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
