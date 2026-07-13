# Weak (Salpeter) and Chugunov (2007) electron screening.

function _screening_composition_factor(network::ReactionNetwork, Y::AbstractVector{<:Real})
    factor = 0.0
    for (i, species) in pairs(network.species_info)
        species.Z <= 0 && continue
        factor += (species.Z^2 + species.Z) * max(Float64(Y[i]), 0.0)
    end
    return factor
end

"""
    weak_screening_multiplier(network, reaction, Y, rho, T9)

Return an approximate weak-screening multiplier for charged-particle reactions.
This is a Salpeter-style diagnostic multiplier using the current abundance
composition. Reactions with fewer than two charged reactants return `1.0`.
"""
function weak_screening_multiplier(
    network::ReactionNetwork,
    reaction::Reaction,
    Y::AbstractVector{<:Real},
    rho::Real,
    T9::Real;
    max_exponent::Real=300.0,
)
    length(reaction.reactants) >= 2 || return 1.0
    T6 = 1000.0 * Float64(T9)
    T6 > 0.0 || throw(ArgumentError("T9 must be positive for screening"))
    rho_value = Float64(rho)
    rho_value > 0.0 || throw(ArgumentError("rho must be positive for screening"))

    zeta = _screening_composition_factor(network, Y)
    zeta > 0.0 || return 1.0

    exponent = 0.0
    reactant_info = [species_from_name(name) for name in reaction.reactants]
    for i in 1:(length(reactant_info)-1)
        Zi = reactant_info[i].Z
        Zi <= 0 && continue
        for j in (i+1):length(reactant_info)
            Zj = reactant_info[j].Z
            Zj <= 0 && continue
            exponent += 0.188 * Zi * Zj * sqrt(rho_value * zeta / T6^3)
        end
    end

    exponent <= 0.0 && return 1.0
    return exp(min(exponent, Float64(max_exponent)))
end

# CGS constants for the Chugunov screening evaluation (CODATA 2018).
const _ELEMENTARY_CHARGE_ESU = 4.80320471257e-10
const _BOLTZMANN_ERG_PER_K = 1.380649e-16
const _HBAR_ERG_S = 1.054571817e-27
const _ATOMIC_MASS_UNIT_G = 1.66053906660e-24

# Half-cosine transition between y=x and y=limit, starting at x=start.
# Ported from pynucastro's screening module.
function _smooth_clip(x::Float64, limit::Float64, start::Float64)
    lower, upper = limit < start ? (limit, x) : (x, limit)
    x < min(limit, start) && return lower
    x > max(limit, start) && return upper
    fraction = (1.0 - cos(pi * (x - min(limit, start)) / (start - limit))) / 2.0
    return (1.0 - fraction) * lower + fraction * upper
end

#=
Composition-dependent plasma quantities for Chugunov screening: the electron
number density and the temperature-independent part of the electron Coulomb
coupling parameter. `n_e = rho * sum(Z_i Y_i) / m_u`.
=#
function _ion_plasma_state(network::ReactionNetwork, Y::AbstractVector{<:Real}, rho::Real)
    charge_abundance = 0.0
    for (i, species) in pairs(network.species_info)
        species.Z <= 0 && continue
        charge_abundance += species.Z * max(Float64(Y[i]), 0.0)
    end
    charge_abundance > 0.0 || return nothing

    n_e = Float64(rho) * charge_abundance / _ATOMIC_MASS_UNIT_G
    gamma_e_fac = _ELEMENTARY_CHARGE_ESU^2 / _BOLTZMANN_ERG_PER_K * cbrt(4.0 * pi / 3.0) * cbrt(n_e)
    return (n_e=n_e, gamma_e_fac=gamma_e_fac)
end

#=
Screening exponent h = ln(f_screen) of one ion pair following Chugunov,
DeWitt & Yakovlev (2007), extended to multi-component plasmas as in Yakovlev
et al. (2006). Ported from pynucastro's `chugunov_2007`, which documents the
substitutions (Z^2 -> Z1*Z2, n_i -> n_e/ztilde^3, m_i -> 2*mu12*m_u). Valid
from the weak-screening limit through the strong-screening regime; the fit
caps at Gamma ~ 600 and T ~ 0.1 T_p.
=#
function _chugunov_pair_exponent(Z1::Float64, A1::Float64, Z2::Float64, A2::Float64, temperature_K::Float64, n_e::Float64, gamma_e_fac::Float64)
    ztilde = (cbrt(Z1) + cbrt(Z2)) / 2.0
    reduced_mass = A1 * A2 / (A1 + A2)
    n_i = n_e / ztilde^3
    m_i = 2.0 * reduced_mass * _ATOMIC_MASS_UNIT_G

    T_p = _HBAR_ERG_S / _BOLTZMANN_ERG_PER_K * _ELEMENTARY_CHARGE_ESU * sqrt(4.0 * pi * Z1 * Z2 * n_i / m_i)
    T_norm = _smooth_clip(temperature_K / T_p, 0.1, 0.2)
    Gamma = _smooth_clip(gamma_e_fac * Z1 * Z2 / (ztilde * T_norm * T_p), 600.0, 590.0)

    zeta = cbrt(4.0 / (3.0 * pi^2 * T_norm^2))
    poly = 1.0 + zeta * (0.022 + zeta * ((0.41 - 0.6 / Gamma) + (0.06 + 2.2 / Gamma) * zeta))
    gamtilde = Gamma / cbrt(poly)
    gamtilde2 = gamtilde^2

    A1_fit = 2.7822
    A2_fit = 98.34
    A3_fit = sqrt(3.0) - A1_fit / sqrt(A2_fit)
    B1_fit = -1.7476
    B2_fit = 66.07
    B3_fit = 1.12
    B4_fit = 65.0

    h = gamtilde^1.5 * (A1_fit / sqrt(A2_fit + gamtilde) + A3_fit / (1.0 + gamtilde)) +
        B1_fit * gamtilde2 / (B2_fit + gamtilde) +
        B3_fit * gamtilde2 / (B4_fit + gamtilde2)
    return max(h, 0.0)
end

function _chugunov_reaction_multiplier(compiled::CompiledReaction, context, T9::Real)
    (!context.plasma_active || isempty(compiled.screening_pairs)) && return 1.0
    temperature_K = 1.0e9 * Float64(T9)
    exponent = 0.0
    for (Z1, A1, Z2, A2) in compiled.screening_pairs
        exponent += _chugunov_pair_exponent(Z1, A1, Z2, A2, temperature_K, context.n_e, context.gamma_e_fac)
    end
    return exp(min(exponent, _SCREENING_MAX_EXPONENT))
end

"""
    chugunov_screening_multiplier(network, reaction, Y, rho, T9)

Chugunov, DeWitt & Yakovlev (2007) screening multiplier for one reaction at
the current composition. Available network-wide with `screening=:chugunov`.
"""
function chugunov_screening_multiplier(
    network::ReactionNetwork,
    reaction::Reaction,
    Y::AbstractVector{<:Real},
    rho::Real,
    T9::Real,
)
    T9 > 0.0 || throw(ArgumentError("T9 must be positive for screening"))
    rho > 0.0 || throw(ArgumentError("rho must be positive for screening"))
    length(reaction.reactants) >= 2 || return 1.0

    plasma = _ion_plasma_state(network, Y, rho)
    plasma === nothing && return 1.0

    temperature_K = 1.0e9 * Float64(T9)
    exponent = 0.0
    accumulated_Z = 0
    accumulated_A = 0
    for (k, name) in pairs(reaction.reactants)
        species = try
            species_from_name(name)
        catch
            Species(String(name), 0, 0)
        end
        if k > 1 && accumulated_Z > 0 && species.Z > 0
            exponent += _chugunov_pair_exponent(
                Float64(accumulated_Z), Float64(accumulated_A),
                Float64(species.Z), Float64(species.A),
                temperature_K, plasma.n_e, plasma.gamma_e_fac,
            )
        end
        accumulated_Z += species.Z
        accumulated_A += species.A
    end

    return exp(min(exponent, _SCREENING_MAX_EXPONENT))
end

function _screening_multiplier(screening, network::ReactionNetwork, reaction::Reaction, Y::AbstractVector{<:Real}, rho::Real, T9::Real)
    if screening === nothing || screening === false
        return 1.0
    elseif screening == :weak
        return weak_screening_multiplier(network, reaction, Y, rho, T9)
    elseif screening == :chugunov
        return chugunov_screening_multiplier(network, reaction, Y, rho, T9)
    elseif screening isa Function
        return Float64(screening(network, reaction, Y, rho, T9))
    end

    throw(ArgumentError("unsupported screening=$screening; use nothing, :weak, :chugunov, or a function `(network, reaction, Y, rho, T9) -> multiplier`"))
end

const _SCREENING_MAX_EXPONENT = 300.0

