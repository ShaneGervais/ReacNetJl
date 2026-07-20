# Reaction fluxes, the network right-hand side, and energy generation.

"""
    reaction_flux(reaction, Y, species_index, rho, T9; rate_multiplier=1.0)

Calculate the abundance flux `F_r` for one reaction `r`, in mol g^-1 s^-1.

`Y` is the abundance vector (mol g^-1). `species_index` maps species names to
positions in `Y`. `rho` is mass density in g cm^-3, and `T9` is temperature in
GK. `rate_multiplier` scales the interpolated rate (used for STARLIB
factor-uncertainty sampling and Monte Carlo studies); `rate_p_value`, if given,
samples the rate at a fixed lognormal quantile via `sampled_interpolate_rate`
instead of using the tabulated central value.

For a reaction with `N_r` reactants (counting repeats) and rate `R_r(T_9)`:

```math
F_r = \\frac{\\rho^{N_r-1} R_r(T_9) \\prod_j Y_j^{\\nu^{\\mathrm{react}}_{j,r}}}{\\prod_j \\nu^{\\mathrm{react}}_{j,r}!}
```

where `\\nu^{\\mathrm{react}}_{j,r}` is the number of times species `j` appears
as a reactant. This reduces to `F_r = R_r Y_i` for a one-body decay/disintegration
(`N_r = 1`, no density factor), `F_r = \\rho R_r Y_i Y_j` for a two-body capture
of distinct species, and includes the `1/2!` symmetry factor from
`_symmetry_factor` for reactions with two identical reactants (e.g. `p+p`,
`α+α`) so the same physical collision isn't double-counted.
"""
function reaction_flux(
    reaction::Reaction,
    Y::AbstractVector{<:Real},
    species_index::AbstractDict,
    rho::Real,
    T9::Real;
    rate_multiplier::Real=1.0,
    rate_p_value=nothing,
)
    nreactants = length(reaction.reactants)
    nreactants >= 1 || throw(ArgumentError("reaction must have at least one reactant"))

    base_rate = rate_p_value === nothing ? interpolate_rate(reaction.rate_table, T9) : sampled_interpolate_rate(reaction.rate_table, T9, rate_p_value)
    rate = rate_multiplier * base_rate
    abundance_product = 1.0
    for name in reaction.reactants
        abundance_product *= Y[_species_index(species_index, name)]
    end

    density_factor = Float64(rho)^(nreactants - 1)
    return density_factor * rate * abundance_product / _symmetry_factor(reaction.reactants)
end

function _reaction_flux(
    network::ReactionNetwork,
    reaction::Reaction,
    compiled::CompiledReaction,
    Y::AbstractVector{<:Real},
    rho::Real,
    T9::Real;
    rate_multiplier::Real=1.0,
    rate_p_value=nothing,
    screening=nothing,
)
    compiled.nreactants >= 1 || throw(ArgumentError("reaction must have at least one reactant"))

    base_rate = rate_p_value === nothing ? interpolate_rate(reaction.rate_table, T9) : sampled_interpolate_rate(reaction.rate_table, T9, rate_p_value)
    rate = rate_multiplier * _screening_multiplier(screening, network, reaction, Y, rho, T9) * base_rate
    abundance_product = 1.0
    for index in compiled.reactant_indices
        abundance_product *= Y[index]
    end

    density_factor = Float64(rho)^(compiled.nreactants - 1)
    return density_factor * rate * abundance_product / compiled.symmetry_factor
end

"""
    network_rhs(Y, reactions, species_index, rho, T9; rate_multipliers=nothing)

Calculate `dY/dt` for a single-zone reaction network at fixed density and
temperature (the label-driven form; see the `ReactionNetwork` method below for
the precompiled, performance-path form used by the solvers).

For each reaction `r` with flux `F_r` (see `reaction_flux`), every reactant
loses `F_r` and every product gains `F_r`, so for species `i`:

```math
\\frac{dY_i}{dt} = \\sum_r \\left(\\nu^{\\mathrm{prod}}_{i,r} - \\nu^{\\mathrm{react}}_{i,r}\\right) F_r
```

where `\\nu^{\\mathrm{prod}}_{i,r}` and `\\nu^{\\mathrm{react}}_{i,r}` count how
many times species `i` appears as a product or reactant of reaction `r`
(usually 0 or 1, but 2 for e.g. `p+p -> d + 2p`-style multiproduct channels).

`rate_multipliers` can be supplied as a vector the same length as `reactions`,
one scalar multiplier per reaction. This is the hook used for STARLIB
factor-uncertainty sampling and Monte Carlo uncertainty studies (`rate_p_values`
does the same for lognormal-quantile sampling).
"""
function network_rhs(
    Y::AbstractVector{<:Real},
    reactions::AbstractVector{Reaction},
    species_index::AbstractDict,
    rho::Real,
    T9::Real;
    rate_multipliers=nothing,
    rate_p_values=nothing,
)
    dYdt = zeros(Float64, length(Y))

    if rate_multipliers !== nothing && length(rate_multipliers) != length(reactions)
        throw(ArgumentError("rate_multipliers must have the same length as reactions"))
    end
    if rate_p_values !== nothing && length(rate_p_values) != length(reactions)
        throw(ArgumentError("rate_p_values must have the same length as reactions"))
    end

    for (reaction_number, reaction) in pairs(reactions)
        multiplier = rate_multipliers === nothing ? 1.0 : rate_multipliers[reaction_number]
        p_value = rate_p_values === nothing ? nothing : rate_p_values[reaction_number]
        flux = reaction_flux(reaction, Y, species_index, rho, T9; rate_multiplier=multiplier, rate_p_value=p_value)

        for name in reaction.reactants
            dYdt[_species_index(species_index, name)] -= flux
        end

        for name in reaction.products
            dYdt[_species_index(species_index, name)] += flux
        end
    end

    return dYdt
end

function network_rhs(
    Y::AbstractVector{<:Real},
    network::ReactionNetwork,
    rho::Real,
    T9::Real;
    rate_multipliers=nothing,
    rate_p_values=nothing,
    screening=nothing,
)
    dYdt = zeros(Float64, length(Y))

    length(Y) == length(network.species) || throw(ArgumentError("Y length must match the number of network species"))
    if rate_multipliers !== nothing && length(rate_multipliers) != length(network.reactions)
        throw(ArgumentError("rate_multipliers must have the same length as reactions"))
    end
    if rate_p_values !== nothing && length(rate_p_values) != length(network.reactions)
        throw(ArgumentError("rate_p_values must have the same length as reactions"))
    end

    for (reaction_number, reaction) in pairs(network.reactions)
        compiled = network.compiled_reactions[reaction_number]
        multiplier = rate_multipliers === nothing ? 1.0 : rate_multipliers[reaction_number]
        p_value = rate_p_values === nothing ? nothing : rate_p_values[reaction_number]
        flux = _reaction_flux(network, reaction, compiled, Y, rho, T9; rate_multiplier=multiplier, rate_p_value=p_value, screening=screening)

        for (index, count) in zip(compiled.reactant_species_indices, compiled.reactant_species_counts)
            dYdt[index] -= count * flux
        end
        for (index, count) in zip(compiled.product_species_indices, compiled.product_species_counts)
            dYdt[index] += count * flux
        end
    end

    return dYdt
end

"""
    reaction_fluxes(network, Y, rho, T9; rate_multipliers=nothing)

Calculate the instantaneous flux `F_r` (see `reaction_flux`) of every reaction
in a precompiled `network`, using the network's `CompiledReaction` bookkeeping
(reactant/product indices, symmetry factors) for speed. The returned vector is
ordered like `network.reactions`.
"""
function reaction_fluxes(
    network::ReactionNetwork,
    Y::AbstractVector{<:Real},
    rho::Real,
    T9::Real;
    rate_multipliers=nothing,
    rate_p_values=nothing,
    screening=nothing,
)
    if rate_multipliers !== nothing && length(rate_multipliers) != length(network.reactions)
        throw(ArgumentError("rate_multipliers must have the same length as network.reactions"))
    end
    if rate_p_values !== nothing && length(rate_p_values) != length(network.reactions)
        throw(ArgumentError("rate_p_values must have the same length as network.reactions"))
    end

    fluxes = zeros(Float64, length(network.reactions))
    for (i, reaction) in pairs(network.reactions)
        multiplier = rate_multipliers === nothing ? 1.0 : rate_multipliers[i]
        p_value = rate_p_values === nothing ? nothing : rate_p_values[i]
        fluxes[i] = _reaction_flux(network, reaction, network.compiled_reactions[i], Y, rho, T9; rate_multiplier=multiplier, rate_p_value=p_value, screening=screening)
    end
    return fluxes
end

"""
    reaction_flux_history(network, history, times, rho, T9; rate_multipliers=nothing)

Calculate `F_r(t_n)` (see `reaction_flux`) at every saved timestep from
`solve_network`/`solve_network_adaptive` output. `rho` and `T9` may be
constants or functions of time (trajectory profiles), evaluated once per saved
time via `_profile_value`.

Returns a matrix where `flux_history[n, r]` is the flux of reaction `r` at
`times[n]`. This is the per-timestep flux CSV data (`reaction_fluxes.csv` via
`write_reaction_flux_csv`) and the input to `integrated_fluxes`.
"""
function reaction_flux_history(
    network::ReactionNetwork,
    history::AbstractMatrix{<:Real},
    times::AbstractVector{<:Real},
    rho,
    T9;
    rate_multipliers=nothing,
    rate_p_values=nothing,
    screening=nothing,
)
    length(times) == size(history, 1) || throw(ArgumentError("times length must match the number of history rows"))
    size(history, 2) == length(network.species) || throw(ArgumentError("history column count must match the number of network species"))

    flux_history = Matrix{Float64}(undef, length(times), length(network.reactions))
    for (n, t) in pairs(times)
        rho_t = _profile_value(rho, t)
        T9_t = _profile_value(T9, t)
        flux_history[n, :] .= reaction_fluxes(network, view(history, n, :), rho_t, T9_t; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)
    end
    return flux_history
end

"""
    integrated_fluxes(times, flux_history)

Integrate reaction flux histories in time using the trapezoid rule, giving the
total abundance processed through each reaction over the run:

```math
\\Phi_r = \\int F_r(t)\\, dt \\approx \\sum_n \\tfrac{1}{2}(t_{n+1}-t_n)\\left(F_r(t_n) + F_r(t_{n+1})\\right)
```

`\\Phi_r` (mol g^-1) is the basis for the "largest abundance residual" flux
reports: multiplying by mass number `A` and the reaction's net stoichiometric
change for a species gives that reaction's total contribution to the species'
mass-fraction change over the run. Returns one integrated flux per reaction,
ordered like `network.reactions`.
"""
function integrated_fluxes(times::AbstractVector{<:Real}, flux_history::AbstractMatrix{<:Real})
    length(times) == size(flux_history, 1) || throw(ArgumentError("times length must match the number of flux-history rows"))
    length(times) >= 2 || throw(ArgumentError("at least two time points are required"))

    totals = zeros(Float64, size(flux_history, 2))
    for n in 1:(length(times)-1)
        dt = Float64(times[n+1] - times[n])
        dt >= 0.0 || throw(ArgumentError("times must be monotonically increasing"))
        totals .+= 0.5 * dt .* (view(flux_history, n, :) .+ view(flux_history, n + 1, :))
    end
    return totals
end

"""
    energy_generation_rate(network, Y, rho, T9; ...)

Return diagnostic nuclear energy generation `\\epsilon_{\\mathrm{nuc}}` in
erg g^-1 s^-1 from reaction Q-values (MeV, from each rate table's `q_value`)
and instantaneous reaction fluxes:

```math
\\epsilon_{\\mathrm{nuc}} = N_A \\cdot \\mathrm{MeV\\_to\\_erg} \\sum_r Q_r F_r
```

This is diagnostic only: it reports how much energy the current abundance
flow *would* release, but it does not feed back into `T9`/`rho` — single-zone
PPN evolves the network along a prescribed trajectory, it does not solve
`dT/dt`.
"""
function energy_generation_rate(
    network::ReactionNetwork,
    Y::AbstractVector{<:Real},
    rho::Real,
    T9::Real;
    rate_multipliers=nothing,
    rate_p_values=nothing,
    screening=nothing,
)
    fluxes = reaction_fluxes(network, Y, rho, T9; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)
    epsilon = 0.0
    for (i, reaction) in pairs(network.reactions)
        epsilon += reaction.rate_table.q_value * fluxes[i]
    end
    return epsilon * AVOGADRO * MEV_TO_ERG
end

"""
    energy_generation_history(network, history, times, rho, T9; ...)

Return diagnostic nuclear energy generation `\\epsilon_{\\mathrm{nuc}}(t_n)`
(see `energy_generation_rate`) in erg g^-1 s^-1 at every saved history row.
"""
function energy_generation_history(
    network::ReactionNetwork,
    history::AbstractMatrix{<:Real},
    times::AbstractVector{<:Real},
    rho,
    T9;
    rate_multipliers=nothing,
    rate_p_values=nothing,
    screening=nothing,
)
    length(times) == size(history, 1) || throw(ArgumentError("times length must match the number of history rows"))
    size(history, 2) == length(network.species) || throw(ArgumentError("history column count must match the number of network species"))

    epsilon = Vector{Float64}(undef, length(times))
    for (n, t) in pairs(times)
        rho_t = _profile_value(rho, t)
        T9_t = _profile_value(T9, t)
        epsilon[n] = energy_generation_rate(network, view(history, n, :), rho_t, T9_t; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)
    end
    return epsilon
end

"""
    integrated_energy_generation(times, epsilon_history)

Integrate diagnostic energy generation in time using the trapezoid rule:

```math
E_{\\mathrm{nuc}} = \\int \\epsilon_{\\mathrm{nuc}}(t)\\, dt
```

Returns specific energy release in erg g^-1 over the run (diagnostic only,
see `energy_generation_rate`).
"""
function integrated_energy_generation(times::AbstractVector{<:Real}, epsilon_history::AbstractVector{<:Real})
    length(times) == length(epsilon_history) || throw(ArgumentError("times length must match energy-history length"))
    length(times) >= 2 || throw(ArgumentError("at least two time points are required"))

    total = 0.0
    for n in 1:(length(times)-1)
        dt = Float64(times[n+1] - times[n])
        dt >= 0.0 || throw(ArgumentError("times must be monotonically increasing"))
        total += 0.5 * dt * (epsilon_history[n] + epsilon_history[n+1])
    end
    return total
end

"""
    species_flux_balance(network, Y, rho, T9; rate_multipliers=nothing)

Calculate instantaneous production, destruction, and net `dY/dt` contributions
for every species, split out by direction rather than summed as in
`network_rhs`:

```math
P_i = \\sum_r \\nu^{\\mathrm{prod}}_{i,r} F_r \\qquad
D_i = \\sum_r \\nu^{\\mathrm{react}}_{i,r} F_r \\qquad
\\left(\\frac{dY_i}{dt}\\right)_{\\mathrm{net}} = P_i - D_i
```

Splitting production from destruction (rather than only the net `dY/dt`) is
what lets `species_flux_contributions`-style diagnostics identify which
individual reactions dominate a species' formation versus depletion, even
when the two nearly cancel.

Returns `(production=..., destruction=..., net=...)`, with vectors ordered like
`network.species`.
"""
function species_flux_balance(
    network::ReactionNetwork,
    Y::AbstractVector{<:Real},
    rho::Real,
    T9::Real;
    rate_multipliers=nothing,
    rate_p_values=nothing,
    screening=nothing,
)
    fluxes = reaction_fluxes(network, Y, rho, T9; rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening)
    production = zeros(Float64, length(network.species))
    destruction = zeros(Float64, length(network.species))

    for (reaction_index, reaction) in pairs(network.reactions)
        compiled = network.compiled_reactions[reaction_index]
        flux = fluxes[reaction_index]

        for (index, count) in zip(compiled.product_species_indices, compiled.product_species_counts)
            production[index] += count * flux
        end

        for (index, count) in zip(compiled.reactant_species_indices, compiled.reactant_species_counts)
            destruction[index] += count * flux
        end
    end

    return (production=production, destruction=destruction, net=production .- destruction)
end

"""
    reaction_edges(network)

Return a flattened graph-like edge list for diagnostics or external plotting
(e.g. drawing a reaction-network chart). A reaction with multiple reactants
and products is really a hyperedge (several inputs to several outputs at
once); this flattens each reaction into the full reactant x product cross
product of ordinary from -> to edges, so downstream graph tools that only
understand simple directed edges can still render it. `reaction_index` lets
callers map an edge back to its originating reaction (and its flux).
"""
function reaction_edges(network::ReactionNetwork)
    edges = NamedTuple[]
    for (reaction_index, reaction) in pairs(network.reactions)
        label = reaction_string(reaction)
        for reactant in reaction.reactants
            for product in reaction.products
                push!(edges, (reaction_index=reaction_index, reaction=label, from=reactant, to=product))
            end
        end
    end
    return edges
end

