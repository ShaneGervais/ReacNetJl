# Reaction fluxes, the network right-hand side, and energy generation.

"""
    reaction_flux(reaction, Y, species_index, rho, T9; rate_multiplier=1.0)

Calculate the abundance flux for one reaction.

`Y` is the abundance vector. `species_index` maps species names to positions in
`Y`. `rho` is mass density in g cm^-3, and `T9` is temperature in GK.

For a one-body reaction, the flux is `rate * Yᵢ`. For a two-body reaction using
STARLIB's usual `N_A <σv>` rate, the flux is `rho * rate * Yᵢ * Yⱼ`, with the
standard symmetry correction for identical reactants.
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

#=
    network_rhs(Y, reactions, species_index, rho, T9; rate_multipliers=nothing)

Calculate `dY/dt` for a single-zone reaction network at fixed density and
temperature.

`rate_multipliers` can be supplied as a vector the same length as `reactions`.
This gives us a simple hook for later STARLIB uncertainty studies.
=#
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

#=
    reaction_fluxes(network, Y, rho, T9; rate_multipliers=nothing)

Calculate instantaneous reaction fluxes for every reaction in a network.
The returned vector is ordered like `network.reactions`.
=#
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

#=
    reaction_flux_history(network, history, times, rho, T9; rate_multipliers=nothing)

Calculate reaction fluxes at every saved timestep from `solve_network` output.
`rho` and `T9` may be constants or functions of time.

Returns a matrix where `flux_history[n, r]` is the flux of reaction `r` at
`times[n]`.
=#
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

#=
    integrated_fluxes(times, flux_history)

Integrate reaction flux histories in time using the trapezoid rule.
Returns one integrated flux per reaction.
=#
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

#=
    energy_generation_rate(network, Y, rho, T9; ...)

Return diagnostic nuclear energy generation in erg g^-1 s^-1 from reaction
Q-values and instantaneous reaction fluxes. This does not feed back into the
temperature trajectory.
=#
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

#=
    energy_generation_history(network, history, times, rho, T9; ...)

Return diagnostic nuclear energy generation in erg g^-1 s^-1 at every saved
history row.
=#
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

#=
    integrated_energy_generation(times, epsilon_history)

Integrate diagnostic energy generation in time using the trapezoid rule.
Returns specific energy release in erg g^-1.
=#
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
for every species.

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

#=
    reaction_edges(network)

Return a flattened graph-like edge list for diagnostics or external plotting.
Each edge connects one reactant species to one product species for a reaction.
This is a graph approximation of the reaction hypergraph.
=#
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

