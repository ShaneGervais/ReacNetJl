# Post-processing weak decay of a frozen composition.

function _decay_generator_matrix(network::ReactionNetwork, T9::Real)
    n = length(network.species)
    generator = zeros(Float64, n, n)
    for (reaction, compiled) in zip(network.reactions, network.compiled_reactions)
        length(compiled.reactant_species_indices) == 1 || continue
        source_index = only(compiled.reactant_species_indices)
        rate = interpolate_rate(reaction.rate_table, T9)
        for row in 1:n
            delta = compiled.stoichiometric_delta[row]
            delta == 0.0 && continue
            generator[row, source_index] += delta * rate
        end
    end
    return generator
end

function decay_mass_fractions(
    tables::AbstractVector{ReactionRateTable},
    X::AbstractDict,
    decay_time_s::Real;
    T9::Real=0.1,
)
    decay_time = Float64(decay_time_s)
    decay_time >= 0.0 || throw(ArgumentError("decay_time_s must be non-negative"))

    normalized_X = Dict(normalize_species_name(string(name)) => Float64(value) for (name, value) in X)
    decay_time == 0.0 && return (
        mass_fractions=copy(normalized_X),
        decay_tables=ReactionRateTable[],
        network=nothing,
        times=Float64[0.0],
    )

    decay_tables = select_decay_reaction_tables(tables, keys(normalized_X))
    isempty(decay_tables) && return (
        mass_fractions=copy(normalized_X),
        decay_tables=decay_tables,
        network=nothing,
        times=Float64[0.0, decay_time],
    )

    species = sort!(collect(union(Set(keys(normalized_X)), Set(_infer_species_from_reactions(Reaction.(decay_tables))))); by=name -> begin
        info = species_from_name(name)
        (info.Z, info.A, name)
    end)
    network = network_from_tables(decay_tables; species=species)
    X0 = Dict(name => get(normalized_X, name, 0.0) for name in network.species)
    Y0 = abundances_from_mass_fractions(network, X0)

    generator = _decay_generator_matrix(network, T9)
    Y_final = exp(decay_time * generator) * Y0
    for i in eachindex(Y_final)
        if Y_final[i] < 0.0 && abs(Y_final[i]) <= 1.0e-30
            Y_final[i] = 0.0
        end
    end
    decayed_X = mass_fractions_from_abundances(network, Y_final)
    for (name, value) in decayed_X
        normalized_X[name] = value
    end

    return (
        mass_fractions=normalized_X,
        decay_tables=decay_tables,
        network=network,
        times=Float64[0.0, decay_time],
    )
end

