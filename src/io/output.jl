# CSV writers for mass fractions, reaction fluxes, and the network definition.


# Quote a CSV field per RFC 4180 if it contains a comma (reaction labels can,
# e.g. `"p(n,γ)d"`), doubling any embedded quotes; fields without a comma are
# left bare for readability.
_csv_quote(s::AbstractString) = occursin(',', s) ? "\"" * replace(s, "\"" => "\"\"") * "\"" : s

"""
    write_mass_fraction_csv(path, network, times, history; profiles=nothing)

Write mass fractions over time to a CSV: one row per saved time, one column
per network species. If `profiles` (as returned by `trajectory_profiles`) is
given, `T9` and `rho` columns are included. Returns `path`.
"""
function write_mass_fraction_csv(
    path::AbstractString,
    network::ReactionNetwork,
    times::AbstractVector{<:Real},
    history::AbstractMatrix{<:Real};
    profiles=nothing,
)
    mkpath(dirname(path))
    open(path, "w") do io
        header = vcat(["time_s"], profiles === nothing ? String[] : ["T9", "rho"], network.species)
        println(io, join(header, ","))
        for (n, t) in pairs(times)
            row = [@sprintf("%.16E", t)]
            if profiles !== nothing
                push!(row, @sprintf("%.16E", profiles.T9(t)))
                push!(row, @sprintf("%.16E", profiles.rho(t)))
            end
            X = mass_fractions_from_abundances(network, view(history, n, :))
            append!(row, [@sprintf("%.16E", X[name]) for name in network.species])
            println(io, join(row, ","))
        end
    end
    return path
end

"""
    write_reaction_flux_csv(path, network, times, flux_history)

Write reaction fluxes over time to a CSV: one row per saved time, one column
per network reaction (labeled with `reaction_string`). Returns `path`.
"""
function write_reaction_flux_csv(
    path::AbstractString,
    network::ReactionNetwork,
    times::AbstractVector{<:Real},
    flux_history::AbstractMatrix{<:Real},
)
    mkpath(dirname(path))
    labels = [_csv_quote(reaction_string(reaction)) for reaction in network.reactions]
    open(path, "w") do io
        println(io, join(vcat(["time_s"], labels), ","))
        for (n, t) in pairs(times)
            row = [@sprintf("%.16E", t)]
            append!(row, [@sprintf("%.16E", flux_history[n, r]) for r in eachindex(network.reactions)])
            println(io, join(row, ","))
        end
    end
    return path
end

"""
    write_network_csv(path, network)

Write the reaction network definition to a CSV: one row per reaction, with
its label, reactants, products, chapter, rate source, and Q-value (MeV).
Returns `path`.
"""
function write_network_csv(path::AbstractString, network::ReactionNetwork)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "index,reaction,reactants,products,chapter,source,q_value_mev")
        for (i, reaction) in pairs(network.reactions)
            row = [
                string(i),
                _csv_quote(reaction_string(reaction)),
                join(reaction.reactants, "+"),
                join(reaction.products, "+"),
                string(reaction.rate_table.chapter),
                reaction.rate_table.source,
                @sprintf("%.6f", reaction.rate_table.q_value),
            ]
            println(io, join(row, ","))
        end
    end
    return path
end

