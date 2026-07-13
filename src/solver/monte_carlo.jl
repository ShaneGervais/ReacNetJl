# STARLIB uncertainty Monte Carlo driver.

"""
    run_monte_carlo(network, Y0, tspan, dt, rho, T9; nruns, seed=nothing, method=:rk4, store_histories=false)

Run repeated single-zone network calculations with STARLIB lognormal rate
sampling. For each run and each reaction, a random `p ~ Normal(0, 1)` is drawn
and held fixed for that reaction throughout the run:

    sampled_rate(T) = recommended_rate(T) * factor_uncertainty(T)^p

Returns a named tuple containing final abundances, sampled `p` values, and
optionally all abundance histories.
"""
function run_monte_carlo(
    network::ReactionNetwork,
    Y0::AbstractVector{<:Real},
    tspan::Tuple{<:Real,<:Real},
    dt::Real,
    rho,
    T9;
    nruns::Integer,
    seed=nothing,
    method::Symbol=:rk4,
    rate_multipliers=nothing,
    clamp_negative::Bool=true,
    store_histories::Bool=false,
    screening=nothing,
)
    nruns > 0 || throw(ArgumentError("nruns must be positive"))

    rng = seed === nothing ? Random.default_rng() : Random.MersenneTwister(seed)
    nreactions = length(network.reactions)
    nspecies = length(network.species)
    sampled_p_values = randn(rng, nruns, nreactions)
    final_abundances = Matrix{Float64}(undef, nruns, nspecies)
    histories = Matrix{Float64}[]
    saved_times = Float64[]

    for run_index in 1:nruns
        p_values = view(sampled_p_values, run_index, :)
        times, history = solve_network(
            network,
            Y0,
            tspan,
            dt,
            rho,
            T9;
            method=method,
            rate_multipliers=rate_multipliers,
            rate_p_values=p_values,
            clamp_negative=clamp_negative,
            screening=screening,
        )

        if run_index == 1
            saved_times = collect(times)
        end

        final_abundances[run_index, :] .= history[end, :]
        store_histories && push!(histories, history)
    end

    return (
        times=saved_times,
        final_abundances=final_abundances,
        rate_p_values=sampled_p_values,
        histories=histories,
        species=copy(network.species),
        reactions=[reaction_string(reaction) for reaction in network.reactions],
    )
end

