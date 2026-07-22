#=
Run ReacNetJl's single-zone post-processing network on a trajectory and
initial-abundance file, writing the standard CSV outputs (mass_fractions.csv,
reaction_fluxes.csv, integrated_fluxes.csv, network.csv) plus a
final_state.csv snapshot that decay_ppn.jl can later decay without
re-running the solve.

Usage:
    julia --project=. run_ppn.jl <trajectory_path> <abundance_path> [output_dir] \
        [--rates starlib|iliadis2002] [--screening weak|chugunov|none] \
        [--factor "label=value"]... [--mc-factor "label"]... [--seed N]

Defaults: output_dir="outputs/run_ppn", rates=iliadis2002, screening=chugunov

Examples:
    # Double one named reaction's rate for the whole run (Iliadis-2002-style
    # deterministic sensitivity factoring):
    julia --project=. run_ppn.jl trajectory.input initial_abundance.dat \
        --factor "22Na(p,g)23Mg=2.0"

    # Sample one named reaction's rate from its STARLIB/NACRE lognormal
    # factor uncertainty instead (fixed p ~ Normal(0,1) for this run; repeat
    # the run with different --seed values to build up a distribution):
    julia --project=. run_ppn.jl trajectory.input initial_abundance.dat \
        --mc-factor "16O(p,g)17F" --seed 7
=#

using ReacNetJl
using Random

function parse_args(args)
    positional = String[]
    rates = :iliadis2002
    screening = :chugunov
    factors = Dict{String,Float64}()
    mc_labels = String[]
    seed = nothing
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--rates"
            rates = Symbol(args[i+1])
            i += 2
        elseif arg == "--screening"
            s = args[i+1]
            screening = s == "none" ? nothing : Symbol(s)
            i += 2
        elseif arg == "--factor"
            spec = args[i+1]
            eq = findlast('=', spec)
            eq === nothing && error("--factor must look like \"label=value\", got: $spec")
            factors[spec[1:prevind(spec, eq)]] = parse(Float64, spec[nextind(spec, eq):end])
            i += 2
        elseif arg == "--mc-factor"
            push!(mc_labels, args[i+1])
            i += 2
        elseif arg == "--seed"
            seed = parse(Int, args[i+1])
            i += 2
        else
            push!(positional, arg)
            i += 1
        end
    end
    return positional, rates, screening, factors, mc_labels, seed
end

const USAGE = """
Usage: julia --project=. run_ppn.jl <trajectory_path> <abundance_path> [output_dir] \\
    [--rates starlib|iliadis2002] [--screening weak|chugunov|none] \\
    [--factor "label=value"]... [--mc-factor "label"]... [--seed N]
"""

positional, rates, screening, factors, mc_labels, seed = parse_args(ARGS)
length(positional) >= 2 || error(USAGE)

trajectory_path = positional[1]
abundance_path = positional[2]
output_dir = length(positional) >= 3 ? positional[3] : "outputs/run_ppn"
rng = seed === nothing ? Random.default_rng() : MersenneTwister(seed)

println("trajectory = ", trajectory_path)
println("abundances = ", abundance_path)
println("rates      = ", rates)
println("screening  = ", screening === nothing ? "none" : screening)
println("output_dir = ", output_dir)
isempty(factors) || println("factors    = ", factors)
isempty(mc_labels) || println("mc-factor  = ", mc_labels, "  (seed=", seed === nothing ? "none (non-reproducible)" : seed, ")")
println()

result = run_ppn(
    trajectory_path, abundance_path;
    rates=rates, screening=screening, output_dir=output_dir,
    rate_factors=isempty(factors) ? nothing : factors,
    rate_sample_labels=isempty(mc_labels) ? nothing : mc_labels,
    rng=rng,
)

println("validated reactions      = ", result.validation.num_reactions)
println("network species          = ", length(result.network.species))
println("solver accepted/rejected = ", result.solver_stats.accepted_steps, " / ", result.solver_stats.rejected_steps)
println("mass drift               = ", result.mass_fraction_drift.initial, " -> ", result.mass_fraction_drift.final)
println("output files              = ", result.output_files)

if !isempty(mc_labels)
    println("sampled p-values:")
    for label in mc_labels
        idx = only(find_reaction_indices(result.network, label))
        println("  ", label, " => p=", result.rate_p_values[idx])
    end
end

# Final-state snapshot (active + inert species, plus the final time/T9) so
# decay_ppn.jl can pick up from here without re-running the network solve.
profiles = trajectory_profiles(result.trajectory)
final_time = last(result.times)
final_T9 = profiles.T9(final_time)
final_all = merge(result.inert_mass_fractions, result.final_mass_fractions)

final_state_path = joinpath(output_dir, "final_state.csv")
open(final_state_path, "w") do io
    println(io, "# final_time_s,", final_time)
    println(io, "# final_T9,", final_T9)
    println(io, "species,mass_fraction")
    for (name, value) in sort(collect(final_all); by=first)
        println(io, name, ",", value)
    end
end
println("final state               = ", final_state_path)
