#=
Run ReacNetJl's single-zone post-processing network on a trajectory and
initial-abundance file, writing the standard CSV outputs (mass_fractions.csv,
reaction_fluxes.csv, integrated_fluxes.csv, network.csv) plus a
final_state.csv snapshot that decay_ppn.jl can later decay without
re-running the solve.

Usage:
    julia --project=. run_ppn.jl <trajectory_path> <abundance_path> [output_dir] \
        [--option 1|2|3] [--screening weak|chugunov|none] \
        [--factor "label=value"]... [--mc-factor "label"]... [--mc-all] [--seed N]

Defaults: output_dir="outputs/run_ppn", option=1, screening=chugunov

Rate library options (see rate_options.jl for the full rationale):
    1 = Iliadis 2002: NACRE (A<20) + Iliadis 2001 (A=20-40)   [default]
    2 = STARLIB v6.10 (full library)
    3 = STARLIB v6.10 + etr25 (2025) targeted update

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

    # Whole-network Monte Carlo trial: every reaction with real STARLIB/NACRE
    # uncertainty data gets its own p ~ Normal(0,1); everything else (no
    # uncertainty digitized) stays nominal. One trial per --seed value --
    # see SensitivityStudy/run_monte_carlo.jl for the multi-trial driver.
    julia --project=. run_ppn.jl trajectory.input initial_abundance.dat \
        --option 3 --mc-all --seed 1

    # Use the most up-to-date rate library instead of the Iliadis-2002 baseline:
    julia --project=. run_ppn.jl trajectory.input initial_abundance.dat --option 3
=#

using ReacNetJl
using Random

include(joinpath(@__DIR__, "rate_options.jl"))

function parse_args(args)
    positional = String[]
    option = 1
    screening = :chugunov
    factors = Dict{String,Float64}()
    mc_labels = String[]
    mc_all = false
    seed = nothing
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--option"
            option = parse(Int, args[i+1])
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
        elseif arg == "--mc-all"
            mc_all = true
            i += 1
        elseif arg == "--seed"
            seed = parse(Int, args[i+1])
            i += 2
        else
            push!(positional, arg)
            i += 1
        end
    end
    return positional, option, screening, factors, mc_labels, mc_all, seed
end

const USAGE = """
Usage: julia --project=. run_ppn.jl <trajectory_path> <abundance_path> [output_dir] \\
    [--option 1|2|3] [--screening weak|chugunov|none] \\
    [--factor "label=value"]... [--mc-factor "label"]... [--mc-all] [--seed N]
"""

positional, option, screening, factors, mc_labels, mc_all, seed = parse_args(ARGS)
length(positional) >= 2 || error(USAGE)
mc_all && !isempty(mc_labels) && error("--mc-all and --mc-factor are mutually exclusive")

trajectory_path = positional[1]
abundance_path = positional[2]
output_dir = length(positional) >= 3 ? positional[3] : "outputs/run_ppn"
rng = seed === nothing ? Random.default_rng() : MersenneTwister(seed)

println("trajectory = ", trajectory_path)
println("abundances = ", abundance_path)
println("option     = ", option, "  (", OPTION_DESCRIPTIONS[option], ")")
println("screening  = ", screening === nothing ? "none" : screening)
println("output_dir = ", output_dir)
isempty(factors) || println("factors    = ", factors)
isempty(mc_labels) || println("mc-factor  = ", mc_labels, "  (seed=", seed === nothing ? "none (non-reproducible)" : seed, ")")
mc_all && println("mc-all     = true  (seed=", seed === nothing ? "none (non-reproducible)" : seed, ")")
println()

# Option 1 goes through the rates=:iliadis2002 path so run_ppn still builds
# and returns a rate_policy_report; options 2/3 pass a prebuilt table list
# directly (rate_options.jl), which skips that label-selection reporting
# since it doesn't apply to a plain STARLIB library.
common_kwargs = (
    screening=screening,
    output_dir=output_dir,
    rate_factors=isempty(factors) ? nothing : factors,
    rate_sample_labels=isempty(mc_labels) ? nothing : mc_labels,
    rate_sample_all=mc_all,
    rng=rng,
)
result = if option == 1
    run_ppn(trajectory_path, abundance_path; rates=:iliadis2002, common_kwargs...)
else
    run_ppn(trajectory_path, abundance_path; tables=rate_tables_for_option(option), common_kwargs...)
end

println("validated reactions      = ", result.validation.num_reactions)
println("network species          = ", length(result.network.species))
println("solver accepted/rejected = ", result.solver_stats.accepted_steps, " / ", result.solver_stats.rejected_steps)
println("mass drift               = ", result.mass_fraction_drift.initial, " -> ", result.mass_fraction_drift.final)
println("output files              = ", result.output_files)

if mc_all && result.rate_sample_report !== nothing
    n_sampled = length(result.rate_sample_report.sampled_reactions)
    n_unsampled = length(result.rate_sample_report.unsampled_reactions)
    println("mc-all: ", n_sampled, " reaction(s) sampled from real rate uncertainty; ",
            n_unsampled, " reaction(s) held at nominal rate for lack of it")
    println("WARNING: reactions without uncertainty data (kept nominal, not perturbed this trial):")
    for label in result.rate_sample_report.unsampled_reactions
        println("  - ", label)
    end
end

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
    println(io, "# option,", option)
    println(io, "species,mass_fraction")
    for (name, value) in sort(collect(final_all); by=first)
        println(io, name, ",", value)
    end
end
println("final state               = ", final_state_path)
