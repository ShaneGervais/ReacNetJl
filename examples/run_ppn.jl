#=
Run ReacNetJl's single-zone post-processing network on a trajectory and
initial-abundance file, writing the standard CSV outputs (mass_fractions.csv,
reaction_fluxes.csv, network.csv) plus a final_state.csv snapshot that
decay_ppn.jl can later decay without re-running the solve.

Usage:
    julia --project=. run_ppn.jl <trajectory_path> <abundance_path> [output_dir] \
        [--rates starlib|iliadis2002] [--screening weak|chugunov|none]

Defaults: output_dir="outputs/run_ppn", rates=iliadis2002, screening=chugunov
=#

using ReacNetJl

function parse_args(args)
    positional = String[]
    rates = :iliadis2002
    screening = :chugunov
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
        else
            push!(positional, arg)
            i += 1
        end
    end
    return positional, rates, screening
end

const USAGE = """
Usage: julia --project=. run_ppn.jl <trajectory_path> <abundance_path> [output_dir] \\
    [--rates starlib|iliadis2002] [--screening weak|chugunov|none]
"""

positional, rates, screening = parse_args(ARGS)
length(positional) >= 2 || error(USAGE)

trajectory_path = positional[1]
abundance_path = positional[2]
output_dir = length(positional) >= 3 ? positional[3] : "outputs/run_ppn"

println("trajectory = ", trajectory_path)
println("abundances = ", abundance_path)
println("rates      = ", rates)
println("screening  = ", screening === nothing ? "none" : screening)
println("output_dir = ", output_dir)
println()

result = run_ppn(trajectory_path, abundance_path; rates=rates, screening=screening, output_dir=output_dir)

println("validated reactions      = ", result.validation.num_reactions)
println("network species          = ", length(result.network.species))
println("solver accepted/rejected = ", result.solver_stats.accepted_steps, " / ", result.solver_stats.rejected_steps)
println("mass drift               = ", result.mass_fraction_drift.initial, " -> ", result.mass_fraction_drift.final)
println("output files              = ", result.output_files)

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
