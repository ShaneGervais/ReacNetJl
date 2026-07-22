#=
Apply post-processing weak decay to a run_ppn.jl final-state snapshot for a
given decay time, without re-running the network solve.

Usage:
    julia --project=. decay_ppn.jl <run_output_dir> <decay_time_seconds> [output_dir] [--option 1|2|3]

Reads <run_output_dir>/final_state.csv (written by run_ppn.jl) and writes
decayed_state_<t>s.csv into [output_dir] (defaults to <run_output_dir> if
not given, so the decayed state lands next to the run it came from unless
you ask for somewhere else).

By default, uses the same rate-library --option the original run_ppn.jl
call used (saved in final_state.csv) for the weak-decay tables, so a run
made with --option 3 doesn't get silently decayed with --option 1's tables.
Pass --option explicitly to override.
=#

using ReacNetJl
using Printf

include(joinpath(@__DIR__, "rate_options.jl"))

function parse_args(args)
    positional = String[]
    option = nothing
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--option"
            option = parse(Int, args[i+1])
            i += 2
        else
            push!(positional, arg)
            i += 1
        end
    end
    return positional, option
end

const USAGE = """
Usage: julia --project=. decay_ppn.jl <run_output_dir> <decay_time_seconds> [output_dir] [--option 1|2|3]
"""

positional, option_override = parse_args(ARGS)
length(positional) >= 2 || error(USAGE)

run_output_dir = positional[1]
decay_time_s = parse(Float64, positional[2])
decay_time_s >= 0.0 || error("decay_time_seconds must be non-negative")
output_dir = length(positional) >= 3 ? positional[3] : run_output_dir

final_state_path = joinpath(run_output_dir, "final_state.csv")
isfile(final_state_path) || error("no final_state.csv found in $run_output_dir; run run_ppn.jl first")

final_T9 = 0.0
final_time = 0.0
saved_option = nothing
final_X = Dict{String,Float64}()
for line in eachline(final_state_path)
    if startswith(line, "# final_time_s,")
        global final_time = parse(Float64, split(line, ",")[2])
    elseif startswith(line, "# final_T9,")
        global final_T9 = parse(Float64, split(line, ",")[2])
    elseif startswith(line, "# option,")
        global saved_option = parse(Int, split(line, ",")[2])
    elseif startswith(line, "species,") || isempty(strip(line))
        continue
    else
        name, value = split(line, ",")
        final_X[name] = parse(Float64, value)
    end
end

option = something(option_override, saved_option, 1)

println("read final state from ", final_state_path)
println("final time = ", final_time, " s   final T9 = ", final_T9)
println("option     = ", option, "  (", OPTION_DESCRIPTIONS[option], ")",
    option_override === nothing ? " [from original run]" : " [override]")
println("decaying for ", decay_time_s, " s")

tables = rate_tables_for_option(option)
result = decay_mass_fractions(tables, final_X, decay_time_s; T9=final_T9)

println("weak reactions used = ", length(result.decay_tables))

mkpath(output_dir)
decay_path = joinpath(output_dir, @sprintf("decayed_state_%.0fs.csv", decay_time_s))
open(decay_path, "w") do io
    println(io, "species,mass_fraction")
    for (name, value) in sort(collect(result.mass_fractions); by=first)
        println(io, name, ",", value)
    end
end
println("wrote ", decay_path)
