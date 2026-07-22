#=
Apply post-processing weak decay to a run_ppn.jl final-state snapshot for a
given decay time, without re-running the network solve.

Usage:
    julia --project=. decay_ppn.jl <run_output_dir> <decay_time_seconds> [--rates starlib|iliadis2002]

Reads <run_output_dir>/final_state.csv (written by run_ppn.jl) and writes
<run_output_dir>/decayed_state_<t>s.csv.
=#

using ReacNetJl
using Printf

function parse_args(args)
    positional = String[]
    rates = :iliadis2002
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--rates"
            rates = Symbol(args[i+1])
            i += 2
        else
            push!(positional, arg)
            i += 1
        end
    end
    return positional, rates
end

const USAGE = """
Usage: julia --project=. decay_ppn.jl <run_output_dir> <decay_time_seconds> [--rates starlib|iliadis2002]
"""

positional, rates = parse_args(ARGS)
length(positional) >= 2 || error(USAGE)

run_output_dir = positional[1]
decay_time_s = parse(Float64, positional[2])
decay_time_s >= 0.0 || error("decay_time_seconds must be non-negative")

final_state_path = joinpath(run_output_dir, "final_state.csv")
isfile(final_state_path) || error("no final_state.csv found in $run_output_dir; run run_ppn.jl first")

final_T9 = 0.0
final_time = 0.0
final_X = Dict{String,Float64}()
for line in eachline(final_state_path)
    if startswith(line, "# final_time_s,")
        global final_time = parse(Float64, split(line, ",")[2])
    elseif startswith(line, "# final_T9,")
        global final_T9 = parse(Float64, split(line, ",")[2])
    elseif startswith(line, "species,") || isempty(strip(line))
        continue
    else
        name, value = split(line, ",")
        final_X[name] = parse(Float64, value)
    end
end

println("read final state from ", final_state_path)
println("final time = ", final_time, " s   final T9 = ", final_T9)
println("decaying for ", decay_time_s, " s")

tables = rates == :iliadis2002 ? iliadis2002_rate_tables().tables : read_starlib()
result = decay_mass_fractions(tables, final_X, decay_time_s; T9=final_T9)

println("weak reactions used = ", length(result.decay_tables))

decay_path = joinpath(run_output_dir, @sprintf("decayed_state_%.0fs.csv", decay_time_s))
open(decay_path, "w") do io
    println(io, "species,mass_fraction")
    for (name, value) in sort(collect(result.mass_fractions); by=first)
        println(io, name, ",", value)
    end
end
println("wrote ", decay_path)
