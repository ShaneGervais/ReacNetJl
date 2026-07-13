# Run a nova case end to end: read inputs, solve, write CSVs. No CLI parsing,
# no diagnostics beyond a one-line summary — for that, use
# single_zone_nova_ppn.jl. This is the "just run it" path.
#
# Usage:
#   julia --project=. examples/run_case.jl [case_name]
#
# `case_name` selects nova_cases/<case_name>/ (default: nova_test), which
# must contain trajectory.input and initial_abundance.dat (or .DAT). Output
# goes to outputs/<case_name>/: mass_fractions.csv, reaction_fluxes.csv, and
# network.csv.

if Base.find_package("ReacNetJl") === nothing
    import Pkg
    Pkg.activate(dirname(@__DIR__); io=devnull)
end

using ReacNetJl

project_root = dirname(@__DIR__)
case = length(ARGS) >= 1 ? ARGS[1] : "nova_test"
case_dir = joinpath(project_root, "nova_cases", case)
abundance_path = isfile(joinpath(case_dir, "initial_abundance.DAT")) ?
    joinpath(case_dir, "initial_abundance.DAT") : joinpath(case_dir, "initial_abundance.dat")
output_dir = joinpath(project_root, "outputs", case)

result = run_ppn(
    joinpath(case_dir, "trajectory.input"),
    abundance_path;
    rates=:iliadis2002,
    screening=:chugunov,
    output_dir=output_dir,
)

println("case             = ", case)
println("species/reactions = ", length(result.network.species), " / ", length(result.network.reactions))
println("solver steps     = ", result.solver_stats.accepted_steps, " accepted, ", result.solver_stats.rejected_steps, " rejected")
println("mass drift       = ", result.mass_fraction_drift.max_abs_drift)
println("wrote            = ", output_dir)
