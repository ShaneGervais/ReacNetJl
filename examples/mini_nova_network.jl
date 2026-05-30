using ReacNetJl

#=
    Mini nova single-zone post-processing example

This example uses `trajectory.input` from the project root when present, falling
back to `examples/fake_nova_trajectory.dat` for tests and demos. It follows an
expanded CNO/NeNa/MgAl-to-Ca reaction list from STARLIB. It is intended as a
workflow and performance/diagnostic example, not a validated nova model.

Monte Carlo uncertainty sampling is intentionally not used here. The goal is to
first verify deterministic single-zone post-processing on a trajectory.
=#

labels = [
    # Hot CNO-ish flow and decays
    "12C(p,γ)13N",
    "13N(β+)13C",
    "13C(p,γ)14N",
    "14N(p,γ)15O",
    "15O(β+)15N",
    "15N(p,α)12C",
    "15N(p,γ)16O",
    "16O(p,γ)17F",
    "17F(β+)17O",
    "17O(p,α)14N",
    "17O(p,γ)18F",
    "18F(β+)18O",
    "18O(p,α)15N",
    "18O(p,γ)19F",
    "18F(p,α)15O",
    "18F(p,γ)19Ne",
    "19Ne(β+)19F",
    "19F(p,α)16O",
    "19F(p,γ)20Ne",

    # NeNa cycle fragments
    "20Ne(p,γ)21Na",
    "21Na(β+)21Ne",
    "21Na(p,γ)22Mg",
    "22Mg(β+)22Na",
    "21Ne(p,γ)22Na",
    "22Na(β+)22Ne",
    "22Na(p,γ)23Mg",
    "23Mg(β+)23Na",
    "22Ne(p,γ)23Na",
    "23Na(p,α)20Ne",
    "23Na(p,γ)24Mg",

    # MgAl fragments
    "24Mg(p,γ)25Al",
    "25Al(β+)25Mg",
    "25Al(p,γ)26Si",
    "26Si(β+)26Al",
    "26Si(β+)26Al*",
    "25Mg(p,γ)26Al",
    "26Mg(p,γ)27Al",
    "26Al*(β+)26Mg",
    "26Al(p,γ)27Si",
    "26Al*(p,γ)27Si",
    "27Si(β+)27Al",
    "27Al(p,α)24Mg",
    "27Al(p,γ)28Si",

    # Si-Ca extension
    "28Si(p,γ)29P",
    "29P(β+)29Si",
    "29Si(p,γ)30P",
    "30P(β+)30Si",
    "30P(p,γ)31S",
    "31S(β+)31P",
    "30Si(p,γ)31P",
    "31P(p,γ)32S",
    "32S(p,γ)33Cl",
    "33Cl(β+)33S",
    "33S(p,γ)34Cl",
    "34Cl(β+)34S",
    "34S(p,γ)35Cl",
    "35Cl(p,γ)36Ar",
    "36Ar(p,γ)37K",
    "37K(β+)37Ar",
    "37Ar(p,γ)38K",
    "38K(β+)38Ar",
    "38Ar(p,γ)39K",
    "39K(p,γ)40Ca",
]

tables = read_starlib()
network = network_from_labels(tables, labels)
validation = network_validation_report(network; throw_on_error=true)

project_trajectory_path = joinpath(dirname(@__DIR__), "trajectory.input")
trajectory_path = isfile(project_trajectory_path) ? project_trajectory_path : joinpath(@__DIR__, "fake_nova_trajectory.dat")
trajectory = read_trajectory(trajectory_path)
profiles = trajectory_profiles(trajectory)
trajectory_duration = last(trajectory.time) - first(trajectory.time)
dt_initial = trajectory_duration > 100.0 ? 1.0 : 0.02
dt_min = trajectory_duration > 100.0 ? 1.0e-4 : 1.0e-6
dt_max = trajectory_duration > 100.0 ? 20.0 : 0.05
convergence_dt_limits = trajectory_duration > 100.0 ? (20.0, 10.0, 5.0) : (0.05, 0.025, 0.01)
screening_model = :weak

X0 = Dict(
    "p" => 0.55,
    "he4" => 0.33,
    "12C" => 0.025,
    "14N" => 0.010,
    "16O" => 0.035,
    "20Ne" => 0.030,
    "22Ne" => 0.010,
    "24Mg" => 0.007,
    "25Mg" => 0.002,
    "26Mg" => 0.001,
)
Y0 = abundances_from_mass_fractions(network, X0; check_sum=true, atol=1.0e-12)

times, history, solver_stats = solve_network_adaptive(
    network,
    Y0,
    (first(trajectory.time), last(trajectory.time)),
    dt_initial,
    profiles.rho,
    profiles.T9;
    method=:backward_euler,
    max_fractional_change=0.50,
    max_absolute_change=1.0e-4,
    abundance_floor=1.0e-8,
    dt_min=dt_min,
    dt_max=dt_max,
    screening=screening_model,
    return_stats=true,
)

flux_history = reaction_flux_history(network, history, times, profiles.rho, profiles.T9; screening=screening_model)
total_fluxes = integrated_fluxes(times, flux_history)
energy_history = energy_generation_history(network, history, times, profiles.rho, profiles.T9; screening=screening_model)
total_energy = integrated_energy_generation(times, energy_history)
initial_X = mass_fractions_from_abundances(network, history[1, :])
final_X = mass_fractions_from_abundances(network, history[end, :])
mass_drift = mass_fraction_drift(network, history)
positivity = abundance_diagnostics(network, history)

println("Mini nova trajectory network")
println("trajectory = ", trajectory_path)
println("validated reactions = ", validation.num_reactions)
println("species = ", length(network.species))
println("screening = ", screening_model)
println("accepted timesteps = ", solver_stats.accepted_steps)
println("rejected timesteps = ", solver_stats.rejected_steps)
println("dt range = ", solver_stats.min_dt, " to ", solver_stats.max_dt, " s")
println("Newton iterations mean/max = ", solver_stats.mean_newton_iterations, " / ", solver_stats.max_newton_iterations)
println("Newton failed steps = ", solver_stats.newton_failed_steps)
println("time = ", first(times), " to ", last(times), " s")
println("total mass fraction = ", mass_drift.initial, " to ", mass_drift.final)
println("peak epsilon_nuc = ", maximum(energy_history), " erg g^-1 s^-1")
println("integrated nuclear energy = ", total_energy, " erg g^-1")
if mass_drift.max_abs_drift > 1.0e-3
    println("WARNING: total mass-fraction drift exceeded 1.0e-3; treat this run as numerically suspect.")
end
println(
    "minimum abundance = ",
    positivity.min_abundance,
    " (",
    positivity.min_abundance_species,
    ", row ",
    positivity.min_abundance_time_index,
    ")",
)
println(
    "minimum mass fraction = ",
    positivity.min_mass_fraction,
    " (",
    positivity.min_mass_fraction_species,
    ", row ",
    positivity.min_mass_fraction_time_index,
    ")",
)
if positivity.min_abundance < -1.0e-20
    println("WARNING: negative abundances below tolerance were produced.")
end
println()

println("Selected final mass fractions")
for name in ["p", "he4", "c12", "n14", "o15", "o16", "f18", "ne20", "na22", "na23", "mg24", "al26", "si28"]
    if haskey(final_X, name)
        println(rpad(name, 8), " initial=", get(initial_X, name, 0.0), " final=", final_X[name])
    end
end

println()
println("Top integrated reaction fluxes")
order = sortperm(total_fluxes; rev=true)
for i in first(order, min(10, length(order)))
    println(rpad(reaction_string(network.reactions[i]), 20), "  ", total_fluxes[i])
end

println()
println("Backward Euler timestep convergence")
println("dt_max | steps | rejected | mean Newton | final X(c12) | final X(o16) | final X(si28) | mass drift | top flux")

function convergence_summary(dt_limit)
    case_times, case_history, case_stats = solve_network_adaptive(
        network,
        Y0,
        (first(trajectory.time), last(trajectory.time)),
        min(dt_initial, dt_limit),
        profiles.rho,
        profiles.T9;
        method=:backward_euler,
        max_fractional_change=0.50,
        max_absolute_change=1.0e-4,
        abundance_floor=1.0e-8,
        dt_min=dt_min,
        dt_max=dt_limit,
        screening=screening_model,
        return_stats=true,
    )
    case_flux_history = reaction_flux_history(network, case_history, case_times, profiles.rho, profiles.T9; screening=screening_model)
    case_total_fluxes = integrated_fluxes(case_times, case_flux_history)
    case_final_X = mass_fractions_from_abundances(network, case_history[end, :])
    case_drift = mass_fraction_drift(network, case_history)
    top_index = argmax(case_total_fluxes)

    return (
        stats=case_stats,
        final_X=case_final_X,
        drift=case_drift,
        top_reaction=reaction_string(network.reactions[top_index]),
        top_flux=case_total_fluxes[top_index],
    )
end

summaries = Dict(
    dt_max => (
        stats=solver_stats,
        final_X=final_X,
        drift=mass_drift,
        top_reaction=reaction_string(network.reactions[argmax(total_fluxes)]),
        top_flux=maximum(total_fluxes),
    ),
)

for dt_limit in convergence_dt_limits[2:end]
    summaries[dt_limit] = convergence_summary(dt_limit)
end

for dt_limit in convergence_dt_limits
    summary = summaries[dt_limit]
    println(join(
        [
            string(dt_limit),
            string(summary.stats.accepted_steps),
            string(summary.stats.rejected_steps),
            string(round(summary.stats.mean_newton_iterations; digits=2)),
            string(get(summary.final_X, "c12", 0.0)),
            string(get(summary.final_X, "o16", 0.0)),
            string(get(summary.final_X, "si28", 0.0)),
            string(summary.drift.drift),
            summary.top_reaction * " = " * string(summary.top_flux),
        ],
        " | ",
    ))
end
