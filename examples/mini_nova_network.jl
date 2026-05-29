using ReacNetJl

#=
    Mini nova single-zone post-processing example

This example uses a fake-but-realistic nova thermodynamic trajectory and a
moderate CNO/NeNa/MgAl reaction list from STARLIB. It is intended as a workflow
and performance/diagnostic example, not a validated nova model.

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

    # NeNa cycle fragments
    "20Ne(p,γ)21Na",
    "21Na(β+)21Ne",
    "21Na(p,γ)22Mg",
    "21Ne(p,γ)22Na",
    "22Na(β+)22Ne",
    "22Na(p,γ)23Mg",
    "22Ne(p,γ)23Na",
    "23Na(p,α)20Ne",
    "23Na(p,γ)24Mg",

    # MgAl fragments
    "24Mg(p,γ)25Al",
    "25Al(β+)25Mg",
    "25Al(p,γ)26Si",
    "25Mg(p,γ)26Al",
    "26Mg(p,γ)27Al",
    "26Al(p,γ)27Si",
    "27Al(p,α)24Mg",
    "27Al(p,γ)28Si",
]

tables = read_starlib()
network = network_from_labels(tables, labels)
validation = network_validation_report(network; throw_on_error=true)

trajectory_path = joinpath(@__DIR__, "fake_nova_trajectory.dat")
trajectory = read_trajectory(trajectory_path)
profiles = trajectory_profiles(trajectory)

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

times, history = solve_network_adaptive(
    network,
    Y0,
    (first(trajectory.time), last(trajectory.time)),
    0.02,
    profiles.rho,
    profiles.T9;
    method=:rk4,
    max_fractional_change=0.10,
    abundance_floor=1.0e-24,
    dt_min=1.0e-8,
    dt_max=0.05,
)

flux_history = reaction_flux_history(network, history, times, profiles.rho, profiles.T9)
total_fluxes = integrated_fluxes(times, flux_history)
initial_X = mass_fractions_from_abundances(network, history[1, :])
final_X = mass_fractions_from_abundances(network, history[end, :])

println("Mini nova trajectory network")
println("trajectory = ", trajectory_path)
println("validated reactions = ", validation.num_reactions)
println("species = ", length(network.species))
println("accepted timesteps = ", length(times) - 1)
println("time = ", first(times), " to ", last(times), " s")
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
