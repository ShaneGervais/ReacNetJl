using ReacNetJl

#=
    Single-zone 18F destruction example

This example builds a tiny nova-relevant network from STARLIB data:

    18F(p,α)15O
    18F(p,γ)19Ne

It evolves abundances at fixed temperature and density using the package's
built-in RK4 solver. This is a learning/example calculation, not yet a
production nova network.
=#

labels = [
    "18F(p,α)15O",
    "18F(p,γ)19Ne",
]

tables = read_starlib()
network = network_from_labels(tables, labels)
validation = network_validation_report(network; throw_on_error=true)

X0 = Dict(
    "p" => 0.70,
    "he4" => 0.28,
    "18F" => 1.0e-5,
    "15O" => 0.0,
    "19Ne" => 0.0,
)

Y0 = abundances_from_mass_fractions(network, X0)

rho = 1.0e3    # g cm^-3
T9 = 0.2       # GK
tspan = (0.0, 1.0e-3)
dt = 1.0e-5

times, history = solve_network(network, Y0, tspan, dt, rho, T9; method=:rk4)
flux_history = reaction_flux_history(network, history, times, rho, T9)
total_fluxes = integrated_fluxes(times, flux_history)

initial_X = mass_fractions_from_abundances(network, history[1, :])
final_X = mass_fractions_from_abundances(network, history[end, :])

println("Single-zone 18F example")
println("T9 = ", T9, " GK")
println("rho = ", rho, " g cm^-3")
println("time = ", first(times), " to ", last(times), " s")
println("validated reactions = ", validation.num_reactions)
println()
println("Species      initial X        final X")
for name in network.species
    println(rpad(name, 8), "  ", initial_X[name], "  ", final_X[name])
end

println()
println("Integrated reaction fluxes")
for (i, reaction) in pairs(network.reactions)
    println(rpad(reaction_string(reaction), 20), "  ", total_fluxes[i])
end
