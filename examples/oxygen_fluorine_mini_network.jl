using ReacNetJl

#=
    Oxygen-fluorine mini-network example

This example demonstrates a slightly larger coupled network:

    17O(p,α)14N
    17O(p,γ)18F
    18F(p,α)15O
    18F(p,γ)19Ne

The purpose is to inspect how flux moves from 17O into 18F and then into its
proton-destruction channels. Plotting is intentionally left to external scripts;
this example prints calculation-ready abundance and flux diagnostics.
=#

labels = [
    "17O(p,α)14N",
    "17O(p,γ)18F",
    "18F(p,α)15O",
    "18F(p,γ)19Ne",
]

tables = read_starlib()
network = network_from_labels(tables, labels)
validation = network_validation_report(network; throw_on_error=true)

X0 = Dict(
    "p" => 0.70,
    "he4" => 0.28,
    "17O" => 1.0e-4,
    "14N" => 0.0,
    "18F" => 0.0,
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
final_balance = species_flux_balance(network, history[end, :], rho, T9)

println("Oxygen-fluorine mini-network example")
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
order = sortperm(total_fluxes; rev=true)
for i in order
    println(rpad(reaction_string(network.reactions[i]), 20), "  ", total_fluxes[i])
end

println()
println("Final species flux balance in abundance units")
println("Species      production       destruction      net")
for (i, name) in pairs(network.species)
    println(
        rpad(name, 8), "  ",
        final_balance.production[i], "  ",
        final_balance.destruction[i], "  ",
        final_balance.net[i],
    )
end

println()
println("Graph-like reaction edges")
for edge in reaction_edges(network)
    println(edge.reaction, ": ", edge.from, " -> ", edge.to)
end
