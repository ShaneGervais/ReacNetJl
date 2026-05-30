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

X0 = Dict(
    "p" => 0.70,
    "he4" => 0.28,
    "18F" => 1.0e-5,
    "15O" => 0.0,
    "19Ne" => 0.0,
)

rho = 1.0e3    # g cm^-3
T9 = 0.2       # GK
tspan = (0.0, 1.0e-3)
dt = 1.0e-5

result = solve_single_zone(tables, labels, X0, tspan, dt, rho, T9; adaptive=false)

println("Single-zone 18F example")
println("T9 = ", T9, " GK")
println("rho = ", rho, " g cm^-3")
println("time = ", first(result.times), " to ", last(result.times), " s")
println("validated reactions = ", result.validation.num_reactions)
println()
println("Species      initial X        final X")
for name in result.network.species
    println(rpad(name, 8), "  ", result.initial_mass_fractions[name], "  ", result.final_mass_fractions[name])
end

println()
println("Integrated reaction fluxes")
for (i, reaction) in pairs(result.network.reactions)
    println(rpad(reaction_string(reaction), 20), "  ", result.integrated_fluxes[i])
end
