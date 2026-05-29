using ReacNetJl

#=
    Small STARLIB uncertainty Monte Carlo example

This runs repeated one-zone calculations for the two-channel 18F destruction
network. Each reaction gets one sampled STARLIB uncertainty parameter p per run:

    rate_sampled(T) = rate_recommended(T) * factor_uncertainty(T)^p

This example prints a simple final 18F mass-fraction range. Larger analysis and
plotting should be done in external scripts/notebooks.
=#

labels = [
    "18F(p,α)15O",
    "18F(p,γ)19Ne",
]

tables = read_starlib()
network = network_from_labels(tables, labels)
network_validation_report(network; throw_on_error=true)

X0 = Dict(
    "p" => 0.70,
    "he4" => 0.28,
    "18F" => 1.0e-5,
    "15O" => 0.0,
    "19Ne" => 0.0,
)
Y0 = abundances_from_mass_fractions(network, X0)

rho = 1.0e3
T9 = 0.2
tspan = (0.0, 1.0e-3)
dt = 1.0e-5
nruns = 20

result = run_monte_carlo(
    network,
    Y0,
    tspan,
    dt,
    rho,
    T9;
    nruns=nruns,
    seed=2026,
)

f18_index = findfirst(==("f18"), result.species)
final_X_f18 = [
    mass_fraction_from_abundance(result.final_abundances[i, f18_index], 18)
    for i in 1:nruns
]

println("18F Monte Carlo example")
println("runs = ", nruns)
println("reactions = ", join(result.reactions, ", "))
println("final X_f18 min = ", minimum(final_X_f18))
println("final X_f18 max = ", maximum(final_X_f18))
println("final X_f18 mean = ", sum(final_X_f18) / length(final_X_f18))
