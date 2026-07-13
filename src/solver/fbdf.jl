# solve_network_fbdf hook, implemented by the OrdinaryDiffEqBDF package extension.

"""
    solve_network_fbdf(network, Y0, tspan, rho, T9; kwargs...)

Solve the network with the variable-order stiff `FBDF` integrator from
OrdinaryDiffEqBDF, using ReacNetJl's cached RHS and analytic Jacobian.
Higher-order than backward Euler, so it takes far fewer steps at the same
accuracy.

This function is provided by a package extension: install and load the solver
package first with `import Pkg; Pkg.add("OrdinaryDiffEqBDF")` and
`using OrdinaryDiffEqBDF`. Returns `(times, history, stats)` like
`solve_network_adaptive`. Also available through `run_ppn(...; method=:fbdf)`.
"""
function solve_network_fbdf(network::ReactionNetwork, Y0, tspan, rho, T9; kwargs...)
    error("solve_network_fbdf requires the OrdinaryDiffEqBDF extension; run `import Pkg; Pkg.add(\"OrdinaryDiffEqBDF\")` and add `using OrdinaryDiffEqBDF` before calling it")
end
