# Benchmark the backward-Euler solver: analytic Jacobian + per-step rate cache
# against the finite-difference Jacobian, on the real H-Ca nova network.
#
# Run from the project root (or anywhere):
#   julia --project=. examples/benchmark_solver.jl
#
# The integration covers the first 1000 s of `trajectory.input` so a full
# comparison stays affordable; both runs solve the identical problem and the
# final mass fractions are compared at the end.

if Base.find_package("ReacNetJl") === nothing
    import Pkg
    Pkg.activate(dirname(@__DIR__); io=devnull)
end

using Printf
using ReacNetJl

project_root = dirname(@__DIR__)
trajectory = read_trajectory(joinpath(project_root, "trajectory.input"))
profiles = trajectory_profiles(trajectory)
abundance_path = isfile(joinpath(project_root, "initial_abundance.DAT")) ?
    joinpath(project_root, "initial_abundance.DAT") : joinpath(project_root, "initial_abundance.dat")
X_file = read_initial_abundances(abundance_path; normalize=true)

tables = read_starlib()
forward_tables = select_h_ca_reaction_tables(tables, keys(X_file); projectiles=("p", "he4", "he3", "d", "n"))
reverse_summary = add_reverse_reaction_tables(tables, forward_tables; generate_detailed_balance=true)
network = network_from_tables(reverse_summary.tables)
@printf("network: %d species, %d reactions\n", length(network.species), length(network.reactions))

X0 = Dict(name => value for (name, value) in X_file if haskey(network.species_index, name))
Y0 = abundances_from_mass_fractions(network, X0)

function run_once(jacobian, duration)
    tspan = (first(trajectory.time), first(trajectory.time) + duration)
    return solve_network_adaptive(
        network, Y0, tspan, 1.0, profiles.rho, profiles.T9;
        method=:backward_euler, screening=:weak, jacobian=jacobian,
        max_fractional_change=0.50, max_absolute_change=1.0e-4,
        abundance_floor=1.0e-8, max_newton_iterations=80,
        dt_min=1.0e-8, dt_max=20.0, return_stats=true,
    )
end

# The comparison leg is short: the finite-difference Jacobian is the slow
# path being replaced, and 200 s of trajectory is enough to measure it.
comparison_duration = 200.0

println("warm-up (compilation, ~1 min)...")
run_once(:analytic, 10.0)
run_once(:finite_difference, 10.0)

@printf("timing analytic Jacobian on %.0f s of trajectory...\n", comparison_duration)
t_analytic = @elapsed times_a, history_a, stats_a = run_once(:analytic, comparison_duration)
@printf("  analytic:          %8.2f s  (%d accepted steps, %.2f mean Newton iterations)\n",
    t_analytic, stats_a.accepted_steps, stats_a.mean_newton_iterations)

@printf("timing finite-difference Jacobian on the same %.0f s (expect minutes)...\n", comparison_duration)
t_fd = @elapsed times_f, history_f, stats_f = run_once(:finite_difference, comparison_duration)
@printf("  finite difference: %8.2f s  (%d accepted steps, %.2f mean Newton iterations)\n",
    t_fd, stats_f.accepted_steps, stats_f.mean_newton_iterations)

@printf("speedup: %.1fx\n", t_fd / t_analytic)

# Agreement metric over physically meaningful abundances only: species below
# the Newton tolerance floor are noise in both runs and carry no information.
X_a = mass_fractions_from_abundances(network, history_a[end, :])
X_f = mass_fractions_from_abundances(network, history_f[end, :])
worst_relative = 0.0
worst_absolute = 0.0
for (name, value) in X_a
    reference = X_f[name]
    global worst_absolute = max(worst_absolute, abs(value - reference))
    if max(abs(value), abs(reference)) > 1.0e-12
        global worst_relative = max(worst_relative, abs(value - reference) / max(abs(reference), 1.0e-12))
    end
end
@printf("final mass fractions: worst relative difference (X > 1e-12) = %.3e, worst absolute = %.3e\n",
    worst_relative, worst_absolute)

println("timing analytic Jacobian on 1000 s of trajectory (headline run)...")
t_full = @elapsed times_full, history_full, stats_full = run_once(:analytic, 1000.0)
@printf("  analytic, 1000 s:  %8.2f s  (%d accepted steps)\n", t_full, stats_full.accepted_steps)
