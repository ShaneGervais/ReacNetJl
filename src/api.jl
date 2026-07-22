# User-facing one-call APIs tying io, physics, and solver together.

"""
    solve_single_zone(tables, labels, X0, tspan, dt, rho, T9; adaptive=true, ...)

Build and run a single-zone post-processing network from user-facing inputs.

This is the convenience workflow for interactive use:
- select STARLIB reactions by label
- infer and validate the network
- convert initial mass fractions `X0` to abundances
- run fixed-step or adaptive integration
- return mass-fraction diagnostics with the raw abundance history
"""
function solve_single_zone(
    tables::AbstractVector{ReactionRateTable},
    labels::AbstractVector{<:AbstractString},
    X0::AbstractDict,
    tspan::Tuple{<:Real,<:Real},
    dt::Real,
    rho,
    T9;
    species=nothing,
    source=nothing,
    on_multiple::Symbol=:error,
    adaptive::Bool=true,
    method::Symbol=:rk4,
    normalize_mass_fractions::Bool=false,
    check_mass_fraction_sum::Bool=false,
    mass_fraction_atol::Real=1.0e-8,
    validate::Bool=true,
    max_fractional_change::Real=0.05,
    max_absolute_change::Real=Inf,
    abundance_floor::Real=1.0e-30,
    dt_min::Real=1.0e-12,
    dt_max::Real=Inf,
    safety::Real=0.8,
    growth_factor::Real=2.0,
    shrink_factor::Real=0.5,
    max_steps::Integer=1_000_000,
    rate_multipliers=nothing,
    rate_p_values=nothing,
    clamp_negative::Bool=true,
    screening=nothing,
    newton_tolerance::Real=1.0e-10,
    max_newton_iterations::Integer=20,
    finite_difference_epsilon::Real=sqrt(eps(Float64)),
    jacobian::Symbol=:analytic,
)
    network = network_from_labels(tables, labels; species=species, source=source, on_multiple=on_multiple)
    validation = validate ? network_validation_report(network; throw_on_error=true) : network_validation_report(network)
    Y0 = abundances_from_mass_fractions(
        network,
        X0;
        normalize=normalize_mass_fractions,
        check_sum=check_mass_fraction_sum,
        atol=mass_fraction_atol,
    )

    times, history, solver_stats = if adaptive
        solve_network_adaptive(
            network,
            Y0,
            tspan,
            dt,
            rho,
            T9;
            method=method,
            max_fractional_change=max_fractional_change,
            max_absolute_change=max_absolute_change,
            abundance_floor=abundance_floor,
            dt_min=dt_min,
            dt_max=dt_max,
            safety=safety,
            growth_factor=growth_factor,
            shrink_factor=shrink_factor,
            max_steps=max_steps,
            rate_multipliers=rate_multipliers,
            rate_p_values=rate_p_values,
            clamp_negative=clamp_negative,
            screening=screening,
            return_stats=true,
            newton_tolerance=newton_tolerance,
            max_newton_iterations=max_newton_iterations,
            finite_difference_epsilon=finite_difference_epsilon,
            jacobian=jacobian,
        )
    else
        solve_network(
            network,
            Y0,
            tspan,
            dt,
            rho,
            T9;
            method=method,
            rate_multipliers=rate_multipliers,
            rate_p_values=rate_p_values,
            clamp_negative=clamp_negative,
            screening=screening,
            newton_tolerance=newton_tolerance,
            max_newton_iterations=max_newton_iterations,
            finite_difference_epsilon=finite_difference_epsilon,
            jacobian=jacobian,
            return_stats=true,
        )
    end

    flux_history = reaction_flux_history(
        network,
        history,
        times,
        rho,
        T9;
        rate_multipliers=rate_multipliers,
        rate_p_values=rate_p_values,
        screening=screening,
    )
    epsilon_history = energy_generation_history(
        network,
        history,
        times,
        rho,
        T9;
        rate_multipliers=rate_multipliers,
        rate_p_values=rate_p_values,
        screening=screening,
    )

    return (
        network=network,
        validation=validation,
        times=times,
        abundances=history,
        mass_fraction_history=mass_fraction_history(network, history),
        mass_fraction_drift=mass_fraction_drift(network, history),
        abundance_diagnostics=abundance_diagnostics(network, history),
        initial_mass_fractions=mass_fractions_from_abundances(network, view(history, 1, :)),
        final_mass_fractions=mass_fractions_from_abundances(network, view(history, size(history, 1), :)),
        reaction_fluxes=flux_history,
        integrated_fluxes=integrated_fluxes(times, flux_history),
        energy_generation=epsilon_history,
        integrated_energy_generation=integrated_energy_generation(times, epsilon_history),
        solver_stats=solver_stats,
        adaptive=adaptive,
    )
end

"""
    run_ppn(trajectory_path, abundance_path;
            rates=:starlib, screening=:weak, neutron_captures=true, kwargs...)

Run a complete single-zone post-processing nucleosynthesis calculation from a
trajectory file and an initial-abundance file, in one call.

The pipeline is: read the trajectory and abundances, build the H-Ca nova
network from the chosen rate library, add reverse rates, validate the
network, and integrate the abundances over the full trajectory with the
adaptive backward-Euler solver and analytic Jacobian.

# Arguments
- `trajectory_path`: trajectory file with `AGEUNIT`/`TUNIT`/`RHOUNIT` metadata.
- `abundance_path`: initial abundance table (`Z name A X` rows).

# Keywords
- `rates`: `:starlib` (default) or `:iliadis2002` for the NACRE (A < 20) plus
  Iliadis et al. 2001 (A = 20-40) REACLIB baseline of the 2002 nova
  sensitivity study.
- `tables`: pass a prebuilt `Vector{ReactionRateTable}` to skip library
  loading (overrides `rates`).
- `screening`: `nothing`, `:weak`, or `:chugunov`.
- `neutron_captures`: include neutron-induced reactions in the network
  selection (default `true`).
- `partition_functions`: `:auto` (default; uses `data/winvne_v2.0.dat` when
  present), `nothing`, or a `PartitionFunctionTable` — applied to generated
  detailed-balance reverse rates.
- `jacobian`: `:analytic` (default) or `:finite_difference`.
- `method`: `:backward_euler` (default), `:euler`, `:rk4`, or `:fbdf` for the
  high-order stiff integrator (requires `using OrdinaryDiffEqBDF`).
- `dt_initial`, `dt_min`, `dt_max`: solver step controls; by default they are
  chosen from the trajectory duration like the nova example driver.
- `max_fractional_change`, `max_absolute_change`, `abundance_floor`,
  `max_newton_iterations`, `max_steps`: adaptive controller settings,
  defaulting to the values validated by the nova example driver.
- `output_dir`: when given, write four CSVs there and create the directory
  if needed: `mass_fractions.csv` (mass fraction of every network species at
  every saved trajectory state), `reaction_fluxes.csv` (the flux of every
  reaction at every saved state), `integrated_fluxes.csv` (one row per
  reaction: total flux integrated over the whole run, for an at-a-glance
  "what dominated" summary), and `network.csv` (the static reaction list:
  reactants, products, chapter, source, Q-value). This is the "just run it"
  path — read two input files, pick a solve, say where to write it.
- `write_reaction_fluxes`: skip the flux-history and integrated-flux CSVs
  even when `output_dir` is given (default `true`; the flux pass is one
  extra RHS evaluation per saved state, cheap next to the solve itself).
- `rate_multipliers`, `rate_p_values`: low-level per-reaction rate
  adjustment vectors ordered like the built network's reactions (see
  `network_rhs`); usually easier to reach via `rate_factors`/
  `rate_sample_labels` below, which resolve reaction labels for you once the
  network exists.
- `rate_factors`: a `Dict` (or iterable of pairs) mapping reaction labels to
  a deterministic multiplicative factor, e.g.
  `Dict("22Na(p,γ)23Mg" => 2.0)` doubles that one reaction's rate for the
  whole run while everything else stays nominal (see
  `rate_multipliers_from_factors`). This is the Iliadis-2002-style
  single-rate sensitivity mechanism (Table 8 varies individual rates by
  2, 0.5, 10, 0.1, ...). Mutually exclusive with `rate_multipliers`.
- `rate_sample_labels`: a collection of reaction labels to sample once each
  from their STARLIB lognormal factor uncertainty (`p ~ Normal(0,1)`,
  see `sample_rate_p_values`), holding every other reaction at its nominal
  rate. This is the setup for a future per-reaction Monte Carlo
  sensitivity-run script: call `run_ppn` repeatedly (with a fresh `rng`
  draw each time) to build up a distribution for one named reaction's
  effect. Mutually exclusive with `rate_p_values`.
- `rng`: the random source for `rate_sample_labels` (default
  `Random.default_rng()`); pass your own `MersenneTwister(seed)` for
  reproducible sampling.

# Returns
A named tuple with the `network`, `trajectory`, solution `times` and
abundance `history`, `initial_mass_fractions` and `final_mass_fractions`
dictionaries, `inert_mass_fractions` for species outside the network,
`solver_stats`, `reverse_summary`, the `rate_policy_report` (for
`rates=:iliadis2002`), the network `validation` report, the actual
`rate_multipliers`/`rate_p_values` vectors used (resolved from
`rate_factors`/`rate_sample_labels` if those were given, else whatever was
passed in directly, else `nothing`), and (when `output_dir` is given)
`output_files` naming the CSVs written.
"""
function run_ppn(
    trajectory_path::AbstractString,
    abundance_path::AbstractString;
    rates::Symbol=:starlib,
    tables=nothing,
    screening=:weak,
    neutron_captures::Bool=true,
    generate_detailed_balance::Bool=true,
    partition_functions=:auto,
    jacobian::Symbol=:analytic,
    method::Symbol=:backward_euler,
    max_fractional_change::Real=0.50,
    max_absolute_change::Real=1.0e-4,
    abundance_floor::Real=1.0e-8,
    max_newton_iterations::Integer=80,
    dt_initial=nothing,
    dt_min=nothing,
    dt_max=nothing,
    max_steps::Integer=1_000_000,
    output_dir=nothing,
    write_reaction_fluxes::Bool=true,
    rate_multipliers=nothing,
    rate_p_values=nothing,
    rate_factors=nothing,
    rate_sample_labels=nothing,
    rng::AbstractRNG=Random.default_rng(),
)
    trajectory = read_trajectory(trajectory_path)
    profiles = trajectory_profiles(trajectory)
    X_raw = read_initial_abundances(abundance_path)
    X_normalized = read_initial_abundances(abundance_path; normalize=true)

    rate_policy_report = nothing
    if tables === nothing
        if rates == :iliadis2002
            selection = iliadis2002_rate_tables(; include_reverse=true)
            tables = selection.tables
            rate_policy_report = selection.report
        elseif rates == :starlib
            tables = read_starlib()
        else
            throw(ArgumentError("unsupported rates=$rates; use :starlib or :iliadis2002"))
        end
    end

    projectiles = neutron_captures ? ("p", "he4", "he3", "d", "n") : ("p", "he4", "he3", "d")
    forward_tables = select_h_ca_reaction_tables(tables, keys(X_raw); projectiles=projectiles)

    pf = partition_functions
    if pf === :auto
        pf = isfile(DEFAULT_WINVNE_PATH) ? read_winvne() : nothing
    end
    reverse_summary = add_reverse_reaction_tables(
        tables,
        forward_tables;
        generate_detailed_balance=generate_detailed_balance,
        partition_functions=pf,
    )
    network = network_from_tables(reverse_summary.tables)
    validation = network_validation_report(network; throw_on_error=true)

    if rate_factors !== nothing
        rate_multipliers === nothing || throw(ArgumentError("cannot supply both rate_multipliers and rate_factors"))
        rate_multipliers = rate_multipliers_from_factors(network, rate_factors)
    end
    if rate_sample_labels !== nothing
        rate_p_values === nothing || throw(ArgumentError("cannot supply both rate_p_values and rate_sample_labels"))
        rate_p_values = sample_rate_p_values(network, rate_sample_labels; rng=rng)
    end

    X0 = Dict(name => value for (name, value) in X_normalized if haskey(network.species_index, name))
    inert_mass_fractions = Dict(name => value for (name, value) in X_normalized if !haskey(network.species_index, name))
    Y0 = abundances_from_mass_fractions(network, X0)

    t_start = first(trajectory.time)
    t_end = last(trajectory.time)
    duration = t_end - t_start
    step_initial = dt_initial === nothing ? (duration > 100.0 ? 1.0 : 0.02) : Float64(dt_initial)
    step_min = dt_min === nothing ? (duration > 100.0 ? 1.0e-8 : 1.0e-10) : Float64(dt_min)
    step_max = dt_max === nothing ? (duration > 100.0 ? 20.0 : 0.05) : Float64(dt_max)

    times, history, solver_stats = if method == :fbdf
        solve_network_fbdf(
            network,
            Y0,
            (t_start, t_end),
            profiles.rho,
            profiles.T9;
            screening=screening,
            rate_multipliers=rate_multipliers,
            rate_p_values=rate_p_values,
        )
    else
        solve_network_adaptive(
            network,
            Y0,
            (t_start, t_end),
            step_initial,
            profiles.rho,
            profiles.T9;
            method=method,
            screening=screening,
            jacobian=jacobian,
            max_fractional_change=max_fractional_change,
            max_absolute_change=max_absolute_change,
            abundance_floor=abundance_floor,
            max_newton_iterations=max_newton_iterations,
            dt_min=step_min,
            dt_max=step_max,
            max_steps=max_steps,
            return_stats=true,
            rate_multipliers=rate_multipliers,
            rate_p_values=rate_p_values,
        )
    end

    output_files = nothing
    if output_dir !== nothing
        mkpath(output_dir)
        mass_fraction_path = write_mass_fraction_csv(
            joinpath(output_dir, "mass_fractions.csv"), network, times, history; profiles=profiles,
        )
        network_path = write_network_csv(joinpath(output_dir, "network.csv"), network)
        flux_path = nothing
        integrated_flux_path = nothing
        if write_reaction_fluxes
            flux_history = reaction_flux_history(
                network, history, times, profiles.rho, profiles.T9;
                rate_multipliers=rate_multipliers, rate_p_values=rate_p_values, screening=screening,
            )
            flux_path = write_reaction_flux_csv(joinpath(output_dir, "reaction_fluxes.csv"), network, times, flux_history)
            integrated_flux_path = write_integrated_flux_csv(
                joinpath(output_dir, "integrated_fluxes.csv"), network, integrated_fluxes(times, flux_history),
            )
        end
        output_files = (
            mass_fractions=mass_fraction_path, network=network_path,
            reaction_fluxes=flux_path, integrated_fluxes=integrated_flux_path,
        )
    end

    return (
        network=network,
        trajectory=trajectory,
        times=times,
        history=history,
        initial_mass_fractions=mass_fractions_from_abundances(network, view(history, 1, :)),
        final_mass_fractions=mass_fractions_from_abundances(network, view(history, size(history, 1), :)),
        inert_mass_fractions=inert_mass_fractions,
        mass_fraction_drift=mass_fraction_drift(network, history),
        solver_stats=solver_stats,
        reverse_summary=(explicit=reverse_summary.explicit, generated=reverse_summary.generated, missing=reverse_summary.missing),
        rate_policy_report=rate_policy_report,
        validation=validation,
        output_files=output_files,
        rate_multipliers=rate_multipliers,
        rate_p_values=rate_p_values,
    )
end

# Precompile the solver hot path on a miniature network so first use in a
# fresh session skips most of the compilation latency.

