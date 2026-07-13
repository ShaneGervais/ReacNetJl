module ReacNetJl

using LinearAlgebra
using Printf
using Random
import Downloads
using PrecompileTools: @setup_workload, @compile_workload

export Species,
    Trajectory,
    ReactionRateTable,
    ReaclibSet,
    Reaction,
    ReactionNetwork,
    species_from_name,
    abundance_from_mass_fraction,
    mass_fraction_from_abundance,
    normalize_species_name,
    parse_reaction_label,
    read_starlib,
    read_reaclib,
    read_winvne,
    read_tabulated_rates,
    read_iliadis2001_rates,
    read_nacre_rates,
    read_ame_masses,
    reaction_q_value,
    fetch_data!,
    PartitionFunctionTable,
    reaclib_rate,
    reaclib_rate_tables,
    iliadis2002_rate_tables,
    read_trajectory,
    read_initial_abundances,
    trajectory_profiles,
    first_cooling_threshold_time,
    first_mass_fraction_threshold_crossing,
    starlib_chapter_report,
    find_rate,
    find_reverse_rate,
    reaction_from_label,
    generated_detailed_balance_reverse_table,
    add_reverse_reaction_tables,
    select_h_ca_reaction_tables,
    select_decay_reaction_tables,
    decay_mass_fractions,
    network_from_tables,
    network_from_labels,
    weak_screening_multiplier,
    chugunov_screening_multiplier,
    abundances_from_mass_fractions,
    mass_fractions_from_abundances,
    mass_fraction_history,
    total_mass_fraction,
    total_mass_fraction_history,
    mass_fraction_drift,
    abundance_diagnostics,
    interpolate_rate,
    interpolate_factor_uncertainty,
    sampled_interpolate_rate,
    reaction_string,
    reaction_flux,
    reaction_fluxes,
    reaction_flux_history,
    integrated_fluxes,
    energy_generation_rate,
    energy_generation_history,
    integrated_energy_generation,
    species_flux_balance,
    reaction_edges,
    reaction_conservation,
    network_validation_report,
    network_rhs,
    solve_network,
    solve_network_adaptive,
    solve_network_fbdf,
    run_monte_carlo,
    solve_single_zone,
    run_ppn,
    write_mass_fraction_csv,
    write_reaction_flux_csv,
    write_network_csv

# `@__DIR__` is resolved per source file, not per module, so it must be
# captured once here (where it correctly means the package's src/ directory)
# rather than recomputed in io/paths.jl, which would instead resolve to
# src/io/.
const _PACKAGE_ROOT = dirname(@__DIR__)

include("physics/species.jl")
include("io/paths.jl")
include("io/starlib.jl")
include("io/rate_tables.jl")
include("io/reaclib.jl")
include("io/winvne.jl")
include("io/ame.jl")
include("io/tabulated_rates.jl")
include("io/iliadis2002_policy.jl")
include("io/trajectory.jl")
include("physics/network.jl")
include("physics/network_selection.jl")
include("physics/abundances.jl")
include("physics/reverse_rates.jl")
include("physics/screening.jl")
include("physics/flux.jl")
include("physics/decay.jl")
include("physics/reaction_string.jl")
include("physics/validation.jl")
include("solver/step_cache.jl")
include("solver/explicit.jl")
include("solver/backward_euler.jl")
include("solver/fixed_step.jl")
include("solver/adaptive.jl")
include("solver/monte_carlo.jl")
include("solver/fbdf.jl")
include("io/output.jl")
include("api.jl")

@setup_workload begin
    _pc_grid = STARLIB_T9_GRID
    _pc_unit = ones(length(_pc_grid))
    _pc_tables = [
        ReactionRateTable(4, ["p", "o17"], ["f18"], "pc", 5.607, _pc_grid, fill(2.0e4, length(_pc_grid)), _pc_unit),
        ReactionRateTable(5, ["p", "f18"], ["he4", "o15"], "pc", 2.882, _pc_grid, fill(6.0e4, length(_pc_grid)), _pc_unit),
        ReactionRateTable(1, ["o15"], ["n15"], "pcw", 2.754, _pc_grid, fill(4.5e-3, length(_pc_grid)), _pc_unit),
    ]

    @compile_workload begin
        _pc_network = network_from_tables(_pc_tables)
        _pc_X0 = Dict("p" => 0.7, "o17" => 0.1, "f18" => 1.0e-4, "he4" => 0.2, "o15" => 0.0, "n15" => 0.0)
        _pc_Y0 = abundances_from_mass_fractions(_pc_network, _pc_X0; normalize=true)
        for _pc_screening in (nothing, :weak, :chugunov)
            network_rhs(_pc_Y0, _pc_network, 500.0, 0.25; screening=_pc_screening)
            solve_network_adaptive(
                _pc_network, _pc_Y0, (0.0, 1.0e-6), 1.0e-7, 500.0, 0.25;
                method=:backward_euler, screening=_pc_screening, return_stats=true,
                abundance_floor=1.0e-8, max_newton_iterations=80,
            )
        end
        solve_network(_pc_network, _pc_Y0, (0.0, 1.0e-7), 1.0e-8, 500.0, 0.25; method=:rk4, screening=:weak)
        mass_fractions_from_abundances(_pc_network, _pc_Y0)
    end
end

end # module
