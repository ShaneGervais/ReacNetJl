using ReacNetJl
using Test

@testset "species names" begin
    @test normalize_species_name("18F") == "f18"
    @test normalize_species_name("4He") == "he4"
    @test normalize_species_name("α") == "he4"
    @test normalize_species_name("γ") == "gamma"
    @test normalize_species_name("p") == "p"
    @test normalize_species_name("26Al*") == "al*6"
    @test normalize_species_name("26Alm") == "al*6"
    @test normalize_species_name("26Alg") == "al26"

    @test species_from_name("18F") == Species("f18", 9, 18)
    @test species_from_name("4He") == Species("he4", 2, 4)
    @test species_from_name("p") == Species("p", 1, 1)
    @test species_from_name("26Al*") == Species("al*6", 13, 26)
    @test species_from_name("n") == Species("n", 0, 1)
    @test species_from_name("d") == Species("d", 1, 2)
    @test species_from_name("t") == Species("t", 1, 3)
    @test_throws ArgumentError species_from_name("notaspecies")
end

@testset "reaction labels" begin
    reactants, products = parse_reaction_label("18F(p,α)15O")
    @test reactants == ["p", "f18"]
    @test products == ["he4", "o15"]

    reactants, products = parse_reaction_label("18F(p,γ)19Ne")
    @test reactants == ["p", "f18"]
    @test products == ["ne19"]

    reactants, products = parse_reaction_label("15O(β+)15N")
    @test reactants == ["o15"]
    @test products == ["n15"]

    reactants, products = parse_reaction_label("13N(beta+)13C")
    @test reactants == ["n13"]
    @test products == ["c13"]

    reactants, products = parse_reaction_label("3He(3He,2p)4He")
    @test reactants == ["he3", "he3"]
    @test products == ["p", "p", "he4"]

    reactants, products = parse_reaction_label("8B(β+)2α")
    @test reactants == ["b8"]
    @test products == ["he4", "he4"]

    reactants, products = parse_reaction_label("p(p,eν)d")
    @test reactants == ["p", "p"]
    @test products == ["d"]

    @test_throws ArgumentError parse_reaction_label("not a reaction")
end

@testset "trajectory input" begin
    trajectory = read_trajectory("test/data/fake_nova_trajectory.dat")
    @test trajectory.time[1] == 0.0
    @test trajectory.time[end] == 3.0
    @test trajectory.T9[5] == 0.40
    @test trajectory.rho[5] == 1.0e4

    profiles = trajectory_profiles(trajectory)
    @test profiles.T9(0.0) == 0.08
    @test profiles.rho(1.0) == 1.0e4
    @test profiles.T9(0.35) ≈ 0.16
    @test profiles.rho(0.35) ≈ 2750.0
    @test first_cooling_threshold_time(trajectory, 0.20) ≈ 1.95
    @test first_cooling_threshold_time(trajectory, 0.05) === nothing
    @test_throws ArgumentError first_cooling_threshold_time(trajectory, 0.0)
    @test_throws ArgumentError profiles.T9(-0.1)

    metadata_path = tempname()
    open(metadata_path, "w") do io
        println(io, "# time T rho")
        println(io, "AGEUNIT = YRS")
        println(io, "TUNIT = T8K")
        println(io, "RHOUNIT = LOG")
        println(io, "0.0 1.0 3.0")
        println(io, "1.0 2.0 4.0")
    end
    metadata_trajectory = read_trajectory(metadata_path)
    rm(metadata_path; force=true)
    @test metadata_trajectory.time[1] == 0.0
    @test metadata_trajectory.time[2] ≈ 365.25 * 24.0 * 60.0 * 60.0
    @test metadata_trajectory.T9 == [0.1, 0.2]
    @test metadata_trajectory.rho == [1.0e3, 1.0e4]

    abundance_path = tempname()
    open(abundance_path, "w") do io
        println(io, "1 PROT 7.0e-1")
        println(io, "2 he 4 2.8e-1")
        println(io, "26 fe56 2.0e-2")
    end
    X = read_initial_abundances(abundance_path)
    X_norm = read_initial_abundances(abundance_path; normalize=true)
    rm(abundance_path; force=true)
    @test X["p"] == 0.7
    @test X["he4"] == 0.28
    @test X["fe56"] == 0.02
    @test sum(values(X_norm)) ≈ 1.0

    bad_path = tempname()
    open(bad_path, "w") do io
        println(io, "0.0 0.1 1000")
        println(io, "0.0 0.2 2000")
    end
    @test_throws ArgumentError read_trajectory(bad_path)
    rm(bad_path; force=true)
end

@testset "abundance conversion" begin
    @test abundance_from_mass_fraction(0.75, 1) == 0.75
    @test abundance_from_mass_fraction(0.25, 4) == 0.0625
    @test mass_fraction_from_abundance(0.0625, 4) == 0.25
end

@testset "rate interpolation" begin
    table = ReactionRateTable(
        2,
        ["p", "f18"],
        ["ne19"],
        "test",
        1.0,
        [0.1, 1.0, 10.0],
        [1.0e-10, 1.0e-5, 1.0],
        [2.0, 2.0, 2.0],
    )

    @test interpolate_rate(table, 1.0) == 1.0e-5
    @test interpolate_rate(table, 0.1) == 1.0e-10
    @test interpolate_rate(table, 10.0) == 1.0
    @test interpolate_rate(table, sqrt(0.1)) ≈ sqrt(1.0e-10 * 1.0e-5)
    @test interpolate_factor_uncertainty(table, 1.0) == 2.0
    @test sampled_interpolate_rate(table, 1.0, 0.0) == interpolate_rate(table, 1.0)
    @test sampled_interpolate_rate(table, 1.0, 1.0) == 2.0 * interpolate_rate(table, 1.0)
    @test sampled_interpolate_rate(table, 1.0, -1.0) == 0.5 * interpolate_rate(table, 1.0)
    @test_throws ArgumentError interpolate_rate(table, 0.01)
end

@testset "reaction flux and network RHS" begin
    species_index = Dict("p" => 1, "f18" => 2, "he4" => 3, "o15" => 4, "c12" => 5, "mg24" => 6)
    Y = zeros(6)
    Y[species_index["p"]] = 0.7
    Y[species_index["f18"]] = 1.0e-5
    Y[species_index["c12"]] = 1.0e-2

    two_body_table = ReactionRateTable(
        4,
        ["p", "f18"],
        ["he4", "o15"],
        "test",
        0.0,
        [0.1, 1.0],
        [2.0, 2.0],
        [2.0, 2.0],
    )
    two_body = Reaction(two_body_table)

    flux = reaction_flux(two_body, Y, species_index, 1000.0, 0.5)
    @test flux ≈ 1000.0 * 2.0 * 0.7 * 1.0e-5
    @test reaction_flux(two_body, Y, species_index, 1000.0, 0.5; rate_p_value=1.0) ≈ 2.0 * flux

    dYdt = network_rhs(Y, [two_body], species_index, 1000.0, 0.5)
    @test dYdt[species_index["p"]] ≈ -flux
    @test dYdt[species_index["f18"]] ≈ -flux
    @test dYdt[species_index["he4"]] ≈ flux
    @test dYdt[species_index["o15"]] ≈ flux

    boosted = network_rhs(Y, [two_body], species_index, 1000.0, 0.5; rate_multipliers=[3.0])
    @test boosted[species_index["p"]] ≈ -3.0 * flux

    one_body_table = ReactionRateTable(
        1,
        ["o15"],
        ["f18"],
        "test",
        0.0,
        [0.1, 1.0],
        [0.25, 0.25],
        [1.0, 1.0],
    )
    one_body = Reaction(one_body_table)
    Y[species_index["o15"]] = 4.0e-3
    @test reaction_flux(one_body, Y, species_index, 1000.0, 0.5) ≈ 0.25 * 4.0e-3

    identical_table = ReactionRateTable(
        4,
        ["c12", "c12"],
        ["mg24"],
        "test",
        0.0,
        [0.1, 1.0],
        [6.0, 6.0],
        [1.0, 1.0],
    )
    identical = Reaction(identical_table)
    identical_flux = reaction_flux(identical, Y, species_index, 10.0, 0.5)
    @test identical_flux ≈ 10.0 * 6.0 * (1.0e-2)^2 / 2.0

    identical_rhs = network_rhs(Y, [identical], species_index, 10.0, 0.5)
    @test identical_rhs[species_index["c12"]] ≈ -2.0 * identical_flux
    @test identical_rhs[species_index["mg24"]] ≈ identical_flux

    identical_network = ReactionNetwork(["c12", "mg24"], [identical])
    identical_balance = species_flux_balance(identical_network, [1.0e-2, 0.0], 10.0, 0.5)
    @test identical_balance.destruction[identical_network.species_index["c12"]] ≈ 2.0 * identical_flux
    @test identical_balance.production[identical_network.species_index["mg24"]] ≈ identical_flux

    p_boosted = network_rhs(Y, [two_body], species_index, 1000.0, 0.5; rate_p_values=[1.0])
    @test p_boosted[species_index["p"]] ≈ -2.0 * flux

    @test_throws ArgumentError reaction_flux(two_body, Y, Dict("p" => 1), 1000.0, 0.5)
    @test_throws ArgumentError network_rhs(Y, [two_body], species_index, 1000.0, 0.5; rate_multipliers=[1.0, 2.0])
    @test_throws ArgumentError network_rhs(Y, [two_body], species_index, 1000.0, 0.5; rate_p_values=[1.0, 2.0])
end

@testset "network construction and time evolution" begin
    table = ReactionRateTable(
        4,
        ["p", "f18"],
        ["he4", "o15"],
        "test",
        0.0,
        [0.1, 1.0],
        [2.0, 2.0],
        [2.0, 2.0],
    )
    reaction = Reaction(table)
    network = ReactionNetwork(["p", "f18", "he4", "o15"], [reaction])

    @test network.species == ["p", "f18", "he4", "o15"]
    @test network.species_info == [Species("p", 1, 1), Species("f18", 9, 18), Species("he4", 2, 4), Species("o15", 8, 15)]
    @test network.species_index["f18"] == 2
    @test length(network.reactions) == 1
    @test length(network.compiled_reactions) == 1
    @test only(network.compiled_reactions).reactant_indices == [1, 2]
    @test only(network.compiled_reactions).product_indices == [3, 4]

    Y_from_X = abundances_from_mass_fractions(
        network,
        Dict("p" => 0.7, "18F" => 1.0e-5, "he4" => 0.28, "o15" => 0.0),
    )
    @test Y_from_X[network.species_index["p"]] == 0.7
    @test Y_from_X[network.species_index["f18"]] ≈ 1.0e-5 / 18
    @test Y_from_X[network.species_index["he4"]] == 0.28 / 4

    X_from_Y = mass_fractions_from_abundances(network, Y_from_X)
    @test X_from_Y["p"] == 0.7
    @test X_from_Y["f18"] ≈ 1.0e-5
    @test X_from_Y["he4"] == 0.28
    @test total_mass_fraction(network, Y_from_X) ≈ 0.98001

    test_history = [Y_from_X'; (2.0 .* Y_from_X)']
    X_history = mass_fraction_history(network, test_history)
    @test size(X_history) == (2, length(network.species))
    @test X_history[1, network.species_index["f18"]] ≈ 1.0e-5

    totals = total_mass_fraction_history(network, test_history)
    @test totals ≈ [0.98001, 1.96002]
    drift = mass_fraction_drift(network, test_history)
    @test drift.initial ≈ 0.98001
    @test drift.final ≈ 1.96002
    @test drift.drift ≈ 0.98001
    @test drift.max_abs_drift ≈ 0.98001
    positivity = abundance_diagnostics(network, test_history)
    @test positivity.min_abundance == 0.0
    @test !positivity.has_negative_abundance
    @test positivity.min_mass_fraction == 0.0

    @test sum(abundances_from_mass_fractions(network, Dict("p" => 7.0, "he4" => 3.0); normalize=true) .* [1, 18, 4, 15]) ≈ 1.0
    @test_throws ArgumentError abundances_from_mass_fractions(network, Dict("ne19" => 1.0))
    @test_throws ArgumentError abundances_from_mass_fractions(network, Dict("p" => 0.5); check_sum=true)
    @test_throws ArgumentError mass_fractions_from_abundances(network, [0.7])
    @test_throws ArgumentError total_mass_fraction(network, [0.7])

    Y0 = [0.7, 1.0e-5, 0.0, 0.0]
    @test reaction_string(reaction) == "18F(p,α)15O"

    fluxes = reaction_fluxes(network, Y0, 1000.0, 0.5)
    @test length(fluxes) == 1
    @test fluxes[1] ≈ reaction_flux(reaction, Y0, network.species_index, 1000.0, 0.5)

    boosted_fluxes = reaction_fluxes(network, Y0, 1000.0, 0.5; rate_multipliers=[2.0])
    @test boosted_fluxes[1] ≈ 2.0 * fluxes[1]
    p_fluxes = reaction_fluxes(network, Y0, 1000.0, 0.5; rate_p_values=[1.0])
    @test p_fluxes[1] ≈ 2.0 * fluxes[1]
    @test_throws ArgumentError reaction_fluxes(network, Y0, 1000.0, 0.5; rate_multipliers=[1.0, 2.0])
    @test_throws ArgumentError reaction_fluxes(network, Y0, 1000.0, 0.5; rate_p_values=[1.0, 2.0])

    rhs = network_rhs(Y0, network, 1000.0, 0.5)
    @test rhs[1] < 0.0
    @test rhs[2] < 0.0
    @test rhs[3] > 0.0
    @test rhs[4] > 0.0

    balance = species_flux_balance(network, Y0, 1000.0, 0.5)
    @test balance.net ≈ rhs
    @test balance.destruction[network.species_index["p"]] ≈ fluxes[1]
    @test balance.production[network.species_index["o15"]] ≈ fluxes[1]

    screening_factor = weak_screening_multiplier(network, reaction, Y0, 1000.0, 0.5)
    @test screening_factor > 1.0
    screened_fluxes = reaction_fluxes(network, Y0, 1000.0, 0.5; screening=:weak)
    @test screened_fluxes[1] ≈ screening_factor * fluxes[1]
    screened_rhs = network_rhs(Y0, network, 1000.0, 0.5; screening=:weak)
    @test screened_rhs[network.species_index["f18"]] ≈ screening_factor * rhs[network.species_index["f18"]]

    edges = reaction_edges(network)
    @test length(edges) == 4
    @test (reaction_index=1, reaction="18F(p,α)15O", from="f18", to="o15") in edges
    @test (reaction_index=1, reaction="18F(p,α)15O", from="p", to="he4") in edges

    conservation = reaction_conservation(reaction)
    @test conservation.conserves_A
    @test conservation.conserves_Z
    @test conservation.reactant_A == 19
    @test conservation.product_A == 19
    @test conservation.reactant_Z == 10
    @test conservation.product_Z == 10

    decay_table = ReactionRateTable(
        1,
        ["o15"],
        ["n15"],
        "test",
        0.0,
        [0.1, 1.0],
        [0.1, 0.1],
        [1.0, 1.0],
    )
    decay = Reaction(decay_table)
    @test reaction_string(decay) == "15O(β+)15N"
    decay_conservation = reaction_conservation(decay)
    @test decay_conservation.conserves_A
    @test !decay_conservation.conserves_Z
    @test decay_conservation.is_weak_decay
    @test decay_conservation.valid_nuclear_bookkeeping
    decay_network = ReactionNetwork(["o15", "n15"], [decay])
    @test network_validation_report(decay_network).valid

    b8_decay_table = ReactionRateTable(
        2,
        ["b8"],
        ["he4", "he4"],
        "testw",
        0.0,
        [0.1, 1.0],
        [0.1, 0.1],
        [1.0, 1.0],
    )
    b8_decay = Reaction(b8_decay_table)
    @test reaction_string(b8_decay) == "8B(β+)2α"
    @test reaction_conservation(b8_decay).valid_nuclear_bookkeeping

    pp_table = ReactionRateTable(
        4,
        ["p", "p"],
        ["d"],
        "nacr",
        1.442,
        [0.1, 1.0],
        [1.0e-20, 1.0e-20],
        [1.0, 1.0],
    )
    pp_reaction = Reaction(pp_table)
    @test reaction_string(pp_reaction) == "p(p,eν)d"
    @test reaction_conservation(pp_reaction).valid_nuclear_bookkeeping

    validation = network_validation_report(network)
    @test validation.valid
    @test validation.num_species == 4
    @test validation.num_reactions == 1
    @test isempty(validation.issues)
    @test only(validation.reaction_reports).valid

    invalid_table = ReactionRateTable(
        4,
        ["p", "f18"],
        ["o15"],
        "test",
        0.0,
        [0.1, 1.0],
        [2.0, 2.0],
        [1.0, 1.0],
    )
    invalid_network = ReactionNetwork(["p", "f18", "o15"], [Reaction(invalid_table)])
    invalid_conservation = reaction_conservation(only(invalid_network.reactions))
    @test !invalid_conservation.conserves_A
    @test !invalid_conservation.conserves_Z
    invalid_report = network_validation_report(invalid_network)
    @test !invalid_report.valid
    @test !only(invalid_report.reaction_reports).valid
    @test any(contains("baryon number"), invalid_report.issues)
    @test any(contains("charge"), invalid_report.issues)
    @test_throws ArgumentError network_validation_report(invalid_network; throw_on_error=true)

    times, history = solve_network(network, Y0, (0.0, 1.0e-3), 1.0e-4, 1000.0, 0.5; method=:euler)
    @test times[1] == 0.0
    @test times[end] == 1.0e-3
    @test size(history) == (length(times), length(Y0))
    @test history[end, network.species_index["f18"]] < Y0[network.species_index["f18"]]
    @test history[end, network.species_index["o15"]] > Y0[network.species_index["o15"]]
    f18_initial_X = mass_fraction_from_abundance(Y0[network.species_index["f18"]], 18)
    f18_final_X = mass_fraction_from_abundance(history[end, network.species_index["f18"]], 18)
    f18_threshold = 0.5 * (f18_initial_X + f18_final_X)
    crossing = first_mass_fraction_threshold_crossing(network, times, history, "f18", f18_threshold)
    @test crossing !== nothing
    @test first(times) < crossing.time < last(times)
    @test mass_fraction_from_abundance(crossing.state[network.species_index["f18"]], 18) ≈ f18_threshold
    @test first_mass_fraction_threshold_crossing(network, times, history, "f18", f18_final_X / 2) === nothing
    @test_throws ArgumentError first_mass_fraction_threshold_crossing(network, times, history, "missing", f18_threshold)

    flux_history = reaction_flux_history(network, history, times, 1000.0, 0.5)
    @test size(flux_history) == (length(times), length(network.reactions))
    @test flux_history[1, 1] ≈ fluxes[1]
    @test integrated_fluxes([0.0, 2.0], [1.0 2.0; 3.0 4.0]) ≈ [4.0, 6.0]
    @test length(integrated_fluxes(times, flux_history)) == length(network.reactions)
    @test integrated_fluxes(times, flux_history)[1] > 0.0
    @test_throws ArgumentError reaction_flux_history(network, history, times[1:end-1], 1000.0, 0.5)
    @test_throws ArgumentError integrated_fluxes([0.0], [1.0 2.0])
    @test_throws ArgumentError integrated_fluxes([0.0, -1.0], reshape([1.0, 2.0], 2, 1))

    epsilon = energy_generation_rate(network, Y0, 1000.0, 0.5)
    @test epsilon ≈ reaction.rate_table.q_value * fluxes[1] * ReacNetJl.AVOGADRO * ReacNetJl.MEV_TO_ERG
    epsilon_history = energy_generation_history(network, history, times, 1000.0, 0.5)
    @test length(epsilon_history) == length(times)
    @test epsilon_history[1] ≈ energy_generation_rate(network, history[1, :], 1000.0, 0.5)
    @test integrated_energy_generation([0.0, 2.0], [1.0, 3.0]) ≈ 4.0
    @test_throws ArgumentError energy_generation_history(network, history, times[1:end-1], 1000.0, 0.5)
    @test_throws ArgumentError integrated_energy_generation([0.0], [1.0])

    _, boosted_history = solve_network(network, Y0, (0.0, 1.0e-3), 1.0e-4, 1000.0, 0.5; method=:rk4, rate_multipliers=[2.0])
    @test boosted_history[end, network.species_index["f18"]] < history[end, network.species_index["f18"]]

    _, p_history = solve_network(network, Y0, (0.0, 1.0e-3), 1.0e-4, 1000.0, 0.5; method=:rk4, rate_p_values=[1.0])
    @test p_history[end, network.species_index["f18"]] < history[end, network.species_index["f18"]]

    be_times, be_history = solve_network(network, Y0, (0.0, 1.0e-3), 1.0e-4, 1000.0, 0.5; method=:backward_euler)
    @test be_times[end] == 1.0e-3
    @test be_history[end, network.species_index["f18"]] < Y0[network.species_index["f18"]]
    @test be_history[end, network.species_index["o15"]] > Y0[network.species_index["o15"]]

    be_stats_times, be_stats_history, be_stats = solve_network(network, Y0, (0.0, 1.0e-3), 1.0e-4, 1000.0, 0.5; method=:backward_euler, return_stats=true)
    @test be_stats_times == be_times
    @test be_stats_history ≈ be_history
    @test be_stats.accepted_steps == length(be_stats_times) - 1
    @test length(be_stats.newton_iterations) == be_stats.accepted_steps
    @test be_stats.max_newton_iterations >= 0
    @test be_stats.newton_failed_steps == 0

    mc = run_monte_carlo(network, Y0, (0.0, 1.0e-4), 1.0e-4, 1000.0, 0.5; nruns=3, seed=1234, store_histories=true)
    mc_repeat = run_monte_carlo(network, Y0, (0.0, 1.0e-4), 1.0e-4, 1000.0, 0.5; nruns=3, seed=1234)
    @test size(mc.final_abundances) == (3, length(network.species))
    @test size(mc.rate_p_values) == (3, length(network.reactions))
    @test length(mc.histories) == 3
    @test mc.species == network.species
    @test mc.reactions == [reaction_string(reaction)]
    @test mc.rate_p_values == mc_repeat.rate_p_values
    @test mc.final_abundances == mc_repeat.final_abundances
    @test_throws ArgumentError run_monte_carlo(network, Y0, (0.0, 1.0e-4), 1.0e-4, 1000.0, 0.5; nruns=0)

    rho_profile(t) = 1000.0 + t
    T9_profile(t) = 0.5
    profile_times, profile_history = solve_network(network, Y0, (0.0, 2.5e-4), 1.0e-4, rho_profile, T9_profile)
    @test profile_times[end] == 2.5e-4
    @test size(profile_history, 1) == length(profile_times)

    adaptive_times, adaptive_history = solve_network_adaptive(
        network,
        Y0,
        (0.0, 1.0e-3),
        5.0e-4,
        1000.0,
        0.5;
        max_fractional_change=0.05,
        dt_min=1.0e-8,
        dt_max=5.0e-4,
    )
    @test adaptive_times[1] == 0.0
    @test adaptive_times[end] ≈ 1.0e-3
    @test size(adaptive_history, 1) == length(adaptive_times)
    @test adaptive_history[end, network.species_index["f18"]] < Y0[network.species_index["f18"]]

    stats_times, stats_history, stats = solve_network_adaptive(
        network,
        Y0,
        (0.0, 1.0e-3),
        5.0e-4,
        1000.0,
        0.5;
        max_fractional_change=0.05,
        max_absolute_change=1.0e-8,
        dt_min=1.0e-8,
        dt_max=5.0e-4,
        return_stats=true,
    )
    @test stats_times[end] ≈ 1.0e-3
    @test size(stats_history, 1) == length(stats_times)
    @test stats.accepted_steps == length(stats_times) - 1
    @test stats.rejected_steps >= 0
    @test stats.min_dt <= stats.max_dt
    @test stats.max_absolute_change >= 0.0

    be_adaptive_times, _, be_adaptive_stats = solve_network_adaptive(
        network,
        Y0,
        (0.0, 1.0e-3),
        5.0e-4,
        1000.0,
        0.5;
        method=:backward_euler,
        max_fractional_change=0.05,
        max_absolute_change=1.0e-8,
        dt_min=1.0e-8,
        dt_max=5.0e-4,
        return_stats=true,
    )
    @test be_adaptive_times[end] ≈ 1.0e-3
    @test length(be_adaptive_stats.newton_iterations) == be_adaptive_stats.accepted_steps
    @test be_adaptive_stats.newton_failed_steps >= 0

    trajectory = read_trajectory("test/data/fake_nova_trajectory.dat")
    profiles = trajectory_profiles(trajectory)
    traj_times, traj_history = solve_network_adaptive(network, Y0, (0.2, 2.2), 0.5, profiles.rho, profiles.T9; max_fractional_change=0.5, dt_min=1.0e-6, dt_max=0.5)
    @test traj_times[end] ≈ 2.2
    @test size(traj_history, 2) == length(network.species)

    @test_throws ArgumentError ReactionNetwork(["p", "p"], [reaction])
    @test_throws ArgumentError ReactionNetwork(["p", "f18", "he4"], [reaction])
    @test_throws ArgumentError solve_network(network, [0.7], (0.0, 1.0), 0.1, 1000.0, 0.5)
    @test_throws ArgumentError solve_network(network, Y0, (1.0, 0.0), 0.1, 1000.0, 0.5)
    @test_throws ArgumentError solve_network(network, Y0, (0.0, 1.0), -0.1, 1000.0, 0.5)
    @test_throws ArgumentError solve_network(network, Y0, (0.0, 1.0), 0.1, 1000.0, 0.5; method=:unknown)
    @test_throws ArgumentError solve_network(network, Y0, (0.0, 1.0), 0.1, 1000.0, 0.5; method=:backward_euler, max_newton_iterations=0)
    @test_throws ArgumentError solve_network_adaptive(network, Y0, (0.0, 1.0), 0.1, 1000.0, 0.5; max_fractional_change=0.0)
    @test_throws ArgumentError solve_network_adaptive(network, Y0, (0.0, 1.0), 0.1, 1000.0, 0.5; max_absolute_change=0.0)
end

@testset "read synthetic STARLIB table" begin
    path = tempname()
    open(path, "w") do io
        println(io, "4 p f18 ne19 testsource 3.529")
        for i in 1:60
            T9 = 0.01 * i
            rate = 1.0e-12 * i
            factor = 2.0
            println(io, T9, " ", rate, " ", factor)
        end
    end

    tables = read_starlib(path)
    rm(path; force=true)

    @test length(tables) == 1
    table = only(tables)
    @test table.chapter == 4
    @test table.reactants == ["p", "f18"]
    @test table.products == ["ne19"]
    @test table.source == "testsource"
    @test table.q_value == 3.529
    @test length(table.T9) == 60
    @test table.T9[1] == 0.01
    @test table.rate[end] == 6.0e-11
    @test table.factor_uncertainty[end] == 2.0
    report = starlib_chapter_report(tables)
    @test report.total == 1
    @test report.supported == 1
    @test report.unsupported == 0

    matches = find_rate(tables, "18F(p,γ)19Ne")
    @test matches == tables
    @test find_rate(tables, "18F(p,γ)19Ne"; source="missing") == ReactionRateTable[]
    @test isempty(find_reverse_rate(tables, "18F(p,γ)19Ne"))

    reaction = reaction_from_label(tables, "18F(p,γ)19Ne")
    @test reaction.reactants == ["p", "f18"]
    @test reaction.products == ["ne19"]

    network = network_from_labels(tables, ["18F(p,γ)19Ne"])
    @test network.species == ["p", "f18", "ne19"]
    @test network.species_info == [Species("p", 1, 1), Species("f18", 9, 18), Species("ne19", 10, 19)]

    network_with_species = network_from_labels(tables, ["18F(p,γ)19Ne"]; species=["p", "18F", "19Ne"])
    @test network_with_species.species == ["p", "f18", "ne19"]

    duplicate_tables = [tables; tables]
    @test_throws ArgumentError reaction_from_label(duplicate_tables, "18F(p,γ)19Ne")
    @test reaction_from_label(duplicate_tables, "18F(p,γ)19Ne"; on_multiple=:first).rate_table == table
    @test_throws ArgumentError reaction_from_label(tables, "18F(p,α)15O")

    multiproduct_path = tempname()
    open(multiproduct_path, "w") do io
        println(io, "2 b8 he4 he4 weaksrcw 18.0717")
        for i in 1:60
            println(io, 0.01 * i, " ", 0.1, " ", 1.0)
        end
        println(io, "6 he3 he3 p p he4 il16 12.8596")
        for i in 1:60
            println(io, 0.01 * i, " ", 0.2, " ", 1.0)
        end
        println(io, "8 he4 he4 he4 c12 nacr 7.27475")
        for i in 1:60
            println(io, 0.01 * i, " ", 0.3, " ", 1.0)
        end
        println(io, "4 p p d nacr 1.44222")
        for i in 1:60
            println(io, 0.01 * i, " ", 1.0e-20, " ", 1.0)
        end
    end
    multiproduct_tables = read_starlib(multiproduct_path)
    rm(multiproduct_path; force=true)
    @test multiproduct_tables[1].reactants == ["b8"]
    @test multiproduct_tables[1].products == ["he4", "he4"]
    @test multiproduct_tables[2].reactants == ["he3", "he3"]
    @test multiproduct_tables[2].products == ["p", "p", "he4"]
    @test multiproduct_tables[3].reactants == ["he4", "he4", "he4"]
    @test multiproduct_tables[3].products == ["c12"]
    @test multiproduct_tables[4].reactants == ["p", "p"]
    @test multiproduct_tables[4].products == ["d"]
    @test network_validation_report(network_from_tables(multiproduct_tables); throw_on_error=true).valid
    @test find_rate(multiproduct_tables, "8B(β+)2α") == [multiproduct_tables[1]]
    @test find_rate(multiproduct_tables, "3He(3He,2p)4He") == [multiproduct_tables[2]]
    @test find_rate(multiproduct_tables, "p(p,eν)d") == [multiproduct_tables[4]]
    selected_multiproduct_h_ca = select_h_ca_reaction_tables(multiproduct_tables, ["b8"])
    @test multiproduct_tables[1] in selected_multiproduct_h_ca

    al_path = tempname()
    open(al_path, "w") do io
        println(io, "4 p al*6 si27 isomer 7.69126")
        for i in 1:60
            println(io, 0.01 * i, " ", 1.0e-12 * i, " ", 2.0)
        end
    end
    al_tables = read_starlib(al_path)
    rm(al_path; force=true)
    al_reaction = reaction_from_label(al_tables, "26Al*(p,γ)27Si")
    @test al_reaction.reactants == ["p", "al*6"]
    @test al_reaction.products == ["si27"]

    reverse_table = ReactionRateTable(2, ["ne19"], ["p", "f18"], "reverse", -3.529, table.T9, table.rate, table.factor_uncertainty)
    reverse_matches = find_reverse_rate([table, reverse_table], "18F(p,γ)19Ne")
    @test reverse_matches == [reverse_table]
    @test normalize_species_name("al-6") == "al26"

    generated_reverse = generated_detailed_balance_reverse_table(table)
    @test generated_reverse.reactants == ["ne19"]
    @test generated_reverse.products == ["p", "f18"]
    @test generated_reverse.q_value == -table.q_value
    @test all(>(0.0), generated_reverse.rate)

    reverse_summary = add_reverse_reaction_tables([table, reverse_table], [table])
    @test reverse_summary.explicit == 1
    @test reverse_summary.generated == 0
    @test length(reverse_summary.tables) == 2

    generated_summary = add_reverse_reaction_tables([table], [table])
    @test generated_summary.explicit == 0
    @test generated_summary.generated == 1
    @test length(generated_summary.tables) == 2
    generated_network = network_from_tables(generated_summary.tables)
    @test length(generated_network.reactions) == 2

    decay_result = decay_mass_fractions(multiproduct_tables, Dict("b8" => 8.0e-8), 10.0)
    @test length(decay_result.decay_tables) == 1
    @test get(decay_result.mass_fractions, "b8", 0.0) < 8.0e-8
    @test get(decay_result.mass_fractions, "he4", 0.0) > 0.0

    long_decay_result = decay_mass_fractions(multiproduct_tables, Dict("b8" => 8.0e-8), 3600.0)
    @test long_decay_result.times == [0.0, 3600.0]
    @test get(long_decay_result.mass_fractions, "b8", 0.0) < 1.0e-30
    @test sum(values(long_decay_result.mass_fractions); init=0.0) ≈ 8.0e-8 rtol=1.0e-10

    selected_h_ca = select_h_ca_reaction_tables(tables, ["p", "f18"])
    @test !isempty(selected_h_ca)
    @test any(t -> t.reactants == ["p", "f18"] && t.products == ["ne19"], selected_h_ca)

    result = solve_single_zone(
        tables,
        ["18F(p,γ)19Ne"],
        Dict("p" => 0.7, "18F" => 1.0e-5, "19Ne" => 0.0),
        (0.0, 1.0e-4),
        1.0e-4,
        1000.0,
        0.5;
        adaptive=false,
    )
    @test result.validation.valid
    @test result.network.species == ["p", "f18", "ne19"]
    @test result.initial_mass_fractions["f18"] ≈ 1.0e-5
    @test result.final_mass_fractions["f18"] < result.initial_mass_fractions["f18"]
    @test size(result.reaction_fluxes) == (length(result.times), length(result.network.reactions))
    @test length(result.integrated_fluxes) == 1
    @test size(result.mass_fraction_history) == size(result.abundances)
    @test result.mass_fraction_drift.initial ≈ total_mass_fraction(result.network, result.abundances[1, :])
    @test !result.abundance_diagnostics.has_negative_abundance
    @test result.solver_stats.accepted_steps == length(result.times) - 1
    @test length(result.energy_generation) == length(result.times)
    @test isfinite(result.integrated_energy_generation)

    unsupported_path = tempname()
    open(unsupported_path, "w") do io
        println(io, "11 p f18 he4 o15 ne19 unsupported 0.0")
        for i in 1:60
            println(io, 0.01 * i, " ", 1.0e-12 * i, " ", 2.0)
        end
    end

    unsupported_tables = read_starlib(unsupported_path)
    rm(unsupported_path; force=true)
    unsupported_report = starlib_chapter_report(unsupported_tables)
    @test unsupported_report.total == 1
    @test unsupported_report.supported == 0
    @test unsupported_report.unsupported == 1
    @test unsupported_report.unsupported_by_chapter[11] == 1
    @test isempty(only(unsupported_tables).products)
end
