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

using Printf

@testset "REACLIB parsing, evaluation, and Iliadis-2002 policy" begin
    reaclib_header(species, label, res, rev, q) =
        " "^5 * join([lpad(s, 5) for s in vcat(species, fill("", 6 - length(species)))], "") *
        " "^8 * lpad(label, 4) * res * rev * " "^3 * @sprintf("%12.5e", q)
    reaclib_coefficients(a) =
        join([@sprintf("%13.6e", x) for x in a[1:4]], "") * "\n" *
        join([@sprintf("%13.6e", x) for x in a[5:7]], "")
    reaclib_set(chapter, species, label, res, rev, q, a) =
        string(chapter) * "\n" * reaclib_header(species, label, res, rev, q) * "\n" * reaclib_coefficients(a)

    fixture = join([
        # p + o17 -> he4 + n14: nacr non-resonant + resonant sets, plus an il01 alternative
        reaclib_set(5, ["p", "o17", "he4", "n14"], "nacr", 'n', ' ', 1.1917, (2.0, 0, 0, 0, 0, 0, 0)),
        reaclib_set(5, ["p", "o17", "he4", "n14"], "nacr", 'r', ' ', 1.1917, (1.0, 0, 0, 0, 0, 0, 0)),
        reaclib_set(5, ["p", "o17", "he4", "n14"], "il01", 'n', ' ', 1.1917, (0.5, 0, 0, 0, 0, 0, 0)),
        # p + ne20 -> na21: il01 and nacr alternatives
        reaclib_set(4, ["p", "ne20", "na21"], "il01", 'n', ' ', 2.431, (1.5, 0, 0, 0, 0, 0, 0)),
        reaclib_set(4, ["p", "ne20", "na21"], "nacr", 'n', ' ', 2.431, (1.0, 0, 0, 0, 0, 0, 0)),
        # reverse fits for both labels; only the chosen il01 one may be imported
        reaclib_set(2, ["na21", "p", "ne20"], "il01", 'n', 'v', -2.431, (3.0, 0, 0, 0, 0, 0, 0)),
        reaclib_set(2, ["na21", "p", "ne20"], "nacr", 'n', 'v', -2.431, (2.5, 0, 0, 0, 0, 0, 0)),
        # weak decay and an uncompiled fallback label
        reaclib_set(1, ["f18", "o18"], "wc12", 'w', ' ', 1.6555, (-9.3, 0, 0, 0, 0, 0, 0)),
        reaclib_set(4, ["p", "f17", "ne18"], "dc11", 'n', ' ', 3.9224, (0.25, 0, 0, 0, 0, 0, 0)),
    ], "\n") * "\n"

    fixture_path = joinpath(mktempdir(), "reaclib_fixture.dat")
    write(fixture_path, fixture)

    sets = read_reaclib(fixture_path)
    @test length(sets) == 9
    @test sets[1].chapter == 5
    @test sets[1].reactants == ["p", "o17"]
    @test sets[1].products == ["he4", "n14"]
    @test sets[1].label == "nacr"
    @test sets[1].resonance == 'n'
    @test !sets[1].reverse
    @test sets[1].q_value ≈ 1.1917 atol = 1.0e-4
    @test sets[1].coefficients[1] ≈ 2.0
    @test sets[6].reverse
    @test sets[8].resonance == 'w'

    # analytic evaluation against the REACLIB parameterization
    a = (1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0)
    full_set = ReaclibSet(4, ["p", "ne20"], ["na21"], "test", 'n', false, 2.431, a)
    T9 = 0.7
    expected_exponent = a[1] + a[2] / T9 + a[3] * T9^(-1 / 3) + a[4] * T9^(1 / 3) +
                        a[5] * T9 + a[6] * T9^(5 / 3) + a[7] * log(T9)
    @test reaclib_rate(full_set, T9) ≈ exp(expected_exponent) rtol = 1.0e-12
    @test_throws ArgumentError reaclib_rate(full_set, 0.0)

    # grouped tables: sets of one reaction/label are summed, reverse excluded
    tables = reaclib_rate_tables(sets)
    @test length(tables) == 6
    @test all(t -> length(t.T9) == 60, tables)
    @test all(t -> all(==(1.0), t.factor_uncertainty), tables)

    o17_nacr = only(find_rate(tables, "17O(p,α)14N"; source="nacr"))
    @test o17_nacr.rate[1] ≈ exp(2.0) + exp(1.0) rtol = 1.0e-12
    @test o17_nacr.q_value ≈ 1.1917 atol = 1.0e-4

    f18_decay = only(find_rate(tables, "18F(β+)18O"))
    @test f18_decay.source == "wc12w"
    @test f18_decay.rate[1] ≈ exp(-9.3) rtol = 1.0e-12

    nacr_only = reaclib_rate_tables(sets; labels=["nacr"])
    @test length(nacr_only) == 2
    @test all(t -> t.source == "nacr", nacr_only)

    # Iliadis-2002 policy: nacr below A=20, il01 at and above, fallbacks reported
    result = iliadis2002_rate_tables(sets)
    @test only(find_rate(result.tables, "17O(p,α)14N")).source == "nacr"
    ne20_capture = only(find_rate(result.tables, "20Ne(p,γ)21Na"))
    @test ne20_capture.source == "il01"
    @test ne20_capture.rate[1] ≈ exp(1.5) rtol = 1.0e-12
    @test result.report.counts[:nacr] == 1
    @test result.report.counts[:il01] == 1
    @test result.report.counts[:weak] == 1
    @test result.report.counts[:other] == 1
    @test length(result.report.fallbacks) == 1
    @test occursin("f17", only(result.report.fallbacks).reaction)
    @test only(result.report.fallbacks).preferred == "nacr"

    without_weak = iliadis2002_rate_tables(sets; include_weak=false)
    @test isempty(find_rate(without_weak.tables, "18F(β+)18O"))

    with_reverse = iliadis2002_rate_tables(sets; include_reverse=true)
    na21_reverse = find_reverse_rate(with_reverse.tables, "20Ne(p,γ)21Na")
    @test length(na21_reverse) == 1
    @test only(na21_reverse).source == "il01v"
    @test only(na21_reverse).rate[1] ≈ exp(3.0) rtol = 1.0e-12

    # duplicated input (overlapping label files) must not double any rate
    doubled = iliadis2002_rate_tables(vcat(sets, sets))
    @test only(find_rate(doubled.tables, "17O(p,α)14N")).rate[1] ≈ exp(2.0) + exp(1.0) rtol = 1.0e-12

    # a REACLIB table drives a network exactly like a STARLIB table
    network = network_from_tables(find_rate(result.tables, "20Ne(p,γ)21Na"))
    @test network.species == ["p", "ne20", "na21"]

    bad_path = joinpath(mktempdir(), "bad_chapter.dat")
    write(bad_path, "99\n" * reaclib_header(["p", "o17", "he4", "n14"], "nacr", 'n', ' ', 1.0) * "\n" * reaclib_coefficients((1.0, 0, 0, 0, 0, 0, 0)) * "\n")
    @test_throws ErrorException read_reaclib(bad_path)
end

@testset "cached RHS and analytic Jacobian equivalence" begin
    grid = ReacNetJl.STARLIB_T9_GRID
    unit_uncertainty = ones(length(grid))
    capture = ReactionRateTable(4, ["p", "o17"], ["f18"], "test", 5.607, grid, fill(2.0e4, length(grid)), unit_uncertainty)
    proton_alpha = ReactionRateTable(5, ["p", "f18"], ["he4", "o15"], "test", 2.882, grid, fill(6.0e4, length(grid)), unit_uncertainty)
    decay = ReactionRateTable(1, ["o15"], ["n15"], "testw", 2.754, grid, fill(4.5e-3, length(grid)), unit_uncertainty)
    triple_alpha = ReactionRateTable(8, ["he4", "he4", "he4"], ["c12"], "test", 7.275, grid, fill(1.0e-9, length(grid)), unit_uncertainty)
    network = network_from_tables([capture, proton_alpha, decay, triple_alpha])

    X0 = Dict(
        "p" => 0.6, "o17" => 0.05, "f18" => 1.0e-4, "he4" => 0.3,
        "o15" => 1.0e-5, "n15" => 0.0, "c12" => 0.05,
    )
    Y = abundances_from_mass_fractions(network, X0; normalize=true)
    rho, T9 = 500.0, 0.25
    multipliers = [1.0, 2.5, 1.0, 0.7]

    for screening in (nothing, :weak), rate_multipliers in (nothing, multipliers)
        reference = network_rhs(Y, network, rho, T9; screening=screening, rate_multipliers=rate_multipliers)
        cache = ReacNetJl._build_step_cache(network, rho, T9; screening=screening, rate_multipliers=rate_multipliers)
        fast = ReacNetJl._cached_network_rhs!(similar(Y), network, cache, Y)
        @test isapprox(fast, reference; rtol=1.0e-12, atol=1.0e-300)
    end

    # custom screening functions cannot be cached and must fall back
    @test ReacNetJl._build_step_cache(network, rho, T9; screening=(n, r, y, d, t) -> 2.0) === nothing

    # analytic Jacobian against a central finite difference of the reference RHS
    n = length(Y)
    function fd_jacobian(screening)
        J = Matrix{Float64}(undef, n, n)
        for j in 1:n
            h = 1.0e-6 * max(abs(Y[j]), 1.0e-8)
            Y_plus = copy(Y); Y_plus[j] += h
            Y_minus = copy(Y); Y_minus[j] -= h
            column = (network_rhs(Y_plus, network, rho, T9; screening=screening) .-
                      network_rhs(Y_minus, network, rho, T9; screening=screening)) ./ (2.0 * h)
            J[:, j] .= column
        end
        return J
    end

    # without screening the analytic Jacobian is exact: only FD truncation error
    cache_unscreened = ReacNetJl._build_step_cache(network, rho, T9)
    J_analytic_unscreened = ReacNetJl._cached_network_jacobian!(Matrix{Float64}(undef, n, n), network, cache_unscreened, Y)
    scale_unscreened = maximum(abs, J_analytic_unscreened)
    @test maximum(abs, J_analytic_unscreened .- fd_jacobian(nothing)) <= 1.0e-6 * scale_unscreened

    # with screening the analytic Jacobian deliberately freezes the screening
    # factor with respect to Y (a valid Newton iteration matrix; the converged
    # answer is set by the exact residual), so the FD reference differs by the
    # small screening-derivative terms.
    cache = ReacNetJl._build_step_cache(network, rho, T9; screening=:weak)
    J_analytic = ReacNetJl._cached_network_jacobian!(Matrix{Float64}(undef, n, n), network, cache, Y)
    scale = maximum(abs, J_analytic)
    @test maximum(abs, J_analytic .- fd_jacobian(:weak)) <= 1.0e-2 * scale

    # repeated-reactant derivative: triple-alpha column d(dY_he4/dt)/dY_he4 = -9 * pref * Y^2
    he4 = network.species_index["he4"]
    c12 = network.species_index["c12"]
    triple_index = findfirst(r -> r.reactants == ["he4", "he4", "he4"], network.reactions)
    prefactor = cache.prefactors[triple_index]
    screen = ReacNetJl._cached_screening_multiplier(
        network.compiled_reactions[triple_index],
        ReacNetJl._screening_zeta_scale(cache, network, Y),
    )
    expected_he4_he4 = -9.0 * prefactor * screen * Y[he4]^2
    @test isapprox(J_analytic[he4, he4], expected_he4_he4; rtol=1.0e-10)
    @test isapprox(J_analytic[c12, he4], 3.0 * prefactor * screen * Y[he4]^2; rtol=1.0e-10)

    # full solves: analytic and finite-difference Jacobians must agree
    tspan = (0.0, 2.0e4)
    for screening in (nothing, :weak)
        t_a, h_a, s_a = solve_network_adaptive(
            network, Y, tspan, 1.0, rho, T9;
            method=:backward_euler, screening=screening, return_stats=true,
            jacobian=:analytic, dt_max=50.0,
        )
        t_f, h_f, s_f = solve_network_adaptive(
            network, Y, tspan, 1.0, rho, T9;
            method=:backward_euler, screening=screening, return_stats=true,
            jacobian=:finite_difference, dt_max=50.0,
        )
        @test isapprox(h_a[end, :], h_f[end, :]; rtol=1.0e-6, atol=1.0e-25)
    end

    # with unbounded dt the controller can push I - dt*J to numerical
    # singularity along conserved directions; the solver must treat that as a
    # rejected step and finish instead of crashing
    times_unbounded, history_unbounded, stats_unbounded = solve_network_adaptive(
        network, Y, tspan, 1.0, rho, T9;
        method=:backward_euler, return_stats=true, jacobian=:analytic,
    )
    @test times_unbounded[end] ≈ tspan[2]
    @test all(isfinite, history_unbounded)

    # fixed-step path and explicit methods still run through the cache; the
    # step sizes respect the fixture's fastest burning timescale (~5e-8 s),
    # since fixed-step integration cannot leap it the way the adaptive
    # controller can
    times_fixed, history_fixed = solve_network(network, Y, (0.0, 1.0e-6), 1.0e-7, rho, T9; method=:backward_euler, screening=:weak)
    @test all(isfinite, history_fixed)
    times_rk4, history_rk4 = solve_network(network, Y, (0.0, 1.0e-7), 1.0e-8, rho, T9; method=:rk4, screening=:weak)
    @test all(isfinite, history_rk4)
    @test_throws ArgumentError solve_network(network, Y, (0.0, 1.0e-6), 1.0e-7, rho, T9; method=:backward_euler, jacobian=:bogus)
end

@testset "partition functions and Chugunov screening" begin
    # synthetic winvne fixture: title, packed grid line, directory, records
    winvne_lines = String[]
    push!(winvne_lines, "")
    push!(winvne_lines, "010015020030040050060070080090100150200250300350400450500600700800900100")
    append!(winvne_lines, ["    p", " ne20", " na21"])
    grid_ones = join(fill("  1.00000E+0", 8), "")
    grid_gs = ["  1.00000E+0", "  1.00000E+0", "  1.00000E+0", "  1.00000E+0",
               "  1.00000E+0", "  1.00000E+0", "  1.00000E+0", "  1.00000E+0"]
    push!(winvne_lines, "    p       1.000   1   0   0.5     7.289 ame11")
    append!(winvne_lines, [grid_ones, grid_ones, grid_ones])
    push!(winvne_lines, " ne20      20.000  10  10   0.0    -7.042 ame11")
    append!(winvne_lines, [grid_ones, grid_ones, grid_ones])
    # na21 with a mildly temperature-dependent G, G(10 GK) = 2
    push!(winvne_lines, " na21      21.000  11  10   1.5    -2.184 ame11")
    append!(winvne_lines, [grid_ones, grid_ones, join(fill("  1.00000E+0", 7), "") * "  2.00000E+0"])

    winvne_path = joinpath(mktempdir(), "winvne_fixture.dat")
    write(winvne_path, join(winvne_lines, "\n") * "\n")
    pf = read_winvne(winvne_path)

    @test pf.g0["p"] == 2.0
    @test pf.g0["ne20"] == 1.0
    @test pf.g0["na21"] == 4.0
    @test ReacNetJl._partition_function_at(pf, "na21", 10.0) == 2.0
    @test ReacNetJl._partition_function_at(pf, "na21", 0.05) == 1.0
    @test ReacNetJl._partition_function_at(pf, "sr90", 1.0) === nothing

    grid = ReacNetJl.STARLIB_T9_GRID
    forward = ReactionRateTable(4, ["p", "ne20"], ["na21"], "test", 2.431, grid, fill(1.0e3, length(grid)), ones(length(grid)))
    plain = generated_detailed_balance_reverse_table(forward)
    corrected = generated_detailed_balance_reverse_table(forward; partition_functions=pf)

    # spin factor ga*gb/gc = 2*1/4 = 0.5 everywhere; pf ratio 1/G_c
    i_low = findfirst(==(0.2), grid)
    @test corrected.rate[i_low] ≈ 0.5 * plain.rate[i_low] rtol = 1.0e-12
    i_high = findfirst(==(10.0), grid)
    @test corrected.rate[i_high] ≈ 0.5 * plain.rate[i_high] / 2.0 rtol = 1.0e-12

    # Chugunov screening: sanity plus slow-path/cached-path equivalence
    proton_alpha = ReactionRateTable(5, ["p", "f18"], ["he4", "o15"], "test", 2.882, grid, fill(6.0e4, length(grid)), ones(length(grid)))
    network = network_from_tables([forward, proton_alpha])
    Y = abundances_from_mass_fractions(
        network,
        Dict("p" => 0.7, "ne20" => 0.1, "na21" => 1.0e-6, "f18" => 1.0e-4, "he4" => 0.2, "o15" => 0.0);
        normalize=true,
    )
    rho, T9 = 1.0e3, 0.2

    multiplier = chugunov_screening_multiplier(network, network.reactions[1], Y, rho, T9)
    @test multiplier > 1.0
    @test isfinite(multiplier)
    # hotter plasma screens less
    @test chugunov_screening_multiplier(network, network.reactions[1], Y, rho, 5.0) < multiplier
    # weak and Chugunov agree in the weak-screening regime (hot, dilute)
    weak_value = weak_screening_multiplier(network, network.reactions[1], Y, 10.0, 2.0)
    chugunov_value = chugunov_screening_multiplier(network, network.reactions[1], Y, 10.0, 2.0)
    @test isapprox(log(weak_value), log(chugunov_value); rtol=0.4)

    reference = network_rhs(Y, network, rho, T9; screening=:chugunov)
    cache = ReacNetJl._build_step_cache(network, rho, T9; screening=:chugunov)
    fast = ReacNetJl._cached_network_rhs!(similar(Y), network, cache, Y)
    @test isapprox(fast, reference; rtol=1.0e-12, atol=1.0e-300)

    # solver accepts the new mode end to end with the analytic Jacobian; the
    # fixed steps respect the fixture's ~1e-6 s burning timescale
    times, history = solve_network(network, Y, (0.0, 1.0e-5), 1.0e-6, rho, T9; method=:backward_euler, screening=:chugunov)
    @test all(isfinite, history)
end

@testset "run_ppn one-call workflow" begin
    dir = mktempdir()
    trajectory_path = joinpath(dir, "trajectory.input")
    write(trajectory_path, """
        AGEUNIT = SEC
        TUNIT   = T9K
        RHOUNIT = CGS
        0.0    0.20   1.0e3
        5.0    0.25   1.2e3
        10.0   0.15   8.0e2
        """)
    abundance_path = joinpath(dir, "initial_abundance.dat")
    write(abundance_path, """
          1 PROT          7.0E-01
          8 o  17         1.0E-01
          9 f  18         1.0E-04
          2 he  4         2.0E-01
        """)

    grid = ReacNetJl.STARLIB_T9_GRID
    unit = ones(length(grid))
    tables = [
        ReactionRateTable(4, ["p", "o17"], ["f18"], "test", 5.607, grid, fill(2.0e2, length(grid)), unit),
        ReactionRateTable(5, ["p", "f18"], ["he4", "o15"], "test", 2.882, grid, fill(6.0e2, length(grid)), unit),
        ReactionRateTable(1, ["o15"], ["n15"], "testw", 2.754, grid, fill(4.5e-3, length(grid)), unit),
    ]

    result = run_ppn(trajectory_path, abundance_path; tables=tables, screening=:weak)
    @test result.times[end] ≈ 10.0
    @test haskey(result.final_mass_fractions, "p")
    @test result.solver_stats.accepted_steps > 0
    @test result.mass_fraction_drift.max_abs_drift < 1.0e-6
    total_final = sum(values(result.final_mass_fractions))
    @test isapprox(total_final, sum(values(result.initial_mass_fractions)); rtol=1.0e-6)
    @test result.reverse_summary.generated + result.reverse_summary.explicit + result.reverse_summary.missing > 0
    @test result.rate_policy_report === nothing
end

const HAS_FBDF = Base.find_package("OrdinaryDiffEqBDF") !== nothing
HAS_FBDF && @eval using OrdinaryDiffEqBDF

@testset "FBDF stiff solver extension" begin
    if !HAS_FBDF
        @info "OrdinaryDiffEqBDF not installed; skipping FBDF extension tests"
        @test_skip false
    else
        grid = ReacNetJl.STARLIB_T9_GRID
        unit = ones(length(grid))
        tables = [
            ReactionRateTable(4, ["p", "o17"], ["f18"], "test", 5.607, grid, fill(2.0e2, length(grid)), unit),
            ReactionRateTable(5, ["p", "f18"], ["he4", "o15"], "test", 2.882, grid, fill(6.0e2, length(grid)), unit),
            ReactionRateTable(1, ["o15"], ["n15"], "testw", 2.754, grid, fill(4.5e-3, length(grid)), unit),
        ]
        network = network_from_tables(tables)
        X0 = Dict("p" => 0.7, "o17" => 0.1, "f18" => 1.0e-4, "he4" => 0.2, "o15" => 0.0, "n15" => 0.0)
        Y0 = abundances_from_mass_fractions(network, X0; normalize=true)
        rho, T9 = 1.0e3, 0.2
        tspan = (0.0, 10.0)

        times_f, history_f, stats_f = solve_network_fbdf(network, Y0, tspan, rho, T9; screening=:weak)
        @test times_f[end] ≈ tspan[2]
        @test all(isfinite, history_f)
        # trace species may dip negative at the tolerance level; that noise
        # must stay far below any physical abundance
        @test all(>=(-1.0e-8), history_f)
        @test stats_f.accepted_steps > 0

        times_b, history_b, stats_b = solve_network_adaptive(
            network, Y0, tspan, 1.0e-4, rho, T9;
            method=:backward_euler, screening=:weak, return_stats=true,
            abundance_floor=1.0e-8, max_newton_iterations=80,
        )
        for (i, name) in pairs(network.species)
            a = history_f[end, i]
            b = history_b[end, i]
            max(a, b) > 1.0e-12 && @test isapprox(a, b; rtol=5.0e-3)
        end

        # chugunov screening works through the ODE path too
        times_c, history_c, stats_c = solve_network_fbdf(network, Y0, tspan, rho, T9; screening=:chugunov)
        @test all(isfinite, history_c)

        @test_throws ArgumentError solve_network_fbdf(network, Y0, tspan, rho, T9; screening=(n, r, y, d, t) -> 2.0)
    end
end

@testset "tabulated paper rates and policy override" begin
    dir = mktempdir()

    # winvne fixture for Q-values (p, ne20, na21 with real mass excesses)
    winvne_lines = String[]
    push!(winvne_lines, "")
    push!(winvne_lines, "010015020030040050060070080090100150200250300350400450500600700800900100")
    grid_ones = join(fill("  1.00000E+0", 8), "")
    for (header, _) in [
        ("    p       1.000   1   0   0.5     7.289 ame11", 1),
        (" ne20      20.000  10  10   0.0    -7.042 ame11", 1),
        (" na21      21.000  11  10   1.5    -2.184 ame11", 1),
        (" f18       18.000   9   9   1.0     0.873 ame11", 1),
        (" o17       17.000   8   9   2.5    -0.809 ame11", 1),
    ]
        push!(winvne_lines, header)
        append!(winvne_lines, [grid_ones, grid_ones, grid_ones])
    end
    winvne_path = joinpath(dir, "winvne.dat")
    write(winvne_path, join(winvne_lines, "\n") * "\n")
    pf = read_winvne(winvne_path)

    @test pf.mass_excess["p"] == 7.289
    @test reaction_q_value(pf, ["p", "ne20"], ["na21"]) ≈ 7.289 - 7.042 + 2.184 atol = 1.0e-9
    @test reaction_q_value(pf, ["p", "sr90"], ["na21"]) === nothing

    # two-column block (il01 style) + four-column block (NACRE style) + total variant
    rates_path = joinpath(dir, "tabulated.dat")
    write(rates_path, """
        # fixture
        reaction: p ne20 -> na21 ; label: 20Ne(p,g) ; table: 3 ; variant: standard
        0.1 1.0E-08
        0.2 5.0E-06
        1.0 2.0E-02
        end
        reaction: p o17 -> f18 ; label: 17O(p,g) ; table: 2 ; variant: standard
        0.1 2.0E-09 1.0E-09 4.0E-09
        0.2 3.0E-06 1.5E-06 6.0E-06
        1.0 4.0E-02 2.0E-02 8.0E-02
        end
        reaction: p mg25 -> al26 ; label: 25Mg(p,g) ; table: 4 ; variant: total
        0.1 1.0E-05
        end
        """)

    tables = read_tabulated_rates(rates_path; source="testtab", partition_functions=pf)
    @test length(tables) == 2
    ne20 = only(find_rate(tables, "20Ne(p,γ)21Na"))
    @test ne20.source == "testtab"
    @test ne20.q_value ≈ 2.431 atol = 1.0e-3
    @test interpolate_rate(ne20, 0.2) ≈ 5.0e-6
    @test all(==(1.0), ne20.factor_uncertainty)

    o17 = only(find_rate(tables, "17O(p,γ)18F"))
    @test interpolate_factor_uncertainty(o17, 0.2) ≈ 2.0
    @test o17.q_value ≈ 7.289 - 0.809 - 0.873 atol = 1.0e-9

    with_total = read_tabulated_rates(rates_path; source="testtab", include_total_variants=true)
    @test length(with_total) == 3

    # policy override: paper tables replace matching REACLIB fits; the first
    # paper table wins when two cover the same reaction
    reaclib_header(species, label, res, rev, q) =
        " "^5 * join([lpad(s, 5) for s in vcat(species, fill("", 6 - length(species)))], "") *
        " "^8 * lpad(label, 4) * res * rev * " "^3 * @sprintf("%12.5e", q)
    reaclib_coefficients(a) =
        join([@sprintf("%13.6e", x) for x in a[1:4]], "") * "\n" *
        join([@sprintf("%13.6e", x) for x in a[5:7]], "")
    fixture = "4\n" * reaclib_header(["p", "ne20", "na21"], "nacr", 'n', ' ', 2.431) * "\n" *
              reaclib_coefficients((1.0, 0, 0, 0, 0, 0, 0)) * "\n"
    sets_path = joinpath(dir, "mini_reaclib.dat")
    write(sets_path, fixture)
    sets = read_reaclib(sets_path)

    duplicate = ReactionRateTable(4, ["p", "ne20"], ["na21"], "secondtab", 2.431,
        ne20.T9, fill(9.9e9, length(ne20.T9)), ones(length(ne20.T9)))
    result = iliadis2002_rate_tables(sets; paper_tables=vcat(tables, [duplicate]))
    @test result.report.paper_overrides.replaced == 1
    @test result.report.paper_overrides.added == 1
    chosen = only(find_rate(result.tables, "20Ne(p,γ)21Na"))
    @test chosen.source == "testtab"
    @test interpolate_rate(chosen, 0.2) ≈ 5.0e-6
end

@testset "AME2020 mass excesses" begin
    dir = mktempdir()
    ame_path = joinpath(dir, "ame_fixture.txt")
    # verbatim rows from mass.mas20.txt, incl. an estimated (#) entry
    write(ame_path, """
        header line one
        1N-Z    N    Z   A  EL    O     MASS EXCESS
        0  1    1    0    1  n         8071.31806     0.00044       0.0        0.0     B-    782.3470     0.0004    1 008664.91590     0.00047
          -1    0    1    1 H          7288.971064    0.000013      0.0        0.0     B-      *                    1 007825.031898    0.000014
        0  0    1    1    2 H         13135.722895    0.000015   1112.2831     0.0002  B-      *                    2 014101.777844    0.000015
          -3    0    3    3 Li  -pp   28667#       2000#        -2267#       667#      B-      *                    3 030775#       2147#
        0  2   11    9   20 F         -17.463         3.116     6987.9782     0.1558   B-   7025.4550     3.1201   19 999981.252      3.345
        """)
    masses = read_ame_masses(ame_path)
    @test masses["n"] ≈ 8.07131806
    @test masses["p"] ≈ 7.288971064
    @test masses["d"] ≈ 13.135722895
    @test masses["li3"] ≈ 28.667          # estimated entry kept by default
    @test masses["f20"] ≈ -0.017463
    strict = read_ame_masses(ame_path; include_estimated=false)
    @test !haskey(strict, "li3")
    @test reaction_q_value(masses, ["p", "p"], ["d"]) ≈ 2 * 7.288971064 - 13.135722895 atol = 1.0e-9
end

@testset "CSV outputs (mass fractions, fluxes, network) and run_ppn(output_dir=...)" begin
    grid = ReacNetJl.STARLIB_T9_GRID
    unit = ones(length(grid))
    tables = [
        ReactionRateTable(4, ["p", "o17"], ["f18"], "test", 5.607, grid, fill(2.0e2, length(grid)), unit),
        ReactionRateTable(5, ["p", "f18"], ["he4", "o15"], "test", 2.882, grid, fill(6.0e2, length(grid)), unit),
        ReactionRateTable(1, ["o15"], ["n15"], "testw", 2.754, grid, fill(4.5e-3, length(grid)), unit),
    ]
    network = network_from_tables(tables)
    Y0 = abundances_from_mass_fractions(
        network, Dict("p" => 0.7, "o17" => 0.1, "f18" => 1.0e-4, "he4" => 0.2, "o15" => 0.0, "n15" => 0.0); normalize=true,
    )
    times, history = solve_network_adaptive(network, Y0, (0.0, 5.0), 1.0e-4, 1.0e3, 0.2; method=:backward_euler, screening=:weak)

    dir = mktempdir()
    mf_path = write_mass_fraction_csv(joinpath(dir, "mf.csv"), network, times, history)
    lines = readlines(mf_path)
    @test split(lines[1], ",") == vcat(["time_s"], network.species)
    @test length(lines) == length(times) + 1
    parsed_first_row = parse.(Float64, split(lines[2], ",")[2:end])
    X_first = mass_fractions_from_abundances(network, view(history, 1, :))
    @test isapprox(parsed_first_row, [X_first[name] for name in network.species]; rtol=1.0e-12)

    flux_history = reaction_flux_history(network, history, times, 1.0e3, 0.2; screening=:weak)
    flux_path = write_reaction_flux_csv(joinpath(dir, "flux.csv"), network, times, flux_history)
    flux_lines = readlines(flux_path)
    @test length(flux_lines) == length(times) + 1
    # a reaction label with a comma (e.g. "17O(p,γ)18F") must round-trip as one quoted field
    header_fields = split(flux_lines[1], r",(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)")
    @test length(header_fields) == length(network.reactions) + 1
    @test any(f -> occursin("(p,", f), header_fields)

    net_path = write_network_csv(joinpath(dir, "network.csv"), network)
    net_lines = readlines(net_path)
    @test length(net_lines) == length(network.reactions) + 1
    @test net_lines[1] == "index,reaction,reactants,products,chapter,source,q_value_mev"

    # end-to-end via run_ppn(output_dir=...)
    traj_path = joinpath(dir, "trajectory.input")
    write(traj_path, "AGEUNIT = SEC\nTUNIT = T9K\nRHOUNIT = CGS\n0.0 0.20 1.0e3\n5.0 0.20 1.0e3\n")
    ab_path = joinpath(dir, "initial_abundance.dat")
    write(ab_path, "  1 PROT 7.0E-01\n  8 o  17 1.0E-01\n  9 f  18 1.0E-04\n  2 he  4 2.0E-01\n")
    out_dir = joinpath(dir, "run_ppn_out")

    result = run_ppn(traj_path, ab_path; tables=tables, screening=:weak, output_dir=out_dir)
    @test result.output_files.mass_fractions == joinpath(out_dir, "mass_fractions.csv")
    @test isfile(result.output_files.mass_fractions)
    @test isfile(result.output_files.reaction_fluxes)
    @test isfile(result.output_files.network)

    no_flux_result = run_ppn(traj_path, ab_path; tables=tables, screening=:weak, output_dir=out_dir, write_reaction_fluxes=false)
    @test no_flux_result.output_files.reaction_fluxes === nothing

    no_output_result = run_ppn(traj_path, ab_path; tables=tables, screening=:weak)
    @test no_output_result.output_files === nothing
end

@testset "unseen isotope requires zero special-casing" begin
    # te132 (tellurium-132) appears nowhere else in this codebase: it is not
    # in any rate table fixture, not in _SPECIAL_SPECIES, not in any prior
    # test. If species_from_name/normalize_species_name ever regress into
    # needing per-isotope entries, this is the isotope that would catch it.
    s = species_from_name("132Te")
    @test s == Species("te132", 52, 132)
    @test normalize_species_name("132Te") == "te132"
    @test species_from_name("te132") == s

    # round-trips for a light, a mid-mass, and a super-heavy element, none of
    # which are referenced by any rate table or fixture in this suite either
    @test species_from_name("6He") == Species("he6", 2, 6)
    @test species_from_name("91Zr") == Species("zr91", 40, 91)
    @test species_from_name("287Fl") == Species("fl287", 114, 287)

    # the only physically-motivated exception to "(Z,A) fully identifies a
    # species" is a long-lived isomer, and it is confined to 26Al
    ground = species_from_name("26Al")
    isomer = species_from_name("26Alm")
    @test (ground.Z, ground.A) == (isomer.Z, isomer.A) == (13, 26)
    @test ground.name != isomer.name
end

const HAS_FORWARDDIFF = Base.find_package("ForwardDiff") !== nothing
HAS_FORWARDDIFF && @eval using ForwardDiff

@testset "type-generic RHS/flux differentiate under ForwardDiff (feature spec Tier 0 #15)" begin
    if !HAS_FORWARDDIFF
        @info "ForwardDiff not installed; skipping AD genericity tests"
        @test_skip false
    else
        # Confirms the *kind* of bug this genericity pass fixes: a bare
        # Float64(::Dual) conversion silently discards derivative
        # information rather than erroring loudly -- ForwardDiff instead
        # disallows it outright, which is why every `Float64(rho)`/
        # `Float64(rate_multipliers[r])` cast that used to sit in the
        # flux/RHS hot path had to go, not just be "made to work".
        @test_throws MethodError Float64(ForwardDiff.Dual(1.0, 1.0))

        grid = ReacNetJl.STARLIB_T9_GRID
        unit = ones(length(grid))
        capture = ReactionRateTable(4, ["p", "o17"], ["f18"], "test", 5.607, grid, fill(2.0e2, length(grid)), unit)
        proton_alpha = ReactionRateTable(5, ["p", "f18"], ["he4", "o15"], "test", 2.882, grid, fill(6.0e2, length(grid)), unit)
        network = network_from_tables([capture, proton_alpha])
        Y = abundances_from_mass_fractions(
            network, Dict("p" => 0.6, "o17" => 0.1, "f18" => 1.0e-4, "he4" => 0.2, "o15" => 0.0); normalize=true,
        )
        rho, T9 = 800.0, 0.22

        # dY/dt is linear-homogeneous in rate_multipliers (each reaction's
        # flux scales with only its own multiplier), so the AD Jacobian at
        # any point must reproduce dYdt(m) exactly at a *different* m -- a
        # model-free correctness check that doesn't re-derive the flux
        # formula by hand.
        f(m) = network_rhs(Y, network, rho, T9; rate_multipliers=m)
        m0 = [1.0, 1.0]
        J = ForwardDiff.jacobian(f, m0)
        m_test = [2.3, 0.4]
        @test isapprox(J * m_test, f(m_test); rtol=1.0e-10)
        @test all(iszero, f(zeros(2)))

        # same check through the cached/compiled path the real solvers use
        # (backward_euler.jl via _step_cache_at/_cached_network_rhs!) -- this
        # is the path where `Float64(rate_multipliers[r])` used to sit.
        function f_cached(m)
            cache = ReacNetJl._build_step_cache(network, rho, T9; rate_multipliers=m)
            buffer = zeros(eltype(m), length(Y))
            return ReacNetJl._cached_network_rhs!(buffer, network, cache, Y)
        end
        J_cached = ForwardDiff.jacobian(f_cached, m0)
        @test isapprox(J_cached, J; rtol=1.0e-10)

        # per-reaction fluxes: off-diagonal terms must vanish exactly, since
        # each reaction's flux depends only on its own multiplier
        flux_of(m) = reaction_fluxes(network, Y, rho, T9; rate_multipliers=m)
        J_flux = ForwardDiff.jacobian(flux_of, m0)
        base_flux = flux_of(m0)
        for r in 1:2, s in 1:2
            expected = r == s ? base_flux[r] : 0.0
            @test isapprox(J_flux[r, s], expected; atol=1.0e-12, rtol=1.0e-10)
        end
    end
end

const HAS_KLU = Base.find_package("KLU") !== nothing
HAS_KLU && @eval using KLU

@testset "sparse Jacobian pattern and KLU-backed solve (feature spec Tier 0 #2/#3)" begin
    using SparseArrays

    grid = ReacNetJl.STARLIB_T9_GRID
    unit_uncertainty = ones(length(grid))
    capture = ReactionRateTable(4, ["p", "o17"], ["f18"], "test", 5.607, grid, fill(2.0e4, length(grid)), unit_uncertainty)
    proton_alpha = ReactionRateTable(5, ["p", "f18"], ["he4", "o15"], "test", 2.882, grid, fill(6.0e4, length(grid)), unit_uncertainty)
    decay = ReactionRateTable(1, ["o15"], ["n15"], "testw", 2.754, grid, fill(4.5e-3, length(grid)), unit_uncertainty)
    triple_alpha = ReactionRateTable(8, ["he4", "he4", "he4"], ["c12"], "test", 7.275, grid, fill(1.0e-9, length(grid)), unit_uncertainty)
    network = network_from_tables([capture, proton_alpha, decay, triple_alpha])

    X0 = Dict(
        "p" => 0.6, "o17" => 0.05, "f18" => 1.0e-4, "he4" => 0.3,
        "o15" => 1.0e-5, "n15" => 0.0, "c12" => 0.05,
    )
    Y = abundances_from_mass_fractions(network, X0; normalize=true)
    rho, T9 = 500.0, 0.25
    n = length(network.species)

    # structural pattern: full diagonal always present, exact nonzero count
    # for this fixture (capture: 3 rows x 2 cols; proton_alpha: 4x2; decay:
    # 2x1; triple_alpha's single distinct reactant column he4: 2x1 -- 18
    # unique (row,col) pairs after deduping against the diagonal), and the
    # pattern must not be a placeholder (all-dense or all-empty).
    pattern = network.jacobian_sparsity
    @test all(pattern[i, i] for i in 1:n)
    @test nnz(pattern) == 18
    @test 0 < nnz(pattern) < n^2

    cache = ReacNetJl._build_step_cache(network, rho, T9)
    J_dense = ReacNetJl._cached_network_jacobian!(Matrix{Float64}(undef, n, n), network, cache, Y)

    # every dense-analytic nonzero must be covered by the precomputed pattern
    # (the reverse -- a few structurally-nonzero-but-numerically-zero entries
    # -- is legitimate and not tested against)
    for i in 1:n, j in 1:n
        J_dense[i, j] == 0.0 || @test pattern[i, j]
    end

    # sparse assembly must match the dense analytic Jacobian exactly -- same
    # per-reaction product-rule formula, only the write target differs
    J_sparse = ReacNetJl._cached_network_jacobian_sparse!(sparse_jacobian_prototype(network), network, cache, Y)
    @test Matrix(J_sparse) == J_dense

    # out-of-pattern (row, col) must error, not silently no-op -- guards
    # against a future sparsity-pattern bug masking a missing Jacobian term
    zero_positions = [(i, j) for i in 1:n, j in 1:n if !pattern[i, j]]
    if !isempty(zero_positions)
        i0, j0 = first(zero_positions)
        @test_throws ArgumentError ReacNetJl._sparse_nzval_index(J_sparse, i0, j0)
    end

    if !HAS_KLU
        @info "KLU not installed; skipping sparse-solve physical-accuracy tests"
        @test_skip false
    else
        # physical-accuracy check: jacobian=:sparse (KLU) must reproduce the
        # already-validated jacobian=:analytic (dense LU) trajectory, not
        # merely "run without erroring" -- both the adaptive controller and
        # the fixed-step driver, since both route through _backward_euler_step.
        tspan = (0.0, 2.0e4)
        t_dense, h_dense = solve_network_adaptive(network, Y, tspan, 1.0, rho, T9; method=:backward_euler, jacobian=:analytic, dt_max=50.0)
        t_sparse, h_sparse = solve_network_adaptive(network, Y, tspan, 1.0, rho, T9; method=:backward_euler, jacobian=:sparse, dt_max=50.0)
        @test t_sparse[end] ≈ t_dense[end]
        @test all(isfinite, h_sparse)
        @test isapprox(h_sparse[end, :], h_dense[end, :]; rtol=1.0e-8, atol=1.0e-25)

        tf_dense, hf_dense = solve_network(network, Y, (0.0, 1.0e-6), 1.0e-7, rho, T9; method=:backward_euler)
        tf_sparse, hf_sparse = solve_network(network, Y, (0.0, 1.0e-6), 1.0e-7, rho, T9; method=:backward_euler, jacobian=:sparse)
        @test isapprox(hf_sparse[end, :], hf_dense[end, :]; rtol=1.0e-8, atol=1.0e-25)

        @test_throws ArgumentError solve_network(network, Y, (0.0, 1.0e-6), 1.0e-7, rho, T9; method=:backward_euler, jacobian=:bogus)
    end
end
