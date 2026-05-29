using ReacNetJl
using Test

@testset "species names" begin
    @test normalize_species_name("18F") == "f18"
    @test normalize_species_name("4He") == "he4"
    @test normalize_species_name("α") == "he4"
    @test normalize_species_name("γ") == "gamma"
    @test normalize_species_name("p") == "p"

    @test species_from_name("18F") == Species("f18", 9, 18)
    @test species_from_name("4He") == Species("he4", 2, 4)
    @test species_from_name("p") == Species("p", 1, 1)
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
    @test_throws ArgumentError profiles.T9(-0.1)

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

    @test sum(abundances_from_mass_fractions(network, Dict("p" => 7.0, "he4" => 3.0); normalize=true) .* [1, 18, 4, 15]) ≈ 1.0
    @test_throws ArgumentError abundances_from_mass_fractions(network, Dict("ne19" => 1.0))
    @test_throws ArgumentError abundances_from_mass_fractions(network, Dict("p" => 0.5); check_sum=true)
    @test_throws ArgumentError mass_fractions_from_abundances(network, [0.7])

    Y0 = [0.7, 1.0e-5, 0.0, 0.0]
    @test reaction_string(reaction) == "f18(p,he4)o15"

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

    edges = reaction_edges(network)
    @test length(edges) == 4
    @test (reaction_index=1, reaction="f18(p,he4)o15", from="f18", to="o15") in edges
    @test (reaction_index=1, reaction="f18(p,he4)o15", from="p", to="he4") in edges

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
    @test reaction_string(decay) == "o15(β+)n15"
    decay_conservation = reaction_conservation(decay)
    @test decay_conservation.conserves_A
    @test !decay_conservation.conserves_Z
    @test decay_conservation.is_weak_decay
    @test decay_conservation.valid_nuclear_bookkeeping
    decay_network = ReactionNetwork(["o15", "n15"], [decay])
    @test network_validation_report(decay_network).valid

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

    flux_history = reaction_flux_history(network, history, times, 1000.0, 0.5)
    @test size(flux_history) == (length(times), length(network.reactions))
    @test flux_history[1, 1] ≈ fluxes[1]
    @test integrated_fluxes([0.0, 2.0], [1.0 2.0; 3.0 4.0]) ≈ [4.0, 6.0]
    @test length(integrated_fluxes(times, flux_history)) == length(network.reactions)
    @test integrated_fluxes(times, flux_history)[1] > 0.0
    @test_throws ArgumentError reaction_flux_history(network, history, times[1:end-1], 1000.0, 0.5)
    @test_throws ArgumentError integrated_fluxes([0.0], [1.0 2.0])
    @test_throws ArgumentError integrated_fluxes([0.0, -1.0], reshape([1.0, 2.0], 2, 1))

    _, boosted_history = solve_network(network, Y0, (0.0, 1.0e-3), 1.0e-4, 1000.0, 0.5; method=:rk4, rate_multipliers=[2.0])
    @test boosted_history[end, network.species_index["f18"]] < history[end, network.species_index["f18"]]

    _, p_history = solve_network(network, Y0, (0.0, 1.0e-3), 1.0e-4, 1000.0, 0.5; method=:rk4, rate_p_values=[1.0])
    @test p_history[end, network.species_index["f18"]] < history[end, network.species_index["f18"]]

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
    @test_throws ArgumentError solve_network_adaptive(network, Y0, (0.0, 1.0), 0.1, 1000.0, 0.5; max_fractional_change=0.0)
end

@testset "read synthetic STARLIB table" begin
    path = tempname()
    open(path, "w") do io
        println(io, "2 p f18 ne19 testsource 3.529")
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
    @test table.chapter == 2
    @test table.reactants == ["p", "f18"]
    @test table.products == ["ne19"]
    @test table.source == "testsource"
    @test table.q_value == 3.529
    @test length(table.T9) == 60
    @test table.T9[1] == 0.01
    @test table.rate[end] == 6.0e-11
    @test table.factor_uncertainty[end] == 2.0

    matches = find_rate(tables, "18F(p,γ)19Ne")
    @test matches == tables
    @test find_rate(tables, "18F(p,γ)19Ne"; source="missing") == ReactionRateTable[]

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
end
