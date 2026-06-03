function usage()
    return """
    Usage:
      julia --project=. examples/single_zone_nova_ppn.jl [--jobs N] [--output-stride N] [--screening weak|none] [--decay-time S] [--stop-time S] [--stop-temperature T9] [--stop-hydrogen X] [--flux-report N]

    Options:
      --jobs N           Number of Julia threads to use for threaded Jacobian and DAT output work.
      --output-stride N  Write every Nth trajectory output state. Default: 1, or NOVA_PPN_OUTPUT_STRIDE.
      --screening MODE   Screening model: weak or none. Default: weak, or NOVA_PPN_SCREENING.
      --no-screening     Alias for --screening none.
      --decay-time S     Evolve weak decays for S seconds after the final state and write iso_massfDECAY.DAT.
      --stop-time S      Stop the network at absolute trajectory time S seconds.
      --stop-temperature T9
                         Stop on the post-peak cooling branch when T9 falls to this value.
      --stop-hydrogen X  Report/write the first state where active X(PROT) falls to X.
      --flux-report N    Print integrated production/destruction channels for N worst residual isotopes. Default: 6.
      --help             Show this message.
    """
end

function screening_model_from_name(value)
    name = lowercase(strip(string(value)))
    if name == "weak"
        return :weak
    elseif name in ("none", "off", "false", "no")
        return nothing
    end

    throw(ArgumentError("unsupported screening '$value'; use weak or none"))
end

function screening_label(screening_model)
    screening_model == :weak && return "weak"
    (screening_model === nothing || screening_model === false) && return "none"
    return string(screening_model)
end

function parse_optional_float_env(name::AbstractString)
    value = strip(get(ENV, name, ""))
    isempty(value) && return nothing
    return parse(Float64, value)
end

function parse_cli_args(args)
    jobs = 1
    output_stride = parse(Int, get(ENV, "NOVA_PPN_OUTPUT_STRIDE", "1"))
    screening = screening_model_from_name(get(ENV, "NOVA_PPN_SCREENING", "weak"))
    decay_time_s = parse(Float64, get(ENV, "NOVA_PPN_DECAY_TIME", "0.0"))
    stop_time_s = parse_optional_float_env("NOVA_PPN_STOP_TIME")
    stop_temperature_T9 = parse_optional_float_env("NOVA_PPN_STOP_TEMPERATURE")
    stop_hydrogen_X = parse_optional_float_env("NOVA_PPN_STOP_HYDROGEN")
    flux_report_count = parse(Int, get(ENV, "NOVA_PPN_FLUX_REPORT", "6"))
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--help" || arg == "-h"
            println(usage())
            exit(0)
        elseif arg == "--jobs"
            i < length(args) || throw(ArgumentError("--jobs requires an integer value"))
            jobs = parse(Int, args[i + 1])
            i += 2
        elseif startswith(arg, "--jobs=")
            jobs = parse(Int, split(arg, "="; limit=2)[2])
            i += 1
        elseif arg == "--output-stride"
            i < length(args) || throw(ArgumentError("--output-stride requires an integer value"))
            output_stride = parse(Int, args[i + 1])
            i += 2
        elseif startswith(arg, "--output-stride=")
            output_stride = parse(Int, split(arg, "="; limit=2)[2])
            i += 1
        elseif arg == "--screening"
            i < length(args) || throw(ArgumentError("--screening requires weak or none"))
            screening = screening_model_from_name(args[i + 1])
            i += 2
        elseif startswith(arg, "--screening=")
            screening = screening_model_from_name(split(arg, "="; limit=2)[2])
            i += 1
        elseif arg == "--no-screening"
            screening = nothing
            i += 1
        elseif arg == "--decay-time"
            i < length(args) || throw(ArgumentError("--decay-time requires a non-negative number of seconds"))
            decay_time_s = parse(Float64, args[i + 1])
            i += 2
        elseif startswith(arg, "--decay-time=")
            decay_time_s = parse(Float64, split(arg, "="; limit=2)[2])
            i += 1
        elseif arg == "--stop-time"
            i < length(args) || throw(ArgumentError("--stop-time requires seconds"))
            stop_time_s = parse(Float64, args[i + 1])
            i += 2
        elseif startswith(arg, "--stop-time=")
            stop_time_s = parse(Float64, split(arg, "="; limit=2)[2])
            i += 1
        elseif arg == "--stop-temperature"
            i < length(args) || throw(ArgumentError("--stop-temperature requires a T9 value"))
            stop_temperature_T9 = parse(Float64, args[i + 1])
            i += 2
        elseif startswith(arg, "--stop-temperature=")
            stop_temperature_T9 = parse(Float64, split(arg, "="; limit=2)[2])
            i += 1
        elseif arg == "--stop-hydrogen"
            i < length(args) || throw(ArgumentError("--stop-hydrogen requires a mass fraction"))
            stop_hydrogen_X = parse(Float64, args[i + 1])
            i += 2
        elseif startswith(arg, "--stop-hydrogen=")
            stop_hydrogen_X = parse(Float64, split(arg, "="; limit=2)[2])
            i += 1
        elseif arg == "--flux-report"
            i < length(args) || throw(ArgumentError("--flux-report requires an integer count"))
            flux_report_count = parse(Int, args[i + 1])
            i += 2
        elseif startswith(arg, "--flux-report=")
            flux_report_count = parse(Int, split(arg, "="; limit=2)[2])
            i += 1
        else
            throw(ArgumentError("unknown argument '$arg'\n$(usage())"))
        end
    end

    jobs >= 1 || throw(ArgumentError("--jobs must be at least 1"))
    output_stride >= 1 || throw(ArgumentError("--output-stride must be at least 1"))
    decay_time_s >= 0.0 || throw(ArgumentError("--decay-time must be non-negative"))
    stop_time_s === nothing || stop_time_s >= 0.0 || throw(ArgumentError("--stop-time must be non-negative"))
    stop_temperature_T9 === nothing || stop_temperature_T9 > 0.0 || throw(ArgumentError("--stop-temperature must be positive"))
    stop_hydrogen_X === nothing || stop_hydrogen_X >= 0.0 || throw(ArgumentError("--stop-hydrogen must be non-negative"))
    flux_report_count >= 0 || throw(ArgumentError("--flux-report must be non-negative"))
    return (
        jobs=jobs,
        output_stride=output_stride,
        screening=screening,
        decay_time_s=decay_time_s,
        stop_time_s=stop_time_s,
        stop_temperature_T9=stop_temperature_T9,
        stop_hydrogen_X=stop_hydrogen_X,
        flux_report_count=flux_report_count,
    )
end

function maybe_relaunch_with_threads(cli)
    cli.jobs <= Base.Threads.nthreads() && return
    get(ENV, "REACNETJL_PPN_RELAUNCHED", "0") == "1" && return

    project = Base.active_project()
    project_dir = project === nothing ? dirname(@__DIR__) : dirname(project)
    env = copy(ENV)
    env["REACNETJL_PPN_RELAUNCHED"] = "1"
    env["JULIA_NUM_THREADS"] = string(cli.jobs)
    cmd = `$(Base.julia_cmd()) --project=$project_dir --threads=$(cli.jobs) $(PROGRAM_FILE) $(ARGS)`
    run(setenv(cmd, env))
    exit()
end

const CLI = parse_cli_args(ARGS)
maybe_relaunch_with_threads(CLI)

using Printf
using ReacNetJl

#=
    Full single-zone nova PPN example

This script evolves the current expanded nova network over the project-root
`trajectory.input` using `initial_abundance.DAT` or `initial_abundance.dat`.
It writes one `iso_massfXXXXX.DAT` mass-fraction file per saved trajectory
state, plus a wide CSV containing every isotope column. With the current
805-row trajectory, the last default file is `iso_massf00804.DAT`.

Set `--jobs 8` to run threaded Jacobian construction and parallel DAT output
writing with eight Julia threads. Set `--output-stride 1000` for quick checks.
The default stride is 1, which writes every trajectory state.
=#

const YEAR_SECONDS = 365.25 * 24.0 * 60.0 * 60.0
const AVOGADRO_CGS = 6.02214076e23
const MASS_FRACTION_FLOOR = 1.0e-99

labels = [
    # Hot CNO-ish flow and decays
    "12C(p,γ)13N",
    "13N(β+)13C",
    "13C(p,γ)14N",
    "14N(p,γ)15O",
    "15O(β+)15N",
    "15N(p,α)12C",
    "15N(p,γ)16O",
    "16O(p,γ)17F",
    "17F(β+)17O",
    "17O(p,α)14N",
    "17O(p,γ)18F",
    "18F(β+)18O",
    "18O(p,α)15N",
    "18O(p,γ)19F",
    "18F(p,α)15O",
    "18F(p,γ)19Ne",
    "19Ne(β+)19F",
    "19F(p,α)16O",
    "19F(p,γ)20Ne",

    # NeNa cycle fragments
    "20Ne(p,γ)21Na",
    "21Na(β+)21Ne",
    "21Na(p,γ)22Mg",
    "22Mg(β+)22Na",
    "21Ne(p,γ)22Na",
    "22Na(β+)22Ne",
    "22Na(p,γ)23Mg",
    "23Mg(β+)23Na",
    "22Ne(p,γ)23Na",
    "23Na(p,α)20Ne",
    "23Na(p,γ)24Mg",

    # MgAl fragments
    "24Mg(p,γ)25Al",
    "25Al(β+)25Mg",
    "25Al(p,γ)26Si",
    "26Si(β+)26Al",
    "26Si(β+)26Al*",
    "25Mg(p,γ)26Al",
    "26Mg(p,γ)27Al",
    "26Al*(β+)26Mg",
    "26Al(p,γ)27Si",
    "26Al*(p,γ)27Si",
    "27Si(β+)27Al",
    "27Al(p,α)24Mg",
    "27Al(p,γ)28Si",

    # Si-Ca extension
    "28Si(p,γ)29P",
    "29P(β+)29Si",
    "29Si(p,γ)30P",
    "30P(β+)30Si",
    "30P(p,γ)31S",
    "31S(β+)31P",
    "30Si(p,γ)31P",
    "31P(p,γ)32S",
    "32S(p,γ)33Cl",
    "33Cl(β+)33S",
    "33S(p,γ)34Cl",
    "34Cl(β+)34S",
    "34S(p,γ)35Cl",
    "35Cl(p,γ)36Ar",
    "36Ar(p,γ)37K",
    "37K(β+)37Ar",
    "37Ar(p,γ)38K",
    "38K(β+)38Ar",
    "38Ar(p,γ)39K",
    "39K(p,γ)40Ca",

    # Ca-Fe/Ni seed extension
    "40Ca(p,γ)41Sc",
    "41Sc(β+)41Ca",
    "41Ca(p,γ)42Sc",
    "42Ca(p,γ)43Sc",
    "42Sc(β+)42Ca",
    "43Sc(β+)43Ca",
    "43Ca(p,γ)44Sc",
    "44Ca(p,γ)45Sc",
    "44Sc(β+)44Ca",
    "45Sc(p,γ)46Ti",
    "45Sc(p,α)42Ca",
    "46Ti(p,γ)47V",
    "47V(β+)47Ti",
    "47Ti(p,γ)48V",
    "48Ti(p,γ)49V",
    "48V(β+)48Ti",
    "49V(β+)49Ti",
    "49Ti(p,γ)50V",
    "50Ti(p,γ)51V",
    "50V(p,γ)51Cr",
    "51V(p,γ)52Cr",
    "51Cr(β+)51V",
    "52Cr(p,γ)53Mn",
    "53Mn(p,γ)54Fe",
    "53Mn(β+)53Cr",
    "54Fe(p,γ)55Co",
    "55Co(β+)55Fe",
    "55Mn(p,γ)56Fe",
    "56Fe(p,γ)57Co",
    "57Co(β+)57Fe",
    "57Fe(p,γ)58Co",
    "58Ni(p,γ)59Cu",
    "59Cu(β+)59Ni",
    "59Co(p,γ)60Ni",
]

function required_file(paths::Vector{String}, description::String)
    index = findfirst(isfile, paths)
    index === nothing && throw(ArgumentError("missing $description; checked: $(join(paths, ", "))"))
    return paths[index]
end

function iso_label(name::AbstractString)
    normalized = normalize_species_name(name)
    normalized == "n" && return "NEUT"
    normalized == "p" && return "PROT"
    normalized == "d" && return "H   2"
    normalized == "t" && return "H   3"
    normalized == "he3" && return "HE  3"
    normalized == "he4" && return "HE  4"

    m = match(r"^([a-z]+)\*(\d+)$", normalized)
    if m !== nothing
        return uppercase(m.captures[1]) * "*" * string(species_from_name(normalized).A)
    end

    info = species_from_name(normalized)
    symbol_match = match(r"^([a-z]+)\d+$", normalized)
    symbol = symbol_match === nothing ? uppercase(normalized) : uppercase(symbol_match.captures[1])
    return @sprintf("%-2s%3d", symbol, info.A)
end

function output_record(name::AbstractString)
    normalized = normalize_species_name(name)
    info = species_from_name(normalized)
    isom = occursin("*", normalized) ? 2 : 1
    return (name=normalized, Z=info.Z, A=info.A, isom=isom, label=iso_label(normalized))
end

function parse_template_records(path::AbstractString)
    records = NamedTuple[]
    open(path, "r") do io
        for raw_line in eachline(io)
            line = strip(raw_line)
            (isempty(line) || startswith(line, "H ") || startswith(line, "#")) && continue
            fields = split(line)
            length(fields) < 6 && continue

            label = join(fields[6:end], "")
            name = if label == "NEUT"
                "n"
            elseif label == "PROT"
                "p"
            else
                normalize_species_name(label)
            end

            push!(
                records,
                (
                    name=name,
                    Z=Int(round(parse(Float64, fields[2]))),
                    A=Int(round(parse(Float64, fields[3]))),
                    isom=Int(round(parse(Float64, fields[4]))),
                    label=iso_label(name),
                ),
            )
        end
    end
    return records
end

function build_output_records(template_path, network::ReactionNetwork, X_file::Dict{String,Float64})
    records = template_path === nothing ? NamedTuple[] : parse_template_records(template_path)
    seen = Set(record.name for record in records)

    names = sort!(collect(union(Set(network.species), Set(keys(X_file)))); by=name -> begin
        info = species_from_name(name)
        (info.Z, info.A, name)
    end)
    for name in names
        if !(name in seen)
            push!(records, output_record(name))
            push!(seen, name)
        end
    end
    return records
end

function saved_indices(nstates::Int, stride::Int)
    stride >= 1 || throw(ArgumentError("NOVA_PPN_OUTPUT_STRIDE must be at least 1"))
    indices = collect(1:stride:nstates)
    indices[end] == nstates || push!(indices, nstates)
    return indices
end

function interpolated_state(times, history, time_s::Float64)
    time_s <= first(times) && return collect(history[1, :])
    time_s >= last(times) && return collect(history[end, :])

    i = searchsortedlast(times, time_s)
    if times[i] == time_s
        return collect(history[i, :])
    end

    weight = (time_s - times[i]) / (times[i + 1] - times[i])
    return (1.0 - weight) .* collect(history[i, :]) .+ weight .* collect(history[i + 1, :])
end

function output_values(network::ReactionNetwork, Y, X_file::Dict{String,Float64}, records)
    active_X = mass_fractions_from_abundances(network, Y)
    return [max(get(active_X, record.name, get(X_file, record.name, 0.0)), MASS_FRACTION_FLOOR) for record in records]
end

function output_values_from_mass_fractions(X::Dict{String,Float64}, records)
    return [max(get(X, record.name, 0.0), MASS_FRACTION_FLOOR) for record in records]
end

function density_markers(rho::Float64, values, records)
    index = Dict(record.name => i for (i, record) in pairs(records))
    x(name) = haskey(index, name) ? values[index[name]] : 0.0
    densn = rho * AVOGADRO_CGS * x("n")
    densp = rho * AVOGADRO_CGS * x("p")
    densa = rho * AVOGADRO_CGS * x("he4") / 4.0
    return densn, densp, densa
end

function write_iso_massf_file(
    path::AbstractString,
    step_number::Int,
    time_s::Float64,
    previous_time_s::Float64,
    T9::Float64,
    rho::Float64,
    epsilon::Float64,
    records,
    values,
)
    mkpath(dirname(path))
    dzeit = (time_s - previous_time_s) / YEAR_SECONDS
    agej = time_s / YEAR_SECONDS
    densn, densp, densa = density_markers(rho, values, records)

    open(path, "w") do io
        println(io, "H NUM     Z    A ISOM  ABUNDANCE_MF  ISOTP")
        @printf(io, " # mod %8d dzeit %9.2E agej %11.4E\n", step_number, dzeit, agej)
        @printf(io, " # t9 = %12.3E rho = %12.3E\n", T9, rho)
        @printf(io, " # densn  %.5E\n", densn)
        @printf(io, " # densp  %.5E\n", densp)
        @printf(io, " # densa  %.5E\n", densa)
        @printf(io, " # total_energy_flux_[erg/(g*s)]  %.5E\n", epsilon)
        for (row, record) in pairs(records)
            @printf(
                io,
                "%5d %5.0f. %4.0f. %3d %12.5E  %s\n",
                row,
                Float64(record.Z),
                Float64(record.A),
                record.isom,
                values[row],
                record.label,
            )
        end
    end
end

function write_csv(
    path::AbstractString,
    indices,
    output_times,
    solver_times,
    history,
    profiles,
    network::ReactionNetwork,
    X_file::Dict{String,Float64},
    records,
    screening_model,
)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(vcat(["step", "time_s", "T9", "rho", "epsilon_nuc"], [record.name for record in records]), ","))
        for output_index in indices
            step_number = output_index - 1
            time_s = output_times[output_index]
            T9 = profiles.T9(time_s)
            rho = profiles.rho(time_s)
            Y = interpolated_state(solver_times, history, time_s)
            epsilon = energy_generation_rate(network, Y, rho, T9; screening=screening_model)
            values = output_values(network, Y, X_file, records)
            row = String[
                string(step_number),
                @sprintf("%.16E", time_s),
                @sprintf("%.16E", T9),
                @sprintf("%.16E", rho),
                @sprintf("%.16E", epsilon),
            ]
            append!(row, [@sprintf("%.16E", value) for value in values])
            println(io, join(row, ","))
        end
    end
end

function clean_previous_outputs(output_dir::AbstractString)
    isdir(output_dir) || return
    for path in readdir(output_dir; join=true)
        base = basename(path)
        if occursin(r"^iso_massf\d+\.DAT$", base) ||
           base == "iso_massfDECAY.DAT" ||
           base == "mass_fractions.csv" ||
           base == "comparison_final_vs_iliadis.csv" ||
           base == "comparison_postdecay_vs_iliadis.csv"
            rm(path; force=true)
        end
    end
end

function write_saved_iso_state(
    output_index::Int,
    output_dir::AbstractString,
    output_times,
    solver_times,
    history,
    profiles,
    network::ReactionNetwork,
    X_file::Dict{String,Float64},
    records,
    screening_model,
)
    step_number = output_index - 1
    time_s = output_times[output_index]
    previous_time_s = output_index == 1 ? time_s : output_times[output_index - 1]
    T9 = profiles.T9(time_s)
    rho = profiles.rho(time_s)
    Y = interpolated_state(solver_times, history, time_s)
    epsilon = energy_generation_rate(network, Y, rho, T9; screening=screening_model)
    values = output_values(network, Y, X_file, records)
    path = joinpath(output_dir, @sprintf("iso_massf%05d.DAT", step_number))
    write_iso_massf_file(path, step_number, time_s, previous_time_s, T9, rho, epsilon, records, values)
end

function read_iso_massf_mass_fractions(path::AbstractString)
    X = Dict{String,Float64}()
    open(path, "r") do io
        for raw_line in eachline(io)
            fields = split(strip(raw_line))
            length(fields) >= 6 || continue
            all(isdigit, fields[1]) || continue

            label = join(fields[6:end], "")
            name = if label == "NEUT"
                "n"
            elseif label == "PROT"
                "p"
            else
                normalize_species_name(label)
            end
            X[name] = parse(Float64, fields[5])
        end
    end
    return X
end

function combined_output_mass_fractions(network::ReactionNetwork, Y, X_file::Dict{String,Float64})
    X = copy(X_file)
    active_X = mass_fractions_from_abundances(network, Y)
    for (name, value) in active_X
        X[name] = value
    end
    return X
end

function reaction_participation_counts(network::ReactionNetwork)
    production = Dict(name => 0 for name in network.species)
    destruction = Dict(name => 0 for name in network.species)
    for compiled in network.compiled_reactions
        for (index, delta) in pairs(compiled.stoichiometric_delta)
            delta > 0.0 && (production[network.species[index]] += 1)
            delta < 0.0 && (destruction[network.species[index]] += 1)
        end
    end
    return production, destruction
end

function print_ppn_report(
    network::ReactionNetwork,
    forward_tables,
    reverse_summary,
    X_file::Dict{String,Float64},
    history,
    mass_drift,
    trajectory::Trajectory,
    reference_path::AbstractString,
    output_dir::AbstractString,
)
    active_initial_mass = sum(get(X_file, name, 0.0) for name in network.species; init=0.0)
    final_X = combined_output_mass_fractions(network, history[end, :], X_file)
    source_counts = Dict{String,Int}()
    for reaction in network.reactions
        source_counts[reaction.rate_table.source] = get(source_counts, reaction.rate_table.source, 0) + 1
    end

    println()
    println("Network coverage report")
    println("forward H-Ca tables = ", length(forward_tables))
    println("network species = ", length(network.species))
    println("network reactions = ", length(network.reactions))
    println("reverse rates: explicit STARLIB = ", reverse_summary.explicit,
        ", generated detailed-balance = ", reverse_summary.generated,
        ", unavailable = ", reverse_summary.missing)
    println("generated reverse reaction tables in solve = ", get(source_counts, "detail_balance", 0))
    println("active initial mass = ", @sprintf("%.6E", active_initial_mass))
    println("inert carried initial mass = ", @sprintf("%.6E", 1.0 - active_initial_mass))
    println("trajectory T9 peak = ", @sprintf("%.5f", maximum(trajectory.T9)))
    println("mass-fraction drift = ", @sprintf("%.3E", mass_drift.drift))

    isfile(reference_path) || return NamedTuple[]

    reference_X = read_iso_massf_mass_fractions(reference_path)
    reference_names = sort(collect(keys(reference_X)); by=species_sort_key)
    active = count(name -> haskey(network.species_index, name), reference_names)
    inert = count(name -> !haskey(network.species_index, name) && haskey(X_file, name), reference_names)
    missing = length(reference_names) - active - inert
    println("Iliadis reference = ", reference_path)
    println("Iliadis species active/inert/missing = ", active, "/", inert, "/", missing)

    return write_and_print_iliadis_comparison(
        "Iliadis comparison (final network state)",
        joinpath(output_dir, "comparison_final_vs_iliadis.csv"),
        reference_X,
        final_X,
        network,
        X_file,
    )
end

function write_postdecay_iso_file(
    output_dir::AbstractString,
    final_step_number::Int,
    final_time_s::Float64,
    decay_time_s::Float64,
    final_T9::Float64,
    final_rho::Float64,
    records,
    decayed_X::Dict{String,Float64},
)
    values = output_values_from_mass_fractions(decayed_X, records)
    path = joinpath(output_dir, "iso_massfDECAY.DAT")
    write_iso_massf_file(
        path,
        final_step_number,
        final_time_s + decay_time_s,
        final_time_s,
        final_T9,
        final_rho,
        0.0,
        records,
        values,
    )
    return path
end

function stopped_output_times(trajectory_times::AbstractVector{<:Real}, stop_time_s::Real)
    stop_time = Float64(stop_time_s)
    times = Float64[t for t in trajectory_times if t <= stop_time]
    if isempty(times) || !isapprox(last(times), stop_time; rtol=0.0, atol=1.0e-9)
        push!(times, stop_time)
    end
    return times
end

function truncate_solution(times::Vector{Float64}, history::Matrix{Float64}, stop_time_s::Float64, stop_state::Vector{Float64})
    keep = findall(t -> t < stop_time_s, times)
    truncated_times = Float64[times[i] for i in keep]
    rows = Vector{Float64}[collect(history[i, :]) for i in keep]

    if isempty(truncated_times) || !isapprox(last(truncated_times), stop_time_s; rtol=0.0, atol=1.0e-9)
        push!(truncated_times, stop_time_s)
        push!(rows, copy(stop_state))
    else
        rows[end] = copy(stop_state)
    end

    truncated_history = Matrix{Float64}(undef, length(rows), size(history, 2))
    for (i, row) in pairs(rows)
        truncated_history[i, :] .= row
    end
    return truncated_times, truncated_history
end

function species_sort_key(name::AbstractString)
    info = species_from_name(name)
    return (info.Z, info.A, name)
end

function comparison_status(name::AbstractString, network::ReactionNetwork, X_file::Dict{String,Float64})
    haskey(network.species_index, name) && return "active"
    haskey(X_file, name) && return "inert"
    return "missing"
end

function iliadis_comparison_rows(
    reference_X::Dict{String,Float64},
    model_X::Dict{String,Float64},
    network::ReactionNetwork,
    X_file::Dict{String,Float64},
)
    names = sort(collect(keys(reference_X)); by=species_sort_key)
    rows = NamedTuple[]
    for name in names
        ref = reference_X[name]
        ours = get(model_X, name, 0.0)
        ratio = ref > 0.0 ? max(ours, MASS_FRACTION_FLOOR) / ref : (ours > 0.0 ? Inf : NaN)
        log10_ratio = ref > 0.0 ? log10(ratio) : NaN
        info = species_from_name(name)
        push!(
            rows,
            (
                name=name,
                label=iso_label(name),
                Z=info.Z,
                A=info.A,
                ref=ref,
                ours=ours,
                residual=ours - ref,
                abs_residual=abs(ours - ref),
                ratio=ratio,
                log10_ratio=log10_ratio,
                abs_log10_ratio=isfinite(log10_ratio) ? abs(log10_ratio) : Inf,
                status=comparison_status(name, network, X_file),
            ),
        )
    end
    return rows
end

function csv_float(value::Real)
    isnan(value) && return "NaN"
    isinf(value) && return value > 0.0 ? "Inf" : "-Inf"
    return @sprintf("%.16E", value)
end

function write_iliadis_comparison_csv(path::AbstractString, rows)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "isotope,Z,A,ref_mass_fraction,model_mass_fraction,residual,abs_residual,ratio,log10_ratio,abs_log10_ratio,status")
        for row in rows
            println(
                io,
                join(
                    (
                        strip(row.label),
                        string(row.Z),
                        string(row.A),
                        csv_float(row.ref),
                        csv_float(row.ours),
                        csv_float(row.residual),
                        csv_float(row.abs_residual),
                        csv_float(row.ratio),
                        csv_float(row.log10_ratio),
                        csv_float(row.abs_log10_ratio),
                        row.status,
                    ),
                    ",",
                ),
            )
        end
    end
    return path
end

function group_name(name::AbstractString)
    info = species_from_name(name)
    info.Z <= 1 && return "H"
    2 <= info.Z <= 5 && return "He-Li-Be-B"
    6 <= info.Z <= 8 && return "CNO"
    9 <= info.Z <= 11 && return "F-Ne-Na"
    12 <= info.Z <= 13 && return "Mg-Al"
    14 <= info.Z <= 20 && return "Si-Ca"
    return "Other"
end

function print_group_comparison(rows)
    groups = ["H", "He-Li-Be-B", "CNO", "F-Ne-Na", "Mg-Al", "Si-Ca", "Other"]
    println("group sums over Iliadis isotope set:")
    for group in groups
        selected = [row for row in rows if group_name(row.name) == group]
        isempty(selected) && continue
        ref = sum(row.ref for row in selected; init=0.0)
        ours = sum(row.ours for row in selected; init=0.0)
        ratio = ref > 0.0 ? ours / ref : Inf
        @printf("  %-10s ref=%10.5E ours=%10.5E delta=%+10.3E ratio=%10.3E\n", group, ref, ours, ours - ref, ratio)
    end
end

function print_iliadis_comparison(title::AbstractString, rows, csv_path::AbstractString; top_n::Int=12)
    isempty(rows) && return
    println(title)
    println("comparison csv = ", csv_path)
    println("Iliadis mass sum = ", @sprintf("%.12f", sum(row.ref for row in rows; init=0.0)))
    println("model mass over Iliadis isotopes = ", @sprintf("%.12f", sum(row.ours for row in rows; init=0.0)))
    print_group_comparison(rows)

    absolute_rows = sort(collect(rows); by=row -> row.abs_residual, rev=true)
    println("largest absolute residuals:")
    for row in Iterators.take(absolute_rows, top_n)
        @printf(
            "  %-6s ref=%10.3E ours=%10.3E delta=%+10.3E ratio=%10.3E %-7s\n",
            row.label,
            row.ref,
            row.ours,
            row.residual,
            row.ratio,
            row.status,
        )
    end

    log_rows = sort([row for row in rows if row.ref > 0.0]; by=row -> row.abs_log10_ratio, rev=true)
    println("largest log-ratio residuals:")
    for row in Iterators.take(log_rows, top_n)
        @printf(
            "  %-6s ref=%10.3E ours=%10.3E ratio=%10.3E log10=%+7.3f %-7s\n",
            row.label,
            row.ref,
            row.ours,
            row.ratio,
            row.log10_ratio,
            row.status,
        )
    end
end

function write_and_print_iliadis_comparison(
    title::AbstractString,
    csv_path::AbstractString,
    reference_X::Dict{String,Float64},
    model_X::Dict{String,Float64},
    network::ReactionNetwork,
    X_file::Dict{String,Float64},
)
    rows = iliadis_comparison_rows(reference_X, model_X, network, X_file)
    written_path = write_iliadis_comparison_csv(csv_path, rows)
    print_iliadis_comparison(title, rows, written_path)
    return rows
end

function species_flux_contributions(network::ReactionNetwork, integrated_fluxes_by_reaction, name::AbstractString)
    haskey(network.species_index, name) || return NamedTuple[]
    index = network.species_index[name]
    A = network.species_info[index].A
    rows = NamedTuple[]
    for (reaction_index, compiled) in pairs(network.compiled_reactions)
        delta = compiled.stoichiometric_delta[index]
        delta == 0.0 && continue
        contribution = A * delta * integrated_fluxes_by_reaction[reaction_index]
        push!(
            rows,
            (
                abs_contribution=abs(contribution),
                contribution=contribution,
                reaction=reaction_string(network.reactions[reaction_index]),
                source=network.reactions[reaction_index].rate_table.source,
            ),
        )
    end
    sort!(rows; by=row -> row.abs_contribution, rev=true)
    return rows
end

function print_top_flux_rows(rows, positive::Bool, top_n::Int)
    printed = 0
    for row in rows
        positive ? row.contribution > 0.0 || continue : row.contribution < 0.0 || continue
        @printf("    %+11.4E  %-32s %-8s\n", row.contribution, row.reaction, row.source)
        printed += 1
        printed >= top_n && break
    end
    printed == 0 && println("    none")
end

function print_integrated_flux_report(
    network::ReactionNetwork,
    times,
    history,
    profiles,
    screening_model,
    comparison_rows;
    species_count::Int=6,
    reactions_per_species::Int=5,
)
    species_count > 0 || return
    active_rows = [row for row in comparison_rows if haskey(network.species_index, row.name)]
    sort!(active_rows; by=row -> row.abs_residual, rev=true)
    isempty(active_rows) && return

    target_rows = active_rows[1:min(species_count, length(active_rows))]
    flux_history = reaction_flux_history(network, history, times, profiles.rho, profiles.T9; screening=screening_model)
    integrated = integrated_fluxes(times, flux_history)

    println()
    println("Integrated flux report for worst residual isotopes")
    println("contributions are approximate integrated mass-fraction changes over the network solve")
    for row in target_rows
        println("  ", row.label, " ref=", @sprintf("%.3E", row.ref),
            " model=", @sprintf("%.3E", row.ours),
            " residual=", @sprintf("%+.3E", row.residual))
        flux_rows = species_flux_contributions(network, integrated, row.name)
        println("  top production:")
        print_top_flux_rows(flux_rows, true, reactions_per_species)
        println("  top destruction:")
        print_top_flux_rows(flux_rows, false, reactions_per_species)
    end
end

project_root = dirname(@__DIR__)
trajectory_path = required_file([joinpath(project_root, "trajectory.input")], "trajectory file")
abundance_path = required_file(
    [joinpath(project_root, "initial_abundance.DAT"), joinpath(project_root, "initial_abundance.dat")],
    "initial abundance file",
)
template_candidates = [joinpath(project_root, "iso_massf00804.DAT")]
template_index = findfirst(isfile, template_candidates)
template_path = template_index === nothing ? nothing : template_candidates[template_index]
output_dir = joinpath(project_root, "outputs", "single_zone_nova_ppn")
output_stride = CLI.output_stride

tables = read_starlib()
trajectory = read_trajectory(trajectory_path)
profiles = trajectory_profiles(trajectory)
raw_X_file = read_initial_abundances(abundance_path)
X_file = read_initial_abundances(abundance_path; normalize=true)
forward_tables = select_h_ca_reaction_tables(tables, keys(raw_X_file))
reverse_summary = add_reverse_reaction_tables(tables, forward_tables; generate_detailed_balance=true)
network = network_from_tables(reverse_summary.tables)
validation = network_validation_report(network; throw_on_error=true)

trajectory_start_time = first(trajectory.time)
trajectory_end_time = last(trajectory.time)
integration_end_time = trajectory_end_time
stop_reasons = String[]

if CLI.stop_time_s !== nothing
    trajectory_start_time < CLI.stop_time_s <= trajectory_end_time ||
        throw(ArgumentError("--stop-time must be inside the trajectory interval ($(trajectory_start_time), $(trajectory_end_time)] s"))
    integration_end_time = min(integration_end_time, CLI.stop_time_s)
    push!(stop_reasons, @sprintf("time %.6E s", CLI.stop_time_s))
end

if CLI.stop_temperature_T9 !== nothing
    temperature_stop_time = first_cooling_threshold_time(trajectory, CLI.stop_temperature_T9)
    temperature_stop_time === nothing &&
        throw(ArgumentError("trajectory never cools to T9=$(CLI.stop_temperature_T9) after peak"))
    integration_end_time = min(integration_end_time, temperature_stop_time)
    push!(stop_reasons, @sprintf("cooling T9 %.5f at %.6E s", CLI.stop_temperature_T9, temperature_stop_time))
end

trajectory_duration = integration_end_time - trajectory_start_time
dt_initial = trajectory_duration > 100.0 ? 1.0 : 0.02
dt_min = trajectory_duration > 100.0 ? 1.0e-8 : 1.0e-10
dt_max = trajectory_duration > 100.0 ? 20.0 : 0.05
screening_model = CLI.screening

X0 = Dict(name => value for (name, value) in X_file if haskey(network.species_index, name))
Y0 = abundances_from_mass_fractions(network, X0)

times, history, solver_stats = solve_network_adaptive(
    network,
    Y0,
    (trajectory_start_time, integration_end_time),
    dt_initial,
    profiles.rho,
    profiles.T9;
    method=:backward_euler,
    max_fractional_change=0.50,
    max_absolute_change=1.0e-4,
    abundance_floor=1.0e-8,
    dt_min=dt_min,
    dt_max=dt_max,
    screening=screening_model,
    return_stats=true,
    max_newton_iterations=80,
)

hydrogen_stop_applied = false
if CLI.stop_hydrogen_X !== nothing
    crossing = first_mass_fraction_threshold_crossing(network, times, history, "p", CLI.stop_hydrogen_X; direction=:down)
    if crossing === nothing
        push!(stop_reasons, @sprintf("hydrogen X(PROT) %.6E not reached", CLI.stop_hydrogen_X))
    elseif crossing.time < last(times) - 1.0e-9
        times, history = truncate_solution(times, history, crossing.time, crossing.state)
        hydrogen_stop_applied = true
        push!(stop_reasons, @sprintf("hydrogen X(PROT) %.6E at %.6E s", CLI.stop_hydrogen_X, crossing.time))
    else
        push!(stop_reasons, @sprintf("hydrogen X(PROT) %.6E at final state", CLI.stop_hydrogen_X))
    end
end

energy_history = energy_generation_history(network, history, times, profiles.rho, profiles.T9; screening=screening_model)
mass_drift = mass_fraction_drift(network, history)
records = build_output_records(template_path, network, X_file)
output_times = stopped_output_times(trajectory.time, last(times))
indices = saved_indices(length(output_times), output_stride)
clean_previous_outputs(output_dir)
mkpath(output_dir)

if CLI.jobs > 1 && Base.Threads.nthreads() > 1
    Base.Threads.@threads for i in eachindex(indices)
        write_saved_iso_state(indices[i], output_dir, output_times, times, history, profiles, network, X_file, records, screening_model)
    end
else
    for output_index in indices
        write_saved_iso_state(output_index, output_dir, output_times, times, history, profiles, network, X_file, records, screening_model)
    end
end

csv_path = joinpath(output_dir, "mass_fractions.csv")
write_csv(csv_path, indices, output_times, times, history, profiles, network, X_file, records, screening_model)

println("Single-zone nova PPN")
println("trajectory = ", trajectory_path)
println("initial abundances = ", abundance_path)
println("template = ", template_path === nothing ? "generated from species inventory" : template_path)
println("output directory = ", output_dir)
println("csv = ", csv_path)
println("validated reactions = ", validation.num_reactions)
println("active species = ", length(network.species))
println("output isotopes = ", length(records))
println("input abundance raw total = ", sum(values(raw_X_file); init=0.0))
println("active network initial mass = ", sum(values(X0); init=0.0))
println("inert/outside-network initial mass = ", sum(values(X_file); init=0.0) - sum(values(X0); init=0.0))
println("screening = ", screening_label(screening_model))
println("post-decay time = ", CLI.decay_time_s, " s")
println("stop controls = ", isempty(stop_reasons) ? "none" : join(stop_reasons, "; "))
hydrogen_stop_applied && println("hydrogen stop truncated output/report state after integration; solver step stats cover the underlying solve interval")
println("jobs requested = ", CLI.jobs)
println("Julia threads = ", Base.Threads.nthreads())
println("accepted timesteps = ", solver_stats.accepted_steps)
println("rejected timesteps = ", solver_stats.rejected_steps)
println("dt range = ", solver_stats.min_dt, " to ", solver_stats.max_dt, " s")
println("Newton iterations mean/max = ", solver_stats.mean_newton_iterations, " / ", solver_stats.max_newton_iterations)
println("time = ", first(times), " to ", last(times), " s")
println("final T9/rho = ", profiles.T9(last(times)), " / ", profiles.rho(last(times)))
println("total active mass fraction = ", mass_drift.initial, " to ", mass_drift.final)
println("peak epsilon_nuc = ", maximum(energy_history), " erg g^-1 s^-1")
println("solver states = ", length(times))
println("written trajectory states = ", length(indices), " of ", length(output_times), " (stride=", output_stride, ")")

iliadis_reference_path = joinpath(project_root, "outputs", "iliadis2002_jch1", "iso_massf00000.DAT")
comparison_rows_for_flux = print_ppn_report(network, forward_tables, reverse_summary, X_file, history, mass_drift, trajectory, iliadis_reference_path, output_dir)

if CLI.decay_time_s > 0.0
    final_X = combined_output_mass_fractions(network, history[end, :], X_file)
    final_T9 = profiles.T9(last(times))
    final_rho = profiles.rho(last(times))
    decay_result = decay_mass_fractions(tables, final_X, CLI.decay_time_s; T9=final_T9)
    decay_path = write_postdecay_iso_file(
        output_dir,
        length(output_times) - 1,
        last(output_times),
        CLI.decay_time_s,
        final_T9,
        final_rho,
        records,
        decay_result.mass_fractions,
    )
    println("post-decay weak reactions = ", length(decay_result.decay_tables))
    println("post-decay output = ", decay_path)

    if isfile(iliadis_reference_path)
        reference_X = read_iso_massf_mass_fractions(iliadis_reference_path)
        comparison_rows_for_flux = write_and_print_iliadis_comparison(
            "Iliadis comparison (post-decay state)",
            joinpath(output_dir, "comparison_postdecay_vs_iliadis.csv"),
            reference_X,
            decay_result.mass_fractions,
            network,
            X_file,
        )
    end
end

if CLI.flux_report_count > 0 && !isempty(comparison_rows_for_flux)
    print_integrated_flux_report(
        network,
        times,
        history,
        profiles,
        screening_model,
        comparison_rows_for_flux;
        species_count=CLI.flux_report_count,
    )
end
