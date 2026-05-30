using Printf
using ReacNetJl

#=
    Full single-zone nova PPN example

This script evolves the current expanded nova network over the project-root
`trajectory.input` using `initial_abundance.DAT` or `initial_abundance.dat`.
It writes one `iso_massfXXXXX.DAT` mass-fraction file per saved trajectory
state, plus a wide CSV containing every isotope column. With the current
805-row trajectory, the last default file is `iso_massf00804.DAT`.

Set `NOVA_PPN_OUTPUT_STRIDE=1000` for quick checks. The default stride is 1,
which writes every accepted solver state.
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
        if occursin(r"^iso_massf\d+\.DAT$", base) || base == "mass_fractions.csv"
            rm(path; force=true)
        end
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
output_stride = parse(Int, get(ENV, "NOVA_PPN_OUTPUT_STRIDE", "1"))

tables = read_starlib()
network = network_from_labels(tables, labels)
validation = network_validation_report(network; throw_on_error=true)
trajectory = read_trajectory(trajectory_path)
profiles = trajectory_profiles(trajectory)

trajectory_duration = last(trajectory.time) - first(trajectory.time)
dt_initial = trajectory_duration > 100.0 ? 1.0 : 0.02
dt_min = trajectory_duration > 100.0 ? 1.0e-4 : 1.0e-6
dt_max = trajectory_duration > 100.0 ? 20.0 : 0.05
screening_model = :weak

raw_X_file = read_initial_abundances(abundance_path)
X_file = read_initial_abundances(abundance_path; normalize=true)
X0 = Dict(name => value for (name, value) in X_file if haskey(network.species_index, name))
Y0 = abundances_from_mass_fractions(network, X0)

times, history, solver_stats = solve_network_adaptive(
    network,
    Y0,
    (first(trajectory.time), last(trajectory.time)),
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
)

energy_history = energy_generation_history(network, history, times, profiles.rho, profiles.T9; screening=screening_model)
mass_drift = mass_fraction_drift(network, history)
records = build_output_records(template_path, network, X_file)
output_times = trajectory.time
indices = saved_indices(length(output_times), output_stride)
clean_previous_outputs(output_dir)

for output_index in indices
    step_number = output_index - 1
    time_s = output_times[output_index]
    previous_time_s = output_index == 1 ? time_s : output_times[output_index - 1]
    T9 = profiles.T9(time_s)
    rho = profiles.rho(time_s)
    Y = interpolated_state(times, history, time_s)
    epsilon = energy_generation_rate(network, Y, rho, T9; screening=screening_model)
    values = output_values(network, Y, X_file, records)
    path = joinpath(output_dir, @sprintf("iso_massf%05d.DAT", step_number))
    write_iso_massf_file(
        path,
        step_number,
        time_s,
        previous_time_s,
        T9,
        rho,
        epsilon,
        records,
        values,
    )
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
println("screening = ", screening_model)
println("accepted timesteps = ", solver_stats.accepted_steps)
println("rejected timesteps = ", solver_stats.rejected_steps)
println("dt range = ", solver_stats.min_dt, " to ", solver_stats.max_dt, " s")
println("Newton iterations mean/max = ", solver_stats.mean_newton_iterations, " / ", solver_stats.max_newton_iterations)
println("time = ", first(times), " to ", last(times), " s")
println("total active mass fraction = ", mass_drift.initial, " to ", mass_drift.final)
println("peak epsilon_nuc = ", maximum(energy_history), " erg g^-1 s^-1")
println("solver states = ", length(times))
println("written trajectory states = ", length(indices), " of ", length(output_times), " (stride=", output_stride, ")")
