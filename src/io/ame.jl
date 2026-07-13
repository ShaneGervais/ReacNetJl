# AME2020 atomic mass evaluation reader.

"""
    read_ame_masses(path=DEFAULT_AME_PATH; include_estimated=true)

Read atomic mass excesses (MeV) from an Atomic Mass Evaluation file
(`mass.mas20.txt` from <https://www-nds.iaea.org/amdc/>), keyed by normalized
species name. AME marks extrapolated (non-experimental) values with `#`;
these are kept unless `include_estimated=false`. AME is the authoritative
source for Q-values; winvne mass excesses derive from older evaluations.
"""
function read_ame_masses(path::AbstractString=DEFAULT_AME_PATH; include_estimated::Bool=true)
    isfile(path) || error("AME mass table not found at $path; run ReacNetJl.fetch_data!() or data/download_rates.sh")
    masses = Dict{String,Float64}()

    for raw_line in eachline(path)
        line = rpad(raw_line, 60)
        Z = tryparse(Int, strip(line[10:14]))
        A = tryparse(Int, strip(line[15:19]))
        (Z === nothing || A === nothing) && continue
        element = strip(line[21:23])
        isempty(element) && continue

        field = strip(line[29:42])
        estimated = occursin('#', field)
        estimated && !include_estimated && continue
        excess_keV = tryparse(Float64, replace(field, '#' => ' ') |> strip)
        excess_keV === nothing && continue

        name = Z == 0 && A == 1 ? "n" : normalize_species_name(lowercase(element) * string(A))
        masses[name] = excess_keV / 1000.0
    end

    isempty(masses) && error("no mass excesses parsed from $path; is this an AME mass table?")
    return masses
end

