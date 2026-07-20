# Trajectory and initial-abundance file readers.

"""
    Trajectory

A single-zone thermodynamic trajectory with time in seconds, temperature `T9` in
GK, and density `rho` in g cm^-3.
"""
struct Trajectory
    time::Vector{Float64}
    T9::Vector{Float64}
    rho::Vector{Float64}
end

# Sanity-check a parsed trajectory: strictly increasing time (required for
# `_linear_profile`'s interpolation and for the solver's forward-only time
# stepping) and strictly positive T9/rho (both appear in rate/flux
# expressions as bases of fractional powers and inside logarithms, so a
# non-positive value would silently produce NaN/Inf deep in the solver
# instead of failing at the input-parsing stage where the problem is obvious).
function _validate_trajectory(time::Vector{Float64}, T9::Vector{Float64}, rho::Vector{Float64})
    length(time) == length(T9) == length(rho) || throw(ArgumentError("trajectory columns must have the same length"))
    length(time) >= 2 || throw(ArgumentError("trajectory must contain at least two rows"))

    for i in 1:(length(time)-1)
        time[i+1] > time[i] || throw(ArgumentError("trajectory time values must be strictly increasing"))
    end

    any(<=(0.0), T9) && throw(ArgumentError("trajectory T9 values must be positive"))
    any(<=(0.0), rho) && throw(ArgumentError("trajectory rho values must be positive"))
end

"""
    read_trajectory(path)

Read a whitespace- or comma-separated trajectory file with columns:

    time_s  T9  rho

Lines beginning with `#` and blank lines are ignored. Metadata assignments,
one per line before the numeric rows, select the units the numeric columns
are given in (all converted to seconds/GK/g cm^-3 on return, regardless of
which units the file uses):

    AGEUNIT = SEC | YRS
    TUNIT   = T9K | T8K
    RHOUNIT = CGS | LOG

`AGEUNIT=YRS` converts using a 365.25-day Julian year; `RHOUNIT=LOG` treats
the column as `log10(rho)`. Defaults to `SEC`/`T9K`/`CGS` for any unit not
explicitly set. Validated with `_validate_trajectory` before being wrapped in
a `Trajectory`.
"""
function read_trajectory(path::AbstractString)
    time = Float64[]
    T9 = Float64[]
    rho = Float64[]
    age_unit = "SEC"
    temperature_unit = "T9K"
    density_unit = "CGS"

    open(path, "r") do io
        for raw_line in eachline(io)
            line = strip(split(raw_line, '#'; limit=2)[1])
            isempty(line) && continue

            if occursin("=", line)
                key_value = split(line, '='; limit=2)
                key = uppercase(strip(key_value[1]))
                value = uppercase(strip(key_value[2]))
                if key == "AGEUNIT"
                    age_unit = value
                elseif key == "TUNIT"
                    temperature_unit = value
                elseif key == "RHOUNIT"
                    density_unit = value
                end
                continue
            end

            fields = split(replace(line, ',' => ' '))
            length(fields) >= 3 || throw(ArgumentError("trajectory row must contain at least three columns: $raw_line"))

            time_value = parse(Float64, fields[1])
            T_value = parse(Float64, fields[2])
            rho_value = parse(Float64, fields[3])

            if age_unit == "SEC"
                push!(time, time_value)
            elseif age_unit == "YRS"
                push!(time, time_value * 365.25 * 24.0 * 60.0 * 60.0)
            else
                throw(ArgumentError("unsupported trajectory AGEUNIT=$age_unit; use SEC or YRS"))
            end

            if temperature_unit == "T9K"
                push!(T9, T_value)
            elseif temperature_unit == "T8K"
                push!(T9, T_value / 10.0)
            else
                throw(ArgumentError("unsupported trajectory TUNIT=$temperature_unit; use T9K or T8K"))
            end

            if density_unit == "CGS"
                push!(rho, rho_value)
            elseif density_unit == "LOG"
                push!(rho, 10.0^rho_value)
            else
                throw(ArgumentError("unsupported trajectory RHOUNIT=$density_unit; use CGS or LOG"))
            end
        end
    end

    _validate_trajectory(time, T9, rho)
    return Trajectory(time, T9, rho)
end

# Parse one initial-abundance row's species name and mass fraction, in any
# of the three accepted layouts (see `read_initial_abundances`): the special
# "1 PROT X" proton row (3 fields), "Z sym A X" split across separate symbol
# and mass-number fields (4 fields), or "Z symA X" with them already joined
# (3 fields, e.g. "26 fe56 6.0e-4").
function _parse_initial_abundance_species(fields::Vector{SubString{String}}, raw_line::AbstractString)
    if length(fields) == 3 && uppercase(fields[2]) == "PROT"
        return "p", parse(Float64, fields[3])
    elseif length(fields) == 4
        return normalize_species_name(string(fields[2], fields[3])), parse(Float64, fields[4])
    elseif length(fields) == 3
        return normalize_species_name(fields[2]), parse(Float64, fields[3])
    end

    throw(ArgumentError("initial abundance row must look like `Z sym A X`, `Z symA X`, or `1 PROT X`: $raw_line"))
end

"""
    read_initial_abundances(path; normalize=false)

Read an initial-abundance mass-fraction file with rows such as:

    1 PROT  3.5e-1
    6 c 12  1.0e-2
    26 fe56 6.0e-4

Returns a dictionary keyed by normalized species name (see
`_parse_initial_abundance_species` for the accepted row layouts). Duplicate
entries for the same species raise an error rather than silently overwriting.
If `normalize=true`, all mass fractions are divided by the file total, so
`\\sum_i X_i = 1` exactly -- useful when the raw file is a rounded/truncated
composition (as digitized tables often are) that doesn't sum to precisely 1.
"""
function read_initial_abundances(path::AbstractString; normalize::Bool=false)
    X = Dict{String,Float64}()

    open(path, "r") do io
        for raw_line in eachline(io)
            line = strip(split(raw_line, '#'; limit=2)[1])
            isempty(line) && continue
            fields = split(line)
            name, value = _parse_initial_abundance_species(fields, raw_line)
            haskey(X, name) && throw(ArgumentError("duplicate initial abundance entry for species '$name'"))
            X[name] = value
        end
    end

    normalize || return X
    total = sum(values(X); init=0.0)
    total > 0.0 || throw(ArgumentError("cannot normalize initial abundances with non-positive total"))
    return Dict(name => value / total for (name, value) in X)
end

# Linear (not log-log, unlike rate interpolation) interpolation of a
# trajectory column at an arbitrary time -- T9(t) and rho(t) are smooth
# physical profiles sampled at whatever resolution the trajectory file used,
# not quantities spanning orders of magnitude between adjacent grid points
# the way rates do, so ordinary linear interpolation is appropriate. Throws
# if `value` (time) falls outside the trajectory's tabulated range, since the
# solver should never need to evaluate T9/rho beyond the prescribed history.
function _linear_profile(x::Vector{Float64}, y::Vector{Float64}, value::Real)
    v = Float64(value)
    first(x) <= v <= last(x) || throw(ArgumentError("profile value $v is outside the trajectory range $(first(x))–$(last(x))"))

    i = searchsortedlast(x, v)
    if i == length(x) || x[i] == v
        return y[i]
    end

    weight = (v - x[i]) / (x[i+1] - x[i])
    return (1.0 - weight) * y[i] + weight * y[i+1]
end

"""
    trajectory_profiles(trajectory)

Return callable density and temperature profiles as a named tuple, closing
over the trajectory data so callers don't need to pass it explicitly at
every solver step:

    profiles = trajectory_profiles(traj)
    rho_t = profiles.rho(t)
    T9_t = profiles.T9(t)

Both close over `_linear_profile`, so evaluating at a time outside the
trajectory's range throws rather than silently extrapolating.
"""
function trajectory_profiles(trajectory::Trajectory)
    rho_profile(t) = _linear_profile(trajectory.time, trajectory.rho, t)
    T9_profile(t) = _linear_profile(trajectory.time, trajectory.T9, t)
    return (rho=rho_profile, T9=T9_profile)
end

"""
    first_cooling_threshold_time(trajectory, threshold_T9)

Find the first time, on the post-peak cooling branch, that `T9` falls to
`threshold_T9`, linearly interpolating between bracketing trajectory rows
(same interpolation approach as `first_mass_fraction_threshold_crossing`).
Searches only from the trajectory's temperature maximum onward, so a
temperature history with an early warm-up phase that happens to pass through
`threshold_T9` on the way *up* to peak isn't mistaken for the cooling
crossing -- this is the basis for `--stop-temperature` in the example
driver, used to stop integration once burning has clearly ended. Returns
`nothing` if the trajectory never cools to `threshold_T9` after its peak.
"""
function first_cooling_threshold_time(trajectory::Trajectory, threshold_T9::Real)
    threshold = Float64(threshold_T9)
    threshold > 0.0 || throw(ArgumentError("threshold_T9 must be positive"))

    peak_index = argmax(trajectory.T9)
    for i in peak_index:(length(trajectory.time)-1)
        T0 = trajectory.T9[i]
        T1 = trajectory.T9[i+1]
        t0 = trajectory.time[i]
        t1 = trajectory.time[i+1]

        T0 == threshold && return t0
        if T0 > threshold && T1 <= threshold
            fraction = (threshold - T0) / (T1 - T0)
            return t0 + fraction * (t1 - t0)
        end
    end

    return nothing
end
