# winvne nuclear data: spins, partition functions, mass excesses.

"""
    PartitionFunctionTable

Ground-state statistical weights `g0 = 2J+1`, temperature-dependent
normalized partition functions `G(T9)`, and mass excesses (MeV) from a JINA
`winvne` file, keyed by normalized species name. `G` is tabulated on the
standard 24-point winvn temperature grid and is defined so that the total
partition function is `g0 * G(T9)`.
"""
struct PartitionFunctionTable
    T9::Vector{Float64}
    g0::Dict{String,Float64}
    G::Dict{String,Vector{Float64}}
    mass_excess::Dict{String,Float64}
end

const _WINVNE_T9_GRID = [
    0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5,
    2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
]

"""
    read_winvne(path=DEFAULT_WINVNE_PATH)

Read a JINA `winvne` nuclear-data file into a `PartitionFunctionTable`.

The file layout is: a title line, a packed temperature-grid line, a species
directory, then one record per nuclide consisting of a header line
`name A Z N spin mass_excess label` followed by 24 partition-function values.
Records are recognized by their header shape, so the directory needs no
special handling. If a normalized name appears more than once (for example
`al-6` normalizing onto `al26`), the first record wins.
"""
function read_winvne(path::AbstractString=DEFAULT_WINVNE_PATH)
    isfile(path) || error("winvne file not found at $path; run data/download_rates.sh to fetch it")

    g0 = Dict{String,Float64}()
    G = Dict{String,Vector{Float64}}()
    mass_excess = Dict{String,Float64}()

    open(path, "r") do io
        line_number = 0
        while !eof(io)
            line = readline(io)
            line_number += 1
            fields = split(line)
            length(fields) >= 6 || continue

            spin = tryparse(Float64, fields[5])
            mass_number = tryparse(Float64, fields[2])
            excess = tryparse(Float64, fields[6])
            (spin === nothing || mass_number === nothing) && continue

            name = normalize_species_name(fields[1])
            values = Float64[]
            sizehint!(values, length(_WINVNE_T9_GRID))
            while length(values) < length(_WINVNE_T9_GRID)
                eof(io) && error("unexpected end of winvne file inside the record of '$name' at line $line_number of $path")
                row = split(readline(io))
                line_number += 1
                for token in row
                    push!(values, parse(Float64, token))
                end
            end
            length(values) == length(_WINVNE_T9_GRID) ||
                error("winvne record of '$name' at line $line_number of $path has $(length(values)) partition-function values, expected $(length(_WINVNE_T9_GRID))")

            if !haskey(g0, name)
                g0[name] = 2.0 * spin + 1.0
                G[name] = values
                excess === nothing || (mass_excess[name] = excess)
            end
        end
    end

    return PartitionFunctionTable(copy(_WINVNE_T9_GRID), g0, G, mass_excess)
end

"""
    reaction_q_value(mass_excess, reactants, products)

Q-value in MeV from a mass-excess table (a `Dict` of MeV values keyed by
normalized species name, or a `PartitionFunctionTable`):
`Q = sum(ME reactants) - sum(ME products)`. Returns `nothing` when any
participant is missing from the table.
"""
function reaction_q_value(mass_excess::AbstractDict, reactants, products)
    total = 0.0
    for name in reactants
        excess = get(mass_excess, normalize_species_name(name), nothing)
        excess === nothing && return nothing
        total += excess
    end
    for name in products
        excess = get(mass_excess, normalize_species_name(name), nothing)
        excess === nothing && return nothing
        total -= excess
    end
    return total
end

reaction_q_value(pf::PartitionFunctionTable, reactants, products) =
    reaction_q_value(pf.mass_excess, reactants, products)

#=
Interpolate the normalized partition function of one species. Outside the
tabulated range the edge values are used; at nova temperatures (below the
0.1 GK grid start) G is 1 to excellent accuracy for all relevant nuclides.
Returns `nothing` when the species is not in the table.
=#
function _partition_function_at(pf::PartitionFunctionTable, name::AbstractString, T9::Real)
    key = normalize_species_name(name)
    haskey(pf.G, key) || return nothing
    values = pf.G[key]
    T = Float64(T9)
    T <= first(pf.T9) && return values[1]
    T >= last(pf.T9) && return values[end]
    return _interpolate_loglog(pf.T9, values, T; value_name="partition function")
end

_statistical_weight(pf::PartitionFunctionTable, name::AbstractString) =
    get(pf.g0, normalize_species_name(name), nothing)

