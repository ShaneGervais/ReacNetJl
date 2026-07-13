# Reader for rate tables extracted from the original NACRE/Iliadis 2001 papers.

"""
    read_tabulated_rates(path; source, partition_functions=nothing,
                         include_total_variants=false)

Read a ReacNetJl tabulated-rate file (as produced by the paper-table
extraction scripts) into `ReactionRateTable` entries. The format is blocks of

    reaction: p ne20 -> na21 ; label: 20Ne(p,g) ; table: 3 ; variant: standard
    <T9> <rate>
    ...
    end

Q-values are computed from winvne mass excesses when `partition_functions`
is given (0 otherwise, which disables generated detailed-balance reverses
for those tables). Rates with `variant: total` duplicate the sum of their
ground/isomer siblings and are skipped unless `include_total_variants=true`.

Rows may carry two columns (`T9 rate`) or four (`T9 rate lower upper`); the
limits become STARLIB-style factor uncertainties `sqrt(upper/lower)`.
Threshold reactions are tabulated only above their onset temperature, where
the papers' omitted rows mean "negligible": each table is padded down to the
STARLIB grid start with the rate floor (and flat up to 10 GK if needed) so
trajectory evaluation never leaves the interpolation domain.
"""
function read_tabulated_rates(
    path::AbstractString;
    source::AbstractString,
    partition_functions=nothing,
    mass_excesses=nothing,
    include_total_variants::Bool=false,
)
    isfile(path) || error("tabulated rate file not found at $path")
    tables = ReactionRateTable[]

    reactants = String[]
    products = String[]
    variant = "standard"
    T9 = Float64[]
    rate = Float64[]
    factor_uncertainty = Float64[]

    header_pattern = r"^reaction:\s+(.+?)\s+->\s+(.+?)\s+;\s+label:.*;\s+variant:\s+(\w+)\s*$"

    for (line_number, raw_line) in enumerate(eachline(path))
        line = strip(raw_line)
        (isempty(line) || startswith(line, "#")) && continue

        if startswith(line, "reaction:")
            m = match(header_pattern, line)
            m === nothing && error("malformed reaction header at line $line_number of $path: $line")
            reactants = normalize_species_name.(split(m.captures[1]))
            products = normalize_species_name.(split(m.captures[2]))
            variant = m.captures[3]
            T9 = Float64[]
            rate = Float64[]
            factor_uncertainty = Float64[]
        elseif line == "end"
            isempty(reactants) && error("'end' without a reaction header at line $line_number of $path")
            if !(variant == "total" && !include_total_variants) && !isempty(T9)
                q = 0.0
                for mass_source in (mass_excesses, partition_functions)
                    mass_source === nothing && continue
                    computed = reaction_q_value(mass_source, reactants, products)
                    computed === nothing && continue
                    q = computed
                    break
                end
                if first(T9) > first(STARLIB_T9_GRID)
                    pushfirst!(T9, first(STARLIB_T9_GRID))
                    pushfirst!(rate, GENERATED_REVERSE_RATE_FLOOR)
                    pushfirst!(factor_uncertainty, 1.0)
                end
                if last(T9) < last(STARLIB_T9_GRID)
                    push!(T9, last(STARLIB_T9_GRID))
                    push!(rate, last(rate))
                    push!(factor_uncertainty, last(factor_uncertainty))
                end
                chapter = length(reactants) == 2 ? (length(products) == 1 ? 4 : 5) : 1
                push!(tables, ReactionRateTable(
                    chapter,
                    copy(reactants),
                    copy(products),
                    String(source),
                    q,
                    copy(T9),
                    copy(rate),
                    copy(factor_uncertainty),
                ))
            end
            reactants = String[]
            products = String[]
        else
            fields = split(line)
            length(fields) in (2, 4) || error("malformed data row at line $line_number of $path: $line")
            push!(T9, parse(Float64, fields[1]))
            push!(rate, parse(Float64, fields[2]))
            if length(fields) == 4
                lower = parse(Float64, fields[3])
                upper = parse(Float64, fields[4])
                (lower > 0.0 && upper >= lower) || error("invalid rate limits at line $line_number of $path: $line")
                push!(factor_uncertainty, sqrt(upper / lower))
            else
                push!(factor_uncertainty, 1.0)
            end
        end
    end

    return tables
end

"""
    read_iliadis2001_rates(path=DEFAULT_ILIADIS2001_PATH; kwargs...)

Read the recommended reaction rates of Iliadis et al. (2001, ApJS 134, 151),
extracted from the paper's Tables 3-9, as `ReactionRateTable`s with source
`"il01tab"`. See `read_tabulated_rates` for the keyword arguments.
"""
read_iliadis2001_rates(path::AbstractString=DEFAULT_ILIADIS2001_PATH; kwargs...) =
    read_tabulated_rates(path; source="il01tab", kwargs...)

"""
    read_nacre_rates(path=DEFAULT_NACRE_PATH; kwargs...)

Read the adopted reaction rates of the NACRE compilation (Angulo et al.
1999), extracted from the paper's rate tables, as `ReactionRateTable`s with
source `"nacrtab"`. See `read_tabulated_rates` for the keyword arguments.
"""
read_nacre_rates(path::AbstractString=DEFAULT_NACRE_PATH; kwargs...) =
    read_tabulated_rates(path; source="nacrtab", kwargs...)

