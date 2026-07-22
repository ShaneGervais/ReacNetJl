# STARLIB v6 tabulated rate reader (ReactionRateTable and its chapter layout rules).

"""
    ReactionRateTable

A temperature-dependent thermonuclear reaction rate from STARLIB.

`T9` is the temperature grid in GK. `rate` is the recommended STARLIB rate.
`factor_uncertainty` is the STARLIB multiplicative factor uncertainty.
"""
struct ReactionRateTable
    chapter::Int
    reactants::Vector{String}
    products::Vector{String}
    source::String
    q_value::Float64
    T9::Vector{Float64}
    rate::Vector{Float64}
    factor_uncertainty::Vector{Float64}
end


function _split_starlib_species(chapter::Int, species::Vector{String})
    arities = Dict(
        1 => (1, 1),
        2 => (1, 2),
        3 => (1, 3),
        4 => (2, 1),
        5 => (2, 2),
        6 => (2, 3),
        7 => (2, 4),
        8 => (3, 1),
        9 => (3, 2),
        10 => (4, 2),
    )

    if haskey(arities, chapter)
        nreactants, nproducts = arities[chapter]
        if length(species) == nreactants + nproducts
            return species[1:nreactants], species[(nreactants + 1):end]
        end
    end

    # Conservative fallback: keep the raw STARLIB order if this chapter is not
    # supported yet. We will expand this as the network grows.
    return species, String[]
end

function _supported_starlib_layout(chapter::Int, nspecies::Int)
    return (chapter == 1 && nspecies == 2) ||
           (chapter in (2, 4) && nspecies == 3) ||
           (chapter in (3, 5, 8) && nspecies == 4) ||
           (chapter in (6, 9) && nspecies == 5) ||
           (chapter in (7, 10) && nspecies == 6)
end

function _supported_rate_table(table::ReactionRateTable)
    return !isempty(table.reactants) && !isempty(table.products)
end

function _valid_nuclear_bookkeeping(table::ReactionRateTable)
    _supported_rate_table(table) || return false
    return try
        reaction_conservation(Reaction(table)).valid_nuclear_bookkeeping
    catch
        false
    end
end

"""
    starlib_chapter_report(tables)

Summarize which STARLIB chapter layouts were parsed into supported reactant and
product bookkeeping. Unsupported rows are kept by `read_starlib`, but they are
not suitable for network construction until their chapter layout is implemented.
"""
function starlib_chapter_report(tables::AbstractVector{ReactionRateTable})
    supported_by_chapter = Dict{Int,Int}()
    unsupported_by_chapter = Dict{Int,Int}()
    unsupported_examples = ReactionRateTable[]

    for table in tables
        if _supported_rate_table(table)
            supported_by_chapter[table.chapter] = get(supported_by_chapter, table.chapter, 0) + 1
        else
            unsupported_by_chapter[table.chapter] = get(unsupported_by_chapter, table.chapter, 0) + 1
            length(unsupported_examples) < 10 && push!(unsupported_examples, table)
        end
    end

    return (
        total=length(tables),
        supported=sum(values(supported_by_chapter); init=0),
        unsupported=sum(values(unsupported_by_chapter); init=0),
        supported_by_chapter=supported_by_chapter,
        unsupported_by_chapter=unsupported_by_chapter,
        unsupported_examples=unsupported_examples,
    )
end


"""
    read_starlib(path=DEFAULT_STARLIB_PATH; warn_unsupported=false, skip_lines=0)

Read a STARLIB v6-format rate library into `ReactionRateTable` entries.

STARLIB stores each reaction as a fixed-size block: one header line
(`chapter reactant/product-symbols... source q_value`) followed by exactly
`STARLIB_ROWS_PER_REACTION` data rows of `T9 rate factor_uncertainty`
(`rate` is the STARLIB-recommended central value at that grid point;
`factor_uncertainty` is its multiplicative lognormal factor uncertainty, the
basis for Monte Carlo sampling via `sampled_interpolate_rate`). The chapter
number encodes the reactant/product counts (`_split_starlib_species`); rows
whose chapter/species-count combination isn't one of the ones
`_supported_starlib_layout` recognizes are still kept (as raw reactants with
empty products) rather than dropped, so `starlib_chapter_report` can later
show exactly what coverage gap exists instead of silently losing data. Pass
`warn_unsupported=true` to also emit an `@warn` summary immediately.

Some targeted STARLIB-derived updates (e.g. `data/starlib/starlib_etr25_2025.txt`)
are otherwise identical to this format but carry one extra leading summary
line (`"n reactions = N, n temp coordinates = M"`); pass `skip_lines=1` to
skip it rather than have it fail as a malformed reaction header. Combine such
a file with a full library via `override_rate_tables`.
"""
# A handful of newer targeted STARLIB updates (e.g. starlib_etr25_2025.txt)
# carry an isolated doubled-minus-sign typo on a Q-value field
# ("--6.98832e-01" where a single "-" was clearly meant); tolerate it rather
# than failing the whole file over two malformed characters.
function _parse_starlib_float(s::AbstractString)
    fixed = startswith(s, "--") ? "-" * s[3:end] : s
    return parse(Float64, fixed)
end

function read_starlib(path::AbstractString=_default_starlib_path(); warn_unsupported::Bool=false, skip_lines::Integer=0)
    tables = ReactionRateTable[]
    unsupported_counts = Dict{Int,Int}()

    open(path, "r") do io
        for _ in 1:skip_lines
            eof(io) && break
            readline(io)
        end
        line_number = skip_lines
        while !eof(io)
            header = strip(readline(io))
            line_number += 1
            isempty(header) && continue

            fields = split(header)
            length(fields) >= 4 || error("Malformed STARLIB header at line $line_number: $header")

            chapter = parse(Int, fields[1])
            source = fields[end-1]
            q_value = _parse_starlib_float(fields[end])
            species = normalize_species_name.(fields[2:end-2])
            if !_supported_starlib_layout(chapter, length(species))
                unsupported_counts[chapter] = get(unsupported_counts, chapter, 0) + 1
            end
            reactants, products = _split_starlib_species(chapter, species)

            T9 = Float64[]
            rate = Float64[]
            factor_uncertainty = Float64[]
            sizehint!(T9, STARLIB_ROWS_PER_REACTION)
            sizehint!(rate, STARLIB_ROWS_PER_REACTION)
            sizehint!(factor_uncertainty, STARLIB_ROWS_PER_REACTION)

            for _ in 1:STARLIB_ROWS_PER_REACTION
                eof(io) && error("Unexpected end of file after STARLIB header at line $line_number")
                row = split(strip(readline(io)))
                line_number += 1
                length(row) >= 3 || error("Malformed STARLIB data row at line $line_number")

                push!(T9, _parse_starlib_float(row[1]))
                push!(rate, _parse_starlib_float(row[2]))
                push!(factor_uncertainty, _parse_starlib_float(row[3]))
            end

            push!(tables, ReactionRateTable(chapter, reactants, products, source, q_value, T9, rate, factor_uncertainty))
        end
    end

    if warn_unsupported && !isempty(unsupported_counts)
        summary = join(["chapter $chapter: $count" for (chapter, count) in sort(collect(unsupported_counts))], ", ")
        @warn "Unsupported STARLIB chapter layouts were kept as raw reactants with empty products: $summary"
    end

    return tables
end

