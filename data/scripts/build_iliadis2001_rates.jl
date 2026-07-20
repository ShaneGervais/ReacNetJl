using ReacNetJl

const SRC_DIR = "/home/sgervais/Documents/NuGrid/nuppn/bckup/NovaProject_0.2/references"
const OUT_PATH = "/home/sgervais/Documents/ReacNetJl/data/iliadis2001_rates.dat"
const TABLE_NUMBERS = 3:9

struct RateColumn
    label::String
    reactants::Vector{String}
    products::Vector{String}
    T9::Vector{Float64}
    rate::Vector{Float64}
end

function parse_table(path::AbstractString, table_number::Int)
    lines = readlines(path)
    header_idx = findfirst(l -> startswith(l, "T9_GK"), lines)
    header_idx === nothing && error("no header row found in $path")
    header_fields = split(lines[header_idx], '\t')
    labels = String.(header_fields[2:end])

    columns = RateColumn[]
    for (col, label) in enumerate(labels)
        endswith(label, "_total") && continue  # redundant sum of g/m siblings, skip
        reactants, products = try
            parse_reaction_label(label)
        catch e
            @warn "could not parse label, skipping" path label exception=e
            continue
        end
        push!(columns, RateColumn(label, reactants, products, Float64[], Float64[]))
    end

    label_to_column = Dict(c.label => c for c in columns)
    for line in lines[(header_idx + 1):end]
        isempty(strip(line)) && continue
        fields = split(line, '\t')
        t9 = parse(Float64, fields[1])
        for (col, label) in enumerate(labels)
            haskey(label_to_column, label) || continue
            raw = strip(fields[col + 1])
            raw == "..." && continue
            value = parse(Float64, raw)
            value > 0.0 || continue
            column = label_to_column[label]
            push!(column.T9, t9)
            push!(column.rate, value)
        end
    end

    return columns
end

open(OUT_PATH, "w") do io
    println(io, "# Iliadis et al. (2001, ApJS 134, 151), Tables 3-9: recommended reaction rates.")
    println(io, "# Reformatted from ~/Documents/NuGrid/nuppn/bckup/NovaProject_0.2/references/iliadis2001_table03-09.dat")
    println(io, "# (digitized from the paper's PDF tables) into ReacNetJl's per-reaction block format.")
    total_reactions = 0
    total_rows = 0
    for table_number in TABLE_NUMBERS
        path = joinpath(SRC_DIR, "iliadis2001_table0$(table_number).dat")
        columns = parse_table(path, table_number)
        for c in columns
            isempty(c.T9) && continue
            reactant_str = join(c.reactants, " ")
            product_str = join(c.products, " ")
            println(io, "reaction: $reactant_str -> $product_str ; label: $(c.label) ; table: $table_number ; variant: standard")
            for (t9, rate) in zip(c.T9, c.rate)
                println(io, "$t9 $rate")
            end
            println(io, "end")
            total_reactions += 1
            total_rows += length(c.T9)
        end
    end
    println("wrote $total_reactions reactions, $total_rows data rows")
end
println("wrote ", OUT_PATH)
