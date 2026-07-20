using ReacNetJl
using Printf

const NETGEN_DIR = "/home/sgervais/Documents/ReacNetJl/data/netgen"
const REC_PATH = joinpath(NETGEN_DIR, "netgen_nacre_log_100_rec.txt")
const LOW_PATH = joinpath(NETGEN_DIR, "netgen_nacre_log_100_low.txt")
const UP_PATH = joinpath(NETGEN_DIR, "netgen_nacre_log_100_up.txt")
const OUT_PATH = "/home/sgervais/Documents/ReacNetJl/data/nacre_rates.dat"

struct NetgenBlock
    reactant_tokens::Vector{Tuple{Int,String}}  # (count, "SYM A") pairs, e.g. (2, "H 1")
    product_tokens::Vector{Tuple{Int,String}}
    T8::Vector{Float64}
    rate::Vector{Float64}
end

function _parse_species_fields(fields::Vector{<:AbstractString})
    length(fields) < 2 && return nothing
    sym = fields[2]
    sym == "OOOOO" && return nothing
    sym == "PROT" && return "p"
    sym == "NEUT" && return "n"
    sym == "DEUT" && return "d"
    sym == "TRIT" && return "t"
    # Al-26 rows come in three variants for the same target+projectile: an
    # unresolved "AL 26" total (ground+isomer combined, like the Iliadis 2001
    # tables' "_total" columns), plus separately-resolved "AL 26g"
    # (ground-state) and "AL 26m" (isomer) columns. normalize_species_name
    # expects mass-number-first input ("26alg"/"26alm") for these, so handle
    # them directly here; the unresolved total is tagged with a sentinel
    # name so the whole block can be skipped in favor of the resolved g/m
    # entries (see is_known_bad_entry).
    if sym == "AL" && length(fields) >= 3
        mass_field = lowercase(fields[3])
        endswith(mass_field, "g") && return "al" * mass_field[1:end-1]
        endswith(mass_field, "m") && return "al*" * last(mass_field[1:end-1])
        mass_field == "26" && return "al26_unresolved_total"
    end

    # Most entries split symbol and mass number into separate fields
    # ("H   1", "HE  3"), but some (observed for He-3/He-4 specifically)
    # combine them into one token ("HE3", "HE4") with no third field.
    combined = length(fields) >= 3 ? string(sym, fields[3]) : sym
    return normalize_species_name(combined)
end

function parse_netgen_file(path::AbstractString)
    lines = readlines(path)
    blocks = NetgenBlock[]
    i = 1
    n = length(lines)
    while i <= n
        line = lines[i]
        if startswith(strip(line), "#") && occursin(r"^\s*\d+\s+\S+", strip(lstrip(line, '#')))
            # header block: 7 lines (reactant1, +, reactant2, =, product1, +, product2)
            header_fields = [split(strip(lstrip(lines[i+k], '#'))) for k in 0:6]
            reactant_tokens = Tuple{Int,String}[]
            product_tokens = Tuple{Int,String}[]
            for fields in (header_fields[1], header_fields[3])
                length(fields) < 2 && continue
                count = parse(Int, fields[1])
                name = _parse_species_fields(fields)
                name === nothing && continue
                push!(reactant_tokens, (count, name))
            end
            for fields in (header_fields[5], header_fields[7])
                length(fields) < 2 && continue
                count = parse(Int, fields[1])
                name = _parse_species_fields(fields)
                name === nothing && continue
                push!(product_tokens, (count, name))
            end

            i += 7
            # skip to "T8" marker then blank comment line, then data rows
            while i <= n && strip(lines[i]) != "#       T8" && !occursin("T8", lines[i])
                i += 1
            end
            i += 1 # skip "#" separator after T8 label
            while i <= n && strip(lines[i]) == "#"
                i += 1
            end
            T8 = Float64[]
            rate = Float64[]
            while i <= n && !isempty(strip(lines[i])) && !startswith(strip(lines[i]), "#")
                fields = split(strip(lines[i]))
                push!(T8, parse(Float64, fields[1]))
                push!(rate, parse(Float64, fields[2]))
                i += 1
            end
            push!(blocks, NetgenBlock(reactant_tokens, product_tokens, T8, rate))
        else
            i += 1
        end
    end
    return blocks
end

println("parsing rec...")
rec_blocks = parse_netgen_file(REC_PATH)
println("parsing low...")
low_blocks = parse_netgen_file(LOW_PATH)
println("parsing up...")
up_blocks = parse_netgen_file(UP_PATH)
println("rec=", length(rec_blocks), " low=", length(low_blocks), " up=", length(up_blocks))

function expand_species(tokens::Vector{Tuple{Int,String}})
    names = String[]
    for (count, name) in tokens
        append!(names, fill(name, count))
    end
    return names
end

function reaction_key(block::NetgenBlock)
    return (sort(expand_species(block.reactant_tokens)), sort(expand_species(block.product_tokens)))
end

# Check baryon-number balance on the *nuclear* participants alone (weak
# processes that emit a lepton, like p+p->d+e+nu, will show a mass-number
# mismatch here since only nucleons are tracked -- these must come from
# REACLIB instead, not this NACRE strong-reaction table).
function mass_balanced(block::NetgenBlock)
    reactant_A = sum(species_from_name(n).A for n in expand_species(block.reactant_tokens); init=0)
    product_A = sum(species_from_name(n).A for n in expand_species(block.product_tokens); init=0)
    return reactant_A == product_A
end

# netgen_nacre's very first entry (Qrad=1.172, Qnu=0.270 MeV, matching
# p+p->d+e+nu's known Q=1.442 MeV) mislabels its product "HE  2" instead of
# "H   2" -- a labeling quirk isolated to this single row (confirmed: no
# other "HE  2" occurrence in the file), not real NACRE strong-reaction data.
# This process is a weak capture (bundled into netgen's NACRE distribution
# via a separate Goriely 1999 weak-rate index, per the file's own readme),
# already correctly present via REACLIB. Exclude it here rather than
# "fixing" the label, since it isn't NACRE data to begin with.
function is_known_bad_entry(block::NetgenBlock)
    species = vcat(expand_species(block.reactant_tokens), expand_species(block.product_tokens))
    return "he2" in species || "al26_unresolved_total" in species
end

function write_nacre_rates(io)
    n_written = 0
    n_skipped_unbalanced = 0
    n_skipped_nomatch = 0
    n_skipped_baddata = 0
    n_skipped_known_bad = 0
    println(io, "# NACRE (Angulo et al. 1999, Nucl. Phys. A 656, 3) recommended reaction rates,")
    println(io, "# with lower/upper bounds as factor uncertainties.")
    println(io, "# Reformatted from data/netgen/netgen_nacre_log_100_{rec,low,up}.txt")
    println(io, "# (the NACRE-project NetGen distribution format) into ReacNetJl's per-reaction block format.")

    for (idx, rec) in enumerate(rec_blocks)
        if is_known_bad_entry(rec)
            n_skipped_known_bad += 1
            continue
        end
        if !mass_balanced(rec)
            n_skipped_unbalanced += 1
            continue
        end
        key = reaction_key(rec)
        low_idx = findfirst(b -> reaction_key(b) == key, low_blocks)
        up_idx = findfirst(b -> reaction_key(b) == key, up_blocks)
        if low_idx === nothing || up_idx === nothing
            n_skipped_nomatch += 1
            continue
        end
        low = low_blocks[low_idx]
        up = up_blocks[up_idx]
        if length(rec.T8) != length(low.T8) || length(rec.T8) != length(up.T8)
            n_skipped_baddata += 1
            continue
        end

        reactants = expand_species(rec.reactant_tokens)
        products = expand_species(rec.product_tokens)
        label = join(reactants, "+") * "_to_" * join(products, "+")

        println(io, "reaction: $(join(reactants, " ")) -> $(join(products, " ")) ; label: $label ; table: nacre ; variant: standard")
        for k in eachindex(rec.T8)
            t9 = rec.T8[k] / 10.0
            r = rec.rate[k]
            lo = low.rate[k]
            hi = up.rate[k]
            if r <= 0.0
                continue
            end
            if lo > 0.0 && hi >= lo
                println(io, "$t9 $r $lo $hi")
            else
                println(io, "$t9 $r")
            end
        end
        println(io, "end")
        n_written += 1
    end

    println("written = ", n_written)
    println("skipped (known-bad he2 label quirk) = ", n_skipped_known_bad)
    println("skipped (mass-unbalanced, likely weak) = ", n_skipped_unbalanced)
    println("skipped (no low/up match) = ", n_skipped_nomatch)
    println("skipped (bad grid alignment) = ", n_skipped_baddata)
end

open(write_nacre_rates, OUT_PATH, "w")
