#=
Network completeness audit (Milestone 2), self-contained: no nuppn dependency.

For the same rate library used in the JCH1 run (iliadis2002_rate_tables()),
re-derive select_h_ca_reaction_tables' selection decision for every candidate
reaction that touches one of our active network's species, and report *why*
each one was kept or dropped: out of the physical Z/A/(N-Z) band, invalid
mass/charge bookkeeping, no seed projectile, or unreachable (a valid
candidate whose reactants never all became available during the selection
closure).

Run from the ReacNetJl project root:
    julia --project=. SensitivityStudy/network_audit.jl
=#

using ReacNetJl
using Printf

const CASE_DIR = joinpath(@__DIR__, "..", "nova_cases", "ne_nova_1.15_12_X_weiss_mixed")
const TRAJECTORY_PATH = joinpath(CASE_DIR, "trajectory.input")
const ABUNDANCE_PATH = joinpath(CASE_DIR, "initial_abundance_jch1.dat")

result = run_ppn(TRAJECTORY_PATH, ABUNDANCE_PATH; rates=:iliadis2002, screening=:chugunov)
active_species = Set(result.network.species)
println("active network species = ", length(active_species))

tables = iliadis2002_rate_tables().tables
seed_species = Set(normalize_species_name.(collect(keys(read_initial_abundances(ABUNDANCE_PATH)))))
push!(seed_species, "n"); push!(seed_species, "p"); push!(seed_species, "he4")

# --- Reproduce select_h_ca_reaction_tables' filtering/closure, but record *why*. ---
physical(name) = ReacNetJl._h_ca_candidate_species(name; zmax=20, amax=46, proton_rich_margin=4, neutron_rich_margin=3)
const PROJECTILES = ("p", "he4", "he3", "d", "n")  # matches run_ppn's default neutron_captures=true

function classify_candidacy(table)
    ReacNetJl._valid_nuclear_bookkeeping(table) || return :invalid_bookkeeping
    all(physical, vcat(table.reactants, table.products)) || return :outside_physical_band
    if length(table.reactants) == 1 && ReacNetJl._weak_source(table.source)
        return :ok
    end
    any(p -> p in table.reactants, PROJECTILES) && return :ok
    return :no_seed_projectile
end

# Only look at reactions that touch our active species at all -- reactions
# entirely among species we never produce are out of scope for this audit.
touches_active(table) = any(s -> s in active_species, vcat(table.reactants, table.products))

candidacy_reason = Dict{Any,Symbol}()
by_reason_examples = Dict{Symbol,Vector{ReactionRateTable}}()
candidates = ReactionRateTable[]
for table in tables
    touches_active(table) || continue
    reason = classify_candidacy(table)
    key = ReacNetJl._reaction_participant_key(table)
    candidacy_reason[key] = reason
    if reason == :ok
        push!(candidates, table)
    else
        vec = get!(by_reason_examples, reason, ReactionRateTable[])
        length(vec) < 12 && push!(vec, table)
    end
end

# Re-run the same fixed-point closure used inside select_h_ca_reaction_tables.
selected_keys = Set{Any}()
species = copy(seed_species)
changed = true
while changed
    global changed = false
    for table in candidates
        all(r -> r in species, table.reactants) || continue
        key = ReacNetJl._reaction_participant_key(table)
        key in selected_keys && continue
        push!(selected_keys, key)
        for product in table.products
            if !(product in species)
                push!(species, product)
                global changed = true
            end
        end
    end
end

unreachable = [t for t in candidates if !(ReacNetJl._reaction_participant_key(t) in selected_keys)]

println("\n=== Candidate reactions touching active species, by exclusion reason ===")
for (reason, examples) in sort(collect(by_reason_examples); by=first)
    total = count(t -> classify_candidacy(t) == reason, filter(touches_active, tables))
    println("\n-- ", reason, "  (", total, " reactions total, showing up to 12) --")
    for t in examples
        @printf("  %-30s source=%-8s reactants=%-20s products=%s\n",
            reaction_string(Reaction(t)), t.source, join(t.reactants, "+"), join(t.products, "+"))
    end
end

println("\n-- unreachable (valid candidate, but reactants never all available)  (", length(unreachable), " reactions, showing up to 12) --")
for t in unreachable[1:min(12, length(unreachable))]
    missing_reactants = [r for r in t.reactants if !(r in species)]
    @printf("  %-30s source=%-8s missing_reactant(s)=%s\n",
        reaction_string(Reaction(t)), t.source, join(missing_reactants, ","))
end

println("\nselected (in final network) = ", length(selected_keys))
println("candidates but unreachable   = ", length(unreachable))
for (reason, examples) in by_reason_examples
    total = count(t -> classify_candidacy(t) == reason, filter(touches_active, tables))
    println(reason, " = ", total)
end
