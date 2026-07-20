# Detailed-balance reverse-rate generation.

# Product of mass numbers prod_i A_i over a reaction side, the `Aa*Ab` /
# `Ac` building block of the detailed-balance mass factor below.
function _mass_number_factor(names::AbstractVector{String})
    factor = 1.0
    for name in names
        factor *= species_from_name(name).A
    end
    return factor
end

# Detailed balance (as implemented) only covers the simplest, most common
# nova-network case: an exothermic two-body radiative capture a+b -> c+gamma
# (Q>0), whose reverse is a one-body photodisintegration c+gamma -> a+b.
# Endothermic forward reactions, or forward reactions with more than one
# product, aren't handled by this generator (their reverse rates must come
# from the rate library directly, or `reverse_summary.missing` reports them
# as unavailable).
function _can_generate_detailed_balance_reverse(table::ReactionRateTable)
    return length(table.reactants) == 2 &&
           length(table.products) == 1 &&
           table.q_value > 0.0
end

const DETAILED_BALANCE_CAPTURE_CONSTANT = 9.86851e9
const DETAILED_BALANCE_Q_FACTOR = 11.605
const GENERATED_REVERSE_RATE_FLOOR = 1.0e-300
"""
    generated_detailed_balance_reverse_table(table; partition_functions=nothing)

Generate a detailed-balance reverse table for a radiative-capture style
two-body forward reaction `a + b -> c + \\gamma`:

```math
\\lambda_{\\mathrm{rev}}(T_9) = C\\, T_9^{3/2} \\left(\\frac{A_a A_b}{A_c}\\right)^{3/2}
    \\frac{g_a g_b}{g_c} \\frac{G_a(T_9) G_b(T_9)}{G_c(T_9)}\\,
    e^{-11.605\\, Q / T_9}\\; N_A\\langle\\sigma v\\rangle_{\\mathrm{fwd}}(T_9)
```

This is the standard nuclear-statistical-equilibrium detailed-balance
relation for a photodisintegration reverse rate: at a given temperature, the
forward capture and reverse photodisintegration rates are locked together by
the reaction Q-value and the partition functions of the participants, so the
reverse doesn't need its own independent measurement. `C` is
`DETAILED_BALANCE_CAPTURE_CONSTANT` and `11.605` is
`DETAILED_BALANCE_Q_FACTOR` (`1/k_B` in `T_9`/MeV units).

With `partition_functions` (a `PartitionFunctionTable` from `read_winvne`),
the spin factors `g = 2J+1` and normalized partition functions `G(T_9)` are
included; participants missing from the table fall back to `g\\,G = 1`.
Without `partition_functions`, only the mass-factor and Boltzmann terms are
used (`g\\,G = 1` throughout), which is the historical approximate behavior --
still qualitatively correct, but less accurate than including the real spin
and excited-state structure. Prefer explicit reverse tables from the rate
library when they exist (see `add_reverse_reaction_tables`); this generator is
the fallback for reactions with no measured/tabulated reverse rate.
"""
function generated_detailed_balance_reverse_table(
    table::ReactionRateTable;
    source::AbstractString="detail_balance",
    partition_functions=nothing,
)
    _can_generate_detailed_balance_reverse(table) || throw(ArgumentError("can only generate detailed-balance reverse rates for exothermic two-body to one-body reactions"))

    product_mass = _mass_number_factor(table.products)
    reactant_mass = _mass_number_factor(table.reactants)
    mass_factor = (reactant_mass / product_mass)^(3.0 / 2.0)

    spin_factor = 1.0
    if partition_functions !== nothing
        weights = [_statistical_weight(partition_functions, name) for name in vcat(table.reactants, table.products)]
        if !any(isnothing, weights)
            spin_factor = (weights[1] * weights[2]) / weights[3]
        end
    end

    reverse_rate = Float64[]
    sizehint!(reverse_rate, length(table.T9))
    for (T9, forward_rate) in zip(table.T9, table.rate)
        boltzmann = exp(-DETAILED_BALANCE_Q_FACTOR * table.q_value / T9)
        pf_ratio = 1.0
        if partition_functions !== nothing
            G_a = _partition_function_at(partition_functions, table.reactants[1], T9)
            G_b = _partition_function_at(partition_functions, table.reactants[2], T9)
            G_c = _partition_function_at(partition_functions, table.products[1], T9)
            if G_a !== nothing && G_b !== nothing && G_c !== nothing
                pf_ratio = (G_a * G_b) / G_c
            end
        end
        value = DETAILED_BALANCE_CAPTURE_CONSTANT * T9^(3.0 / 2.0) * mass_factor * spin_factor * pf_ratio * forward_rate * boltzmann
        push!(reverse_rate, max(value, GENERATED_REVERSE_RATE_FLOOR))
    end

    return ReactionRateTable(
        1,
        copy(table.products),
        copy(table.reactants),
        string(source),
        -table.q_value,
        copy(table.T9),
        reverse_rate,
        copy(table.factor_uncertainty),
    )
end

"""
    add_reverse_reaction_tables(all_tables, forward_tables;
                                generate_detailed_balance=true,
                                partition_functions=nothing)

Return a named tuple containing forward tables plus reverse tables. Exact
reverse entries from the rate library are preferred. If no explicit reverse
exists and the forward table is a supported radiative-capture-style `2 -> 1`
reaction, a detailed-balance reverse table is generated; pass
`partition_functions` from `read_winvne` to include spin factors and
partition-function ratios in the generated rates.
"""
function add_reverse_reaction_tables(
    all_tables::AbstractVector{ReactionRateTable},
    forward_tables::AbstractVector{ReactionRateTable};
    generate_detailed_balance::Bool=true,
    partition_functions=nothing,
)
    reverse_lookup = _reverse_table_lookup(all_tables)
    selected = _unique_reaction_tables(forward_tables)
    seen = Set{Any}(_reaction_participant_key(table) for table in selected)
    explicit = 0
    generated = 0
    missing = 0

    for table in copy(selected)
        reverse_key = (Tuple(table.products), Tuple(table.reactants))
        if haskey(reverse_lookup, reverse_key)
            reverse_table = reverse_lookup[reverse_key]
            key = _reaction_participant_key(reverse_table)
            if !(key in seen)
                push!(selected, reverse_table)
                push!(seen, key)
            end
            explicit += 1
        elseif generate_detailed_balance && _can_generate_detailed_balance_reverse(table)
            reverse_table = generated_detailed_balance_reverse_table(table; partition_functions=partition_functions)
            key = _reaction_participant_key(reverse_table)
            if !(key in seen)
                push!(selected, reverse_table)
                push!(seen, key)
            end
            generated += 1
        else
            missing += 1
        end
    end

    return (tables=selected, explicit=explicit, generated=generated, missing=missing)
end

