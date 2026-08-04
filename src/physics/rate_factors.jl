# Named-reaction rate factoring: deterministic multipliers and STARLIB
# lognormal Monte Carlo sampling for one or a few specific reactions.
#
# The solver stack (network_rhs, reaction_flux, solve_network_adaptive, ...)
# already accepts `rate_multipliers`/`rate_p_values` vectors ordered like
# `network.reactions`; what's missing without this file is a convenient way
# to say "factor *this* reaction, by name" instead of having to know its
# integer index. This is the building block for the deterministic
# `--factor label=value` runs and the (later) per-reaction Monte Carlo
# sensitivity script -- `run_monte_carlo` already does whole-network
# lognormal sampling; this targets one named reaction at a time instead.

"""
    find_reaction_indices(network, reactants, products)
    find_reaction_indices(network, label)

Find the index/indices into `network.reactions` (and `network.compiled_reactions`)
of the reaction(s) matching a given reactant/product participant list, or a
human-readable label parsed with `parse_reaction_label` (e.g.
`"22Na(p,γ)23Mg"`). Matching ignores reactant/product order but not
multiplicity (`he4+he4` only matches a reaction with exactly two `he4`
reactants). Usually returns zero or one index; STARLIB/REACLIB occasionally
carry more than one rate source for the same reaction, in which case every
matching index is returned so the caller can decide (see
`rate_multipliers_from_factors`, which treats more than one match as an
error rather than guessing).
"""
function find_reaction_indices(
    network::ReactionNetwork,
    reactants::AbstractVector{<:AbstractString},
    products::AbstractVector{<:AbstractString},
)
    wanted_reactants = sort(normalize_species_name.(reactants))
    wanted_products = sort(normalize_species_name.(products))
    return [
        i for (i, reaction) in pairs(network.reactions)
        if sort(reaction.reactants) == wanted_reactants && sort(reaction.products) == wanted_products
    ]
end

function find_reaction_indices(network::ReactionNetwork, label::AbstractString)
    reactants, products = parse_reaction_label(label)
    return find_reaction_indices(network, reactants, products)
end

function _single_reaction_index(network::ReactionNetwork, label::AbstractString)
    indices = find_reaction_indices(network, label)
    isempty(indices) && throw(ArgumentError("no reaction in the network matches '$label'"))
    length(indices) > 1 && throw(ArgumentError(
        "label '$label' matches $(length(indices)) reactions in the network; ambiguous -- narrow it down or factor by index directly",
    ))
    return only(indices)
end

"""Fast Dynamical Modelling of Milky Way Globular Clusters and their Black Hole Popula
    rate_multipliers_from_factors(network, factors)

Build a `rate_multipliers` vector (one entry per `network.reactions`, default
`1.0`) for deterministic reaction-rate factoring: `factors` is a `Dict` (or
any iterable of `label => value` pairs) mapping reaction labels to a
multiplicative factor, e.g. `Dict("22Na(p,γ)23Mg" => 2.0)` doubles that one
reaction's rate while leaving every other reaction at its nominal rate. This
is the mechanism behind Iliadis-style rate-sensitivity studies (Iliadis et
al. 2002 Table 8 varies individual rates by factors like 2, 0.5, 10, 0.1).
Throws if a label matches zero or more than one reaction in the network.
"""
function rate_multipliers_from_factors(network::ReactionNetwork, factors)
    multipliers = ones(Float64, length(network.reactions))
    for (label, factor) in factors
        multipliers[_single_reaction_index(network, label)] = Float64(factor)
    end
    return multipliers
end

"""
    sample_rate_p_values(network, labels; rng=Random.default_rng())

Build a `rate_p_values` vector (one entry per `network.reactions`, default
`0.0` = nominal rate) for STARLIB lognormal Monte Carlo sampling of one or a
few *named* reactions, rather than the whole network (see `run_monte_carlo`
for whole-network sampling). For each reaction label in `labels`, draws a
fresh `p ~ Normal(0, 1)` and holds it fixed for the run; every other
reaction keeps `p = 0` (its tabulated central/recommended rate). At
evaluation time, `sampled_interpolate_rate` turns this into
`rate(T9) * factor_uncertainty(T9)^p`, i.e. the reaction's own STARLIB
factor uncertainty at the current temperature -- unlike
`rate_multipliers_from_factors`'s fixed, temperature-independent factor.

This is the setup step for a future per-reaction sensitivity-run script: call
this once per Monte Carlo trial (with a fresh `rng` draw each time) to get a
new `rate_p_values` vector for that trial's `run_ppn(...; rate_p_values=...)`
call.
"""
function sample_rate_p_values(network::ReactionNetwork, labels; rng::AbstractRNG=Random.default_rng())
    p_values = zeros(Float64, length(network.reactions))
    for label in labels
        p_values[_single_reaction_index(network, label)] = randn(rng)
    end
    return p_values
end

"""
    has_rate_uncertainty(table::ReactionRateTable)

Whether `table` carries real, usable rate-uncertainty information, as
opposed to the structural default every REACLIB-derived
(`_reaclib_group_table`) and 2-column paper-tabulated
(`read_tabulated_rates` without a lower/upper column pair, e.g. our digitized
Iliadis-2001 table) `ReactionRateTable` gets: `factor_uncertainty` uniformly
`1.0`, which makes `sampled_rate = rate * 1.0^p == rate` regardless of `p` --
i.e. "no uncertainty data" and "structurally impossible to perturb" are the
same condition here. Two kinds of source carry real values: a genuine
STARLIB rate-grid override (`mc10`/`mc13`/`etr25`/`v6.10`/... via
`read_starlib`), and NACRE's own paper-tabulated low/recommended/high bounds
(`nacrtab`, when the 4-column `data/nacre_rates.dat` format was digitized for
that reaction) -- NACRE's whole methodology is built around publishing those
three values per reaction, independent of any STARLIB grid.
"""
has_rate_uncertainty(table::ReactionRateTable) = any(!=(1.0), table.factor_uncertainty)

"""
    sample_rate_p_values_all(network; rng=Random.default_rng())

Whole-network STARLIB lognormal Monte Carlo sampling: draws a fresh
`p ~ Normal(0, 1)` for every reaction that `has_rate_uncertainty`, and holds
`p = 0` (nominal/recommended rate, unperturbed) for every other reaction --
typically the majority of a network built from a mixed REACLIB/paper-table/
STARLIB library, since only reactions with a genuine STARLIB rate-grid
override carry a real factor uncertainty at all.

This is the network-selection-aware counterpart to the lower-level
`run_monte_carlo` solver driver (which samples every reaction unconditionally
at fixed T9/rho -- harmless there too, since `factor_uncertainty=1` makes the
unsampled reactions' draws inert, just without the accounting this function
provides). Returns a named tuple:

  - `p_values`: the vector, ready for `run_ppn(...; rate_p_values=...)`.
  - `sampled_reactions` / `unsampled_reactions`: `reaction_string.(...)` labels,
    so a caller can report -- once, not per trial -- which reactions actually
    have STARLIB uncertainty data to draw from and which were left at their
    nominal rate for lack of it.
"""
function sample_rate_p_values_all(network::ReactionNetwork; rng::AbstractRNG=Random.default_rng())
    p_values = zeros(Float64, length(network.reactions))
    sampled_reactions = String[]
    unsampled_reactions = String[]
    for (i, reaction) in pairs(network.reactions)
        label = reaction_string(reaction)
        if has_rate_uncertainty(reaction.rate_table)
            p_values[i] = randn(rng)
            push!(sampled_reactions, label)
        else
            push!(unsampled_reactions, label)
        end
    end
    return (p_values=p_values, sampled_reactions=sampled_reactions, unsampled_reactions=unsampled_reactions)
end
