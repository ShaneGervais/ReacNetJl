# Human-readable reaction labels.

# Render one species' internal name (e.g. "f18", "al*6", "he4") as the
# conventional nuclear-astrophysics display form (e.g. "18F", "26Al*", "α"),
# the inverse direction of `normalize_species_name`. Light particles get
# their usual shorthand (p, n, d, t, α, 3He); isomers keep their trailing
# `*`; everything else is `<A><Symbol>` with the element symbol capitalized.
function _display_species_name(name::AbstractString)
    normalized = normalize_species_name(name)
    normalized == "p" && return "p"
    normalized == "n" && return "n"
    normalized == "d" && return "d"
    normalized == "t" && return "t"
    normalized == "he4" && return "α"
    normalized == "he3" && return "3He"

    m = match(r"^([a-z]+)\*(\d+)$", normalized)
    if m !== nothing
        species = species_from_name(normalized)
        return string(species.A, uppercasefirst(m.captures[1]), "*")
    end

    species = species_from_name(normalized)
    symbol_match = match(r"^([a-z]+)\d+$", normalized)
    symbol = symbol_match === nothing ? normalized : symbol_match.captures[1]
    return string(species.A, uppercasefirst(symbol))
end

# Render a list of species as a display group joined by "+", collapsing
# repeats into a "<count><name>" prefix (e.g. two protons -> "2p", not
# "p+p"). An empty list (a pure disintegration with no light-particle
# ejectile) displays as the emitted photon, "γ".
function _format_species_group(names::AbstractVector{String})
    isempty(names) && return "γ"

    parts = String[]
    seen = Set{String}()
    for name in names
        normalized = normalize_species_name(name)
        normalized in seen && continue
        multiplicity = Base.count(==(normalized), normalize_species_name.(names))
        label = _display_species_name(normalized)
        push!(parts, multiplicity == 1 ? label : string(multiplicity, label))
        push!(seen, normalized)
    end
    return join(parts, "+")
end

function _weak_symbol(delta_Z::Int)
    delta_Z == -1 && return "β+"
    delta_Z == 1 && return "β-"
    return "weak"
end

"""
    reaction_string(reaction)

Return a compact, human-readable reaction label for display or diagnostics,
in conventional nuclear-astrophysics notation.

Dispatches on `reaction_conservation(reaction)` to pick the right notation:
- weak decays (`is_weak_decay`) as `"13N(β+)13C"`, or `"target(projectile,eν)products"`
  for two-reactant weak captures (e.g. electron capture);
- ordinary two-body captures as `"18F(p,γ)19Ne"` (single product) or
  `"18F(p,α)15O"` (light-particle ejectile plus one heavy product, the last
  product listed is treated as the "main" residual nucleus);
- anything else (three-body reactions, multi-product disintegrations) falls
  back to the plain `"reactants->products"` form.

This is the inverse-ish counterpart to `parse_reaction_label`: labels this
function prints are not guaranteed to round-trip through
`parse_reaction_label` byte-for-byte (Greek/Unicode glyphs vs. the parser's
ASCII shorthand), but they describe the same reaction.
"""
function reaction_string(reaction::Reaction)
    conservation = reaction_conservation(reaction)

    if conservation.is_weak_decay
        weak = _weak_symbol(conservation.delta_Z)
        if length(reaction.reactants) == 1
            return string(_format_species_group(reaction.reactants), "(", weak, ")", _format_species_group(reaction.products))
        elseif length(reaction.reactants) == 2
            projectile, target = reaction.reactants
            return string(
                _display_species_name(target),
                "(",
                _display_species_name(projectile),
                ",eν)",
                _format_species_group(reaction.products),
            )
        end
        return string(_format_species_group(reaction.reactants), "(", weak, ")->", _format_species_group(reaction.products))
    end

    if length(reaction.reactants) == 2
        projectile, target = reaction.reactants
        if length(reaction.products) == 1
            return string(_display_species_name(target), "(", _display_species_name(projectile), ",γ)", _format_species_group(reaction.products))
        elseif length(reaction.products) >= 2
            ejectiles = reaction.products[1:end-1]
            product = reaction.products[end:end]
            return string(
                _display_species_name(target),
                "(",
                _display_species_name(projectile),
                ",",
                _format_species_group(ejectiles),
                ")",
                _format_species_group(product),
            )
        end
    end

    return _format_species_group(reaction.reactants) * "->" * _format_species_group(reaction.products)
end

