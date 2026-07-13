# Human-readable reaction labels.

#=
    reaction_string(reaction)

Return a compact, normalized reaction label for display or diagnostics.

For common two-body STARLIB reactions, this returns labels like
`"f18(p,he4)o15"` or `"f18(p,γ)ne19"`.
=#
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

