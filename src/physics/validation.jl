# Conservation checks and network validation reports.

# Sum mass number A or charge Z (property=:A or :Z) over a list of species
# names, e.g. sum_i A_i over a reaction's reactants -- the building block of
# the baryon/charge conservation check in `reaction_conservation`.
function _species_property_sum(names::AbstractVector{String}, property::Symbol)
    total = 0
    for name in names
        species = species_from_name(name)
        if property == :A
            total += species.A
        elseif property == :Z
            total += species.Z
        else
            throw(ArgumentError("unsupported species property $property"))
        end
    end
    return total
end

"""
    reaction_conservation(reaction)

Check baryon-number (`A`) and charge (`Z`) conservation for one reaction.
Returns a named tuple with reactant/product totals and boolean conservation flags.

Photons, if present, are treated as `A=0`, `Z=0`.
"""
function reaction_conservation(reaction::Reaction)
    reactant_A = _species_property_sum(reaction.reactants, :A)
    product_A = _species_property_sum(reaction.products, :A)
    reactant_Z = _species_property_sum(reaction.reactants, :Z)
    product_Z = _species_property_sum(reaction.products, :Z)
    delta_A = product_A - reactant_A
    delta_Z = product_Z - reactant_Z
    is_weak_decay = delta_A == 0 && abs(delta_Z) == 1

    return (
        reaction=join(reaction.reactants, "+") * "->" * join(reaction.products, "+"),
        reactant_A=reactant_A,
        product_A=product_A,
        delta_A=delta_A,
        conserves_A=reactant_A == product_A,
        reactant_Z=reactant_Z,
        product_Z=product_Z,
        delta_Z=delta_Z,
        conserves_Z=reactant_Z == product_Z,
        is_weak_decay=is_weak_decay,
        valid_nuclear_bookkeeping=(reactant_A == product_A) && ((reactant_Z == product_Z) || is_weak_decay),
    )
end

function _network_species_issues(network::ReactionNetwork)
    issues = String[]

    if length(unique(network.species)) != length(network.species)
        push!(issues, "network species list contains duplicates")
    end

    if length(network.species_info) != length(network.species)
        push!(issues, "network species_info length does not match species length")
    end

    for (i, name) in pairs(network.species)
        haskey(network.species_index, name) || push!(issues, "species '$name' is missing from species_index")
        get(network.species_index, name, 0) == i || push!(issues, "species_index for '$name' does not match species order")

        if i <= length(network.species_info)
            info = network.species_info[i]
            info.name == name || push!(issues, "species_info[$i] name '$(info.name)' does not match species '$name'")
        end
    end

    return issues
end

function _reaction_validation_report(network::ReactionNetwork, reaction_index::Int, reaction::Reaction)
    issues = String[]

    for name in Iterators.flatten((reaction.reactants, reaction.products))
        haskey(network.species_index, name) || push!(issues, "species '$name' is not present in the network")
    end

    conservation = try
        reaction_conservation(reaction)
    catch err
        push!(issues, "could not check conservation: $(err)")
        nothing
    end

    if conservation !== nothing
        conservation.conserves_A || push!(issues, "baryon number is not conserved: ΔA=$(conservation.delta_A)")
        if !conservation.conserves_Z && !conservation.is_weak_decay
            push!(issues, "charge is not conserved: ΔZ=$(conservation.delta_Z)")
        end
    end

    return (
        reaction_index=reaction_index,
        reaction=reaction_string(reaction),
        reactants=copy(reaction.reactants),
        products=copy(reaction.products),
        conservation=conservation,
        valid=isempty(issues),
        issues=issues,
    )
end

"""
    network_validation_report(network; throw_on_error=false)

Validate basic reaction-network bookkeeping and physics constraints.

Checks currently include:
- consistency of `species`, `species_info`, and `species_index`
- all reaction species are present in the network
- baryon-number conservation for every reaction
- charge conservation for every reaction

Returns a named tuple containing summary information, per-reaction reports, and
all issues. If `throw_on_error=true`, invalid networks raise an `ArgumentError`.
"""
function network_validation_report(network::ReactionNetwork; throw_on_error::Bool=false)
    issues = _network_species_issues(network)
    reaction_reports = NamedTuple[]

    for (reaction_index, reaction) in pairs(network.reactions)
        report = _reaction_validation_report(network, reaction_index, reaction)
        push!(reaction_reports, report)
        for issue in report.issues
            push!(issues, "reaction $(reaction_index) ($(reaction_string(reaction))): $issue")
        end
    end

    report = (
        valid=isempty(issues),
        num_species=length(network.species),
        num_reactions=length(network.reactions),
        species=copy(network.species),
        reaction_reports=reaction_reports,
        issues=issues,
    )

    if throw_on_error && !report.valid
        throw(ArgumentError("network validation failed:\n" * join(report.issues, "\n")))
    end

    return report
end

