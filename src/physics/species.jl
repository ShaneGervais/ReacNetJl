# Species: the fundamental Z/A-driven particle and isotope model, plus name/label parsing.

"""
    Species

A nuclear species identified by charge number `Z`, mass number `A`, and a display
`name`, such as `Species("f18", 9, 18)`.
"""
struct Species
    name::String
    Z::Int
    A::Int
end

const AVOGADRO = 6.02214076e23
const MEV_TO_ERG = 1.602176634e-6

const _PARTICLE_ALIASES = Dict(
    "p" => "p",
    "proton" => "p",
    "h1" => "p",
    "n" => "n",
    "neutron" => "n",
    "d" => "d",
    "h2" => "d",
    "t" => "t",
    "h3" => "t",
    "a" => "he4",
    "alpha" => "he4",
    "he4" => "he4",
    "g" => "gamma",
    "gamma" => "gamma",
)

const _SPECIAL_SPECIES = Dict(
    "n" => Species("n", 0, 1),
    "p" => Species("p", 1, 1),
    "d" => Species("d", 1, 2),
    "t" => Species("t", 1, 3),
    "gamma" => Species("gamma", 0, 0),
)

const _ELEMENT_SYMBOLS = [
    "h", "he", "li", "be", "b", "c", "n", "o", "f", "ne",
    "na", "mg", "al", "si", "p", "s", "cl", "ar", "k", "ca",
    "sc", "ti", "v", "cr", "mn", "fe", "co", "ni", "cu", "zn",
    "ga", "ge", "as", "se", "br", "kr", "rb", "sr", "y", "zr",
    "nb", "mo", "tc", "ru", "rh", "pd", "ag", "cd", "in", "sn",
    "sb", "te", "i", "xe", "cs", "ba", "la", "ce", "pr", "nd",
    "pm", "sm", "eu", "gd", "tb", "dy", "ho", "er", "tm", "yb",
    "lu", "hf", "ta", "w", "re", "os", "ir", "pt", "au", "hg",
    "tl", "pb", "bi", "po", "at", "rn", "fr", "ra", "ac", "th",
    "pa", "u", "np", "pu", "am", "cm", "bk", "cf", "es", "fm",
    "md", "no", "lr", "rf", "db", "sg", "bh", "hs", "mt", "ds",
    "rg", "cn", "nh", "fl", "mc", "lv", "ts", "og",
]

const _ELEMENT_Z = Dict(symbol => Z for (Z, symbol) in pairs(_ELEMENT_SYMBOLS))

"""
    normalize_species_name(name)

Normalize common isotope notation to STARLIB-style lowercase names.

Examples:
- `"18F" -> "f18"`
- `"4He" -> "he4"`
- `"α" -> "he4"`
- `"γ" -> "gamma"`
"""
function normalize_species_name(name::AbstractString)
    raw = lowercase(strip(name))
    raw == "al-6" && return "al26"

    s = raw
    s = replace(s, " " => "", "-" => "", "_" => "")
    s = replace(s, "α" => "alpha", "γ" => "gamma")

    if haskey(_PARTICLE_ALIASES, s)
        return _PARTICLE_ALIASES[s]
    end

    m = match(r"^(\d+)([a-z]+)\*$", s)
    if m !== nothing
        mass_number = m.captures[1]
        symbol = m.captures[2]
        return symbol * "*" * last(mass_number)
    end

    m = match(r"^(\d+)alm$", s)
    if m !== nothing
        mass_number = m.captures[1]
        return "al*" * last(mass_number)
    end

    m = match(r"^(\d+)alg$", s)
    if m !== nothing
        return "al" * m.captures[1]
    end

    m = match(r"^(\d+)([a-z]+)$", s)
    if m !== nothing
        return m.captures[2] * m.captures[1]
    end

    return s
end

"""
    species_from_name(name)

Create a `Species` object from a normalized or common isotope name.

Examples:
- `species_from_name("18F") == Species("f18", 9, 18)`
- `species_from_name("he4") == Species("he4", 2, 4)`
- `species_from_name("p") == Species("p", 1, 1)`
"""
function species_from_name(name::AbstractString)
    normalized = normalize_species_name(name)

    if haskey(_SPECIAL_SPECIES, normalized)
        return _SPECIAL_SPECIES[normalized]
    end

    m = match(r"^([a-z]+)(\d+)$", normalized)
    if m === nothing
        star_match = match(r"^([a-z]+)\*(\d+)$", normalized)
        if star_match !== nothing
            symbol = star_match.captures[1]
            mass_suffix = star_match.captures[2]
            haskey(_ELEMENT_Z, symbol) || throw(ArgumentError("unknown element symbol '$symbol' in species '$name'"))
            A = symbol == "al" && mass_suffix == "6" ? 26 : parse(Int, mass_suffix)
            return Species(normalized, _ELEMENT_Z[symbol], A)
        end
    end
    m === nothing && throw(ArgumentError("could not infer species information from '$name'"))

    symbol = m.captures[1]
    mass_suffix = m.captures[2]
    A = symbol == "al" && mass_suffix == "6" ? 26 : parse(Int, mass_suffix)
    haskey(_ELEMENT_Z, symbol) || throw(ArgumentError("unknown element symbol '$symbol' in species '$name'"))

    return Species(normalized, _ELEMENT_Z[symbol], A)
end

# Recognize the ejectile slot of a decay-style label, e.g. "18F(β+)18O": the
# token inside the parentheses names the decay mode, not a real projectile
# particle, so parse_reaction_label must special-case it before falling
# through to the general target(projectile,ejectile)product grammar.
function _is_beta_token(token::AbstractString)
    s = lowercase(strip(token))
    s = replace(s, " " => "", "β" => "beta")
    return s in ("beta+", "b+", "β+", "betaplus", "ec", "electroncapture", "beta-", "b-", "β-", "betaminus")
end

function _parse_particle_multiplicity(token::AbstractString)
    s = lowercase(strip(token))
    s = replace(s, " " => "", "α" => "alpha")
    m = match(r"^([2-9])([a-z\+]+)$", s)
    m === nothing && return nothing

    count = parse(Int, m.captures[1])
    particle = m.captures[2]
    particle in ("p", "n", "d", "t", "alpha") || return nothing
    normalized = normalize_species_name(particle)
    return fill(normalized, count)
end

function _is_non_nuclear_ejectile_token(token::AbstractString)
    s = lowercase(strip(token))
    s = replace(s, " " => "", "ν" => "nu")
    return s in ("enu", "e+nu", "e-nu", "nu", "neutrino")
end

function _parse_reaction_token(token::AbstractString; allow_multiplicity::Bool=false)
    _is_non_nuclear_ejectile_token(token) && return String[]

    if allow_multiplicity
        particles = _parse_particle_multiplicity(token)
        particles !== nothing && return particles
    end

    normalized = normalize_species_name(token)
    normalized == "gamma" && return String[]
    return [normalized]
end

"""
    parse_reaction_label(label)

Parse a reaction label of the form `target(projectile,ejectile)product` (the
usual nuclear-astrophysics shorthand, e.g. `"18F(p,α)15O"`), or a decay-style
label `parent(β+)daughter` / `parent(β-)daughter`. Returns `(reactants,
products)`, where both entries are normalized species-name vectors (see
`normalize_species_name`).

This is the inverse-ish counterpart to `reaction_string`. It accepts several
notational variants for convenience:
- Greek or ASCII particle symbols (`α`/`a`/`alpha`, `γ`/`g`/`gamma`).
- Neutrino-emitting weak captures like `"p(p,eν)d"` (the `eν` ejectile
  token is dropped -- it carries no baryon number to place in `products`).
- Multiplicity-prefixed ejectiles like `"2n"` or `"3α"` via
  `_parse_particle_multiplicity`, expanding to repeated entries in the
  returned vector (matching how `_symmetry_factor` expects repeated
  reactants/products to be represented).
- `26Alg`/`26Alm` ground-state/isomer notation for the `Species`
  book-keeping.

Examples:
- `parse_reaction_label("18F(p,α)15O") == (["p", "f18"], ["he4", "o15"])`
- `parse_reaction_label("18F(p,γ)19Ne") == (["p", "f18"], ["ne19"])`
- `parse_reaction_label("p(p,eν)d") == (["p", "p"], ["d"])`
"""
function parse_reaction_label(label::AbstractString)
    decay_match = match(r"^\s*([^\(]+)\(([^,\)]+)\)(.+?)\s*$", label)
    if decay_match !== nothing && _is_beta_token(decay_match.captures[2])
        parent = normalize_species_name(decay_match.captures[1])
        products = _parse_reaction_token(decay_match.captures[3]; allow_multiplicity=true)
        return [parent], products
    end

    m = match(r"^\s*([^\(]+)\(([^,]+),([^\)]+)\)(.+?)\s*$", label)
    m === nothing && throw(ArgumentError("reaction label must look like target(projectile,ejectile)product or parent(β+)daughter"))

    target = normalize_species_name(m.captures[1])
    projectile = normalize_species_name(m.captures[2])
    ejectiles = _parse_reaction_token(m.captures[3]; allow_multiplicity=true)
    product = normalize_species_name(m.captures[4])

    reactants = [projectile, target]
    products = vcat(ejectiles, [product])
    return reactants, products
end

"""
    abundance_from_mass_fraction(X, A)

Convert mass fraction `Xᵢ` to abundance `Yᵢ = Xᵢ / Aᵢ` (mol g^-1) for one
species of mass number `A`. The atomic building block of
`abundances_from_mass_fractions`.
"""
abundance_from_mass_fraction(X::Real, A::Integer) = X / A

"""
    mass_fraction_from_abundance(Y, A)

Convert abundance `Yᵢ` to mass fraction `Xᵢ = Aᵢ Yᵢ` for one species of mass
number `A` -- the inverse of `abundance_from_mass_fraction`, and the atomic
building block of `mass_fractions_from_abundances`.
"""
mass_fraction_from_abundance(Y::Real, A::Integer) = A * Y

