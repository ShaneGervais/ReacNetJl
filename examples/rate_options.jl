#=
Rate library selection for the --option {1,2,3} CLI flag, shared by
run_ppn.jl and decay_ppn.jl so both scripts resolve the same library for the
same option number.

Cumulative design: each option layers the *next era's* targeted rate
remeasurements on top of the previous option, rather than swapping in an
unrelated full library. This keeps the reaction *set* identical across all
three options (Iliadis (2002)'s own NACRE+Iliadis2001 network scope) while
letting the rate *values* evolve -- exactly what a "how much did rates
improve 2002->2013->2026" sensitivity comparison needs; mixing in STARLIB
v6.10's tens of thousands of unrelated reactions would only dilute that
comparison.

Option 1: Iliadis et al. (2002) baseline -- NACRE (A<20) + Iliadis et al.
          2001 (A=20-40), the same REACLIB label sets the JCH1 sensitivity
          study itself used.
Option 2: Option 1, with `data/starlib/starlib_mc10_mc13_082022.txt` (59
          reactions, the "mc10"/"mc13" targeted remeasurements) layered on
          top via override_rate_tables -- the ~2013 (STARLIB-era) update to
          the reactions Iliadis (2002) covered.
Option 3: Option 2, with `data/starlib/starlib_etr25_2025.txt` (78
          reactions) layered on top -- the most up-to-date rates we have
          today for those same reactions.

NOTE: `starlib_taly_012025.txt` (also under data/starlib/) is NOT
incorporated -- at 32202 reactions it's a TALYS Hauser-Feshbach theoretical
set, almost certainly meant as a last-resort filler for reactions with no
experimental rate at all, not a targeted update to layer in the same way as
mc10/mc13/etr25. Revisit if/when that's confirmed.
=#

const MC10_MC13_PATH = joinpath(dirname(@__DIR__), "data", "starlib", "starlib_mc10_mc13_082022.txt")
const ETR25_PATH = joinpath(dirname(@__DIR__), "data", "starlib", "starlib_etr25_2025.txt")

function rate_tables_for_option(option::Integer)
    tables = iliadis2002_rate_tables(; include_reverse=true).tables  # option 1 base
    option == 1 && return tables

    if isfile(MC10_MC13_PATH)
        tables = override_rate_tables(tables, read_starlib(MC10_MC13_PATH; skip_lines=1, extra_trailing_fields=1))
    end
    option == 2 && return tables

    if option == 3
        isfile(ETR25_PATH) && (tables = override_rate_tables(tables, read_starlib(ETR25_PATH; skip_lines=1)))
        return tables
    end

    throw(ArgumentError("unsupported --option $option; use 1 (Iliadis 2002), 2 (+ mc10/mc13, ~2013), or 3 (+ etr25, most up to date)"))
end

const OPTION_DESCRIPTIONS = Dict(
    1 => "Iliadis 2002: NACRE (A<20) + Iliadis 2001 (A=20-40)",
    2 => "Option 1 + mc10/mc13 (~2013 STARLIB-era update)",
    3 => "Option 2 + etr25 (2025, most up-to-date rates)",
)
