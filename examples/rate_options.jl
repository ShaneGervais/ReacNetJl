#=
Rate library selection for the --option {1,2,3} CLI flag, shared by
run_ppn.jl and decay_ppn.jl so both scripts resolve the same library for the
same option number.

Option 1: Iliadis et al. (2002) baseline -- NACRE (A<20) + Iliadis et al.
          2001 (A=20-40), the same REACLIB label sets the JCH1 sensitivity
          study itself used.
Option 2: STARLIB v6.10, the full modern STARLIB snapshot we have on disk
          (data/starlib.dat). NOTE: this is the 2022 point-release, not
          literally the original 2013 STARLIB release (Sallaska et al.) --
          say so if you need that exact vintage snapshot for historical
          fidelity; we don't currently have it on disk.
Option 3: STARLIB v6.10 with the etr25 (2025) targeted update layered on top
          via override_rate_tables -- the most up-to-date rates we have
          today. NOTE: mc10/mc13/taly (also under data/starlib/) aren't
          incorporated yet -- they need their own reader (an extra trailing
          field STARLIB v6 doesn't have) and a priority rule for how they'd
          combine with etr25 (taly in particular is a 32000-reaction
          Hauser-Feshbach theoretical set, presumably a last-resort filler
          rather than a blanket override).
=#

const ETR25_PATH = joinpath(dirname(@__DIR__), "data", "starlib", "starlib_etr25_2025.txt")

function rate_tables_for_option(option::Integer)
    if option == 1
        return iliadis2002_rate_tables(; include_reverse=true).tables
    elseif option == 2
        return read_starlib()
    elseif option == 3
        base = read_starlib()
        isfile(ETR25_PATH) || return base
        return override_rate_tables(base, read_starlib(ETR25_PATH; skip_lines=1))
    end
    throw(ArgumentError("unsupported --option $option; use 1 (Iliadis 2002), 2 (STARLIB v6.10), or 3 (STARLIB + etr25)"))
end

const OPTION_DESCRIPTIONS = Dict(
    1 => "Iliadis 2002: NACRE (A<20) + Iliadis 2001 (A=20-40)",
    2 => "STARLIB v6.10 (full library)",
    3 => "STARLIB v6.10 + etr25 (2025) targeted update",
)
