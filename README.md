# ReacNetJl

ReacNetJl is an early-stage Julia project for building single-zone and eventually multi-zone post-processing nuclear reaction networks for thermonuclear nova environments.

The long-term goal is to provide a user-friendly, fast, and physically accurate tool for nuclear astrophysics reaction-rate uncertainty studies using experimental STARLIB data.

## Current status

The project currently supports the single-zone network workflow:

```text
STARLIB table -> Reaction -> ReactionNetwork -> dY/dt -> time evolution
```

Implemented so far:

- STARLIB-style rate table reading
- reaction label parsing, e.g. `18F(p,α)15O`, `8B(β+)2α`, and `p(p,eν)d`
- isotope/species name normalization
- `ReactionRateTable`
- `Reaction`
- `ReactionNetwork`
- rate interpolation in `log(rate)` versus `log(T9)`
- reaction flux calculation
- network right-hand side calculation:

```math
\frac{dY_i}{dt} = \sum_r \left(\nu^{\mathrm{prod}}_{i,r} - \nu^{\mathrm{react}}_{i,r}\right) F_r
```

## Reaction-rate data

All rate libraries live in `data/`, which is entirely gitignored (the raw
files are too large and too numerous for Git) except for a single tracked
archive, `data.zip`. Get everything in one step:

```sh
unzip data.zip
```

This gives you the full STARLIB v6.10 library, the JINA REACLIB snapshots,
the raw NPDATA/NACRE-NetGen source material, and the small derived
paper-tabulated rate overrides (`data/iliadis2001_rates.dat`,
`data/nacre_rates.dat`) with the scripts that built them
(`data/scripts/build_iliadis2001_rates.jl`, `build_nacre_rates.jl`) --
reproducible from their original sources if the raw material is ever updated.

If you only need the essentials (no NPDATA/NetGen extras, no paper-tabulated
overrides), fetch just the directly-downloadable files instead:

```julia
using ReacNetJl
ReacNetJl.fetch_data!()   # REACLIB (ReaclibV1.0 + current default), winvne, AME2020, STARLIB v6.10
```

Rate sources supported:

- **STARLIB** (`data/starlib.dat`): tabulated rates with factor uncertainties,
  used for Monte Carlo uncertainty sampling. Targeted STARLIB-derived updates
  (e.g. a newer remeasurement of a handful of reactions) can be read with
  `read_starlib(path; skip_lines=..., extra_trailing_fields=...)` and layered
  onto a full library with `override_rate_tables`.
- **JINA REACLIB**: analytic fits evaluated onto the STARLIB temperature grid.
  `data/reaclib_v1.0.dat` is the frozen ReaclibV1.0 snapshot, which already
  contains the complete Iliadis 2001 and NACRE fit sets used below.
- **Paper-tabulated overrides** (`data/iliadis2001_rates.dat`,
  `data/nacre_rates.dat`, read via `read_iliadis2001_rates`/
  `read_nacre_rates`): the literal published rate values from Iliadis et al.
  (2001) and NACRE (Angulo et al. 1999), used to override the REACLIB fits of
  the same reactions where the two differ -- REACLIB's own fit for a given
  label doesn't always match its source paper closely at nova temperatures.

The Iliadis et al. (2002, ApJS 142, 105) nova baseline is built with:

```julia
using ReacNetJl

result = iliadis2002_rate_tables()   # NACRE A<20, Iliadis 2001 A=20-40
result.report.counts                 # tables per category: :nacr, :il01, :weak, :other
result.report.fallbacks              # reactions outside both compilations
```

Every returned table records its REACLIB set label in `table.source`
(`"nacr"`, `"il01"`, `"wc12w"`, ...), so rate provenance stays explicit all the
way into the network.

## Quick single-zone workflow

The highest-level API is `run_ppn`, which takes a trajectory file and an
initial-abundance file and runs the full pipeline (network selection, reverse
rates, validation, adaptive implicit solve):

```julia
using ReacNetJl

result = run_ppn("trajectory.input", "initial_abundance.dat";
                 rates=:iliadis2002, screening=:chugunov,
                 output_dir="outputs/run")

result.final_mass_fractions
result.solver_stats
result.rate_policy_report.counts
result.output_files   # mass_fractions.csv, reaction_fluxes.csv, integrated_fluxes.csv, network.csv
```

Pass `rate_factors=Dict("22Na(p,γ)23Mg" => 2.0)` for a deterministic
sensitivity-study factor on one named reaction, or
`rate_sample_labels=["16O(p,γ)17F"]` to sample that reaction's rate from its
own STARLIB/NACRE lognormal factor uncertainty instead (fixed for the run;
pass your own `rng=MersenneTwister(seed)` for reproducibility). Both leave
every other reaction at its nominal rate.

`examples/run_ppn.jl` and `examples/decay_ppn.jl` wrap this into ready-to-run
CLI scripts:

```sh
julia --project=. examples/run_ppn.jl trajectory.input initial_abundance.dat outputs/run \
    --option 1 --screening chugunov --factor "22Na(p,g)23Mg=2.0"

julia --project=. examples/decay_ppn.jl outputs/run 7200 outputs/decay
```

`--option` selects the rate library cumulatively: `1` = Iliadis (2002)'s own
NACRE+Iliadis-2001 baseline, `2` = option 1 with a targeted ~2013-era STARLIB
update layered on top, `3` = option 2 with the most up-to-date targeted
update layered on top again (see `examples/rate_options.jl`) -- the reaction
*set* stays identical across all three so the comparison isolates how much
the rate *values* have changed over time. `run_ppn.jl` saves which option it
used inside `outputs/run/final_state.csv`, so `decay_ppn.jl` matches it
automatically unless you pass `--option` explicitly to override.

For label-driven interactive experiments, use `solve_single_zone`:

```julia
using ReacNetJl

tables = read_starlib()

result = solve_single_zone(
    tables,
    ["18F(p,α)15O", "18F(p,γ)19Ne"],
    Dict("p" => 0.70, "he4" => 0.28, "18F" => 1.0e-5, "15O" => 0.0, "19Ne" => 0.0),
    (0.0, 1.0e-3),
    1.0e-5,
    1.0e3,
    0.2;
    adaptive=true,
    method=:backward_euler,
    screening=:weak,
    max_fractional_change=0.05,
    max_absolute_change=1.0e-10,
)

println(result.final_mass_fractions)
println(result.integrated_fluxes)
println(result.integrated_energy_generation)
println(result.mass_fraction_drift)
println(result.abundance_diagnostics)
println(result.solver_stats)
```

This builds the network from STARLIB labels, validates reaction bookkeeping,
converts mass fractions to abundances, evolves the one-zone ODE, and returns
raw abundance histories, mass-fraction diagnostics, positivity diagnostics, and
solver statistics.

## Single-Zone PPN Scope

The current target is standard single-zone post-processing nucleosynthesis
(PPN). In this mode ReacNetJl takes prescribed thermodynamic histories:

```text
T9(t), rho(t), initial composition -> network solve -> abundances and diagnostics
```

The network does not change the trajectory temperature or density. Energy
generation is diagnostic only:

```text
epsilon_nuc(t) = N_A * MeV_to_erg * sum_r Q_r F_r(t)
```

This is the right scope for validating rates, reaction flow, abundance output,
screening choices, and STARLIB uncertainty propagation against a given nova
trajectory. Convection, mixing between zones, and hydrodynamic feedback belong
to later MPPN/TPPN stages, not the current single-zone PPN core.

## Current Simulation Example

`examples/mini_nova_network.jl` now uses `trajectory.input` from the project
root when present. The file is parsed with:

- `AGEUNIT = YRS`, converted to seconds
- `TUNIT = T9K`
- `RHOUNIT = CGS`

The current mini nova example uses `initial_abundance.dat` or
`initial_abundance.DAT` from the project root when present. The abundance file
is normalized before use, then species outside the active network are reported
as inert/outside-network mass.

The current mini nova example uses a compact trajectory-driven network over
CNO, NeNa, MgAl, Si-Ca, and a Ca-Fe/Ni seed extension, with explicit `26Al*`
isomer channels. It runs with backward Euler and weak screening:

```text
validated reactions = example-dependent
species = example-dependent
screening = weak
active network initial mass is reported
inert/outside-network initial mass is reported
total active-network mass is conserved
Newton accepted/rejected step statistics are reported
peak epsilon_nuc and integrated nuclear energy are printed as diagnostics
```

This is now a real trajectory-driven post-processing calculation if
`trajectory.input` is present, but it is still not a production nova model until
the network, screening model, reverse rates, and energy feedback are validated
for the target nova regime.

## Full Single-Zone PPN Output

`examples/single_zone_nova_ppn.jl` is the current user-facing single-zone PPN
output driver. It reads:

- `trajectory.input` from the project root
- `initial_abundance.DAT` or `initial_abundance.dat` from the project root
- `iso_massf00804.DAT` as an optional output-format template when present

Run it with:

```sh
julia --project=. examples/single_zone_nova_ppn.jl
```

For a threaded run, pass `--jobs`:

```sh
julia --project=. examples/single_zone_nova_ppn.jl --jobs 8
```

The script relaunches itself with eight Julia threads when needed. This
parallelizes backward-Euler finite-difference Jacobian columns and parallelizes
the per-timestep `.DAT` output writing. The one-zone ODE still advances forward
in time sequentially, so this is not equivalent to splitting timesteps across
workers.

The script builds an expanded H-Ca nova network from STARLIB, adds available
reverse rates, solves it with backward Euler, and interpolates the abundance
history onto the trajectory time grid. For the current 805-row trajectory it
writes:

```text
outputs/single_zone_nova_ppn/iso_massf00000.DAT
...
outputs/single_zone_nova_ppn/iso_massf00804.DAT
outputs/single_zone_nova_ppn/mass_fractions.csv
```

The `.DAT` files use mass fractions under the `ABUNDANCE_MF` column and carry
diagnostic header values for timestep, age, `T9`, `rho`, number-density markers,
and nuclear energy generation. Isotopes outside the active reaction network are
kept at their normalized initial mass fractions; active-network isotopes are
evolved by the solver. The generated `outputs/` directory is ignored by Git.

For a quick formatting check without writing every trajectory state:

```sh
julia --project=. examples/single_zone_nova_ppn.jl --output-stride 100
```

To run with the Iliadis-2002 REACLIB baseline instead of STARLIB:

```sh
julia --project=. examples/single_zone_nova_ppn.jl --rates iliadis2002
```

Neutron-induced reactions are included in the network selection by default;
pass `--no-neutron-captures` to exclude them. The driver prints the rate
library, the projectile list, the policy category counts, and the source-label
histogram of the reactions actually selected into the active network.

Useful comparison and diagnostic options:

```sh
julia --project=. examples/single_zone_nova_ppn.jl --jobs 8 --screening none --decay-time 3600
```

This can write `iso_massfDECAY.DAT`, comparison CSVs against the Iliadis JCH1
reference when present, and integrated flux reports for the largest abundance
residuals.

## Physics Roadmap

Current physics scope:

- Explicit reverse-rate lookup with `find_reverse_rate`, for rates already
  present in STARLIB.
- Generated detailed-balance reverse rates for supported capture reactions.
- Approximate weak charged-particle screening as a multiplicative rate factor,
  with `--screening none` available for baseline comparisons.
- `26Al*`/`26Alm` parsing to STARLIB's `al*6` isomer entries, including isomer
  beta decay and proton capture in the mini network.
- Post-processing decay over a user-supplied decay time. This is written as a
  separate `iso_massfDECAY.DAT` file and does not change the hydrodynamic
  trajectory solve.

Known limits:

- Generated detailed-balance reverse rates include spin factors and
  partition-function ratios when `data/winvne_v2.0.dat` is present; without
  it they fall back to the mass-factor/Boltzmann approximation. Explicit
  reverse rates from the libraries are always preferred.
- `screening=:chugunov` covers the weak-to-strong regimes following Chugunov,
  DeWitt & Yakovlev (2007); `screening=:weak` remains the Salpeter-style
  multiplier for comparisons. Electron degeneracy corrections beyond that
  prescription are not modeled.
- Energy feedback coupling is out of scope for standard single-zone PPN. It
  would require a separate self-heating one-zone mode with `dT/dt =
  (epsilon_nuc - losses) / c_P`, an equation of state, heat capacity, and an
  expansion/cooling model.

## Core equations

The abundance of species `i` is

```math
Y_i = \frac{X_i}{A_i}
```

where `X_i` is mass fraction and `A_i` is mass number.

For a one-body reaction,

```math
F_r = R_r(T_9)Y_i
```

For a two-body reaction,

```math
F_r = \rho R_r(T_9)Y_iY_j
```

For repeated reactants, a symmetry factor is applied:

```math
F_r = \frac{\rho^{N_r-1}R_r(T_9)\prod_j Y_j^{\nu^{\mathrm{react}}_{j,r}}}{\prod_j \nu^{\mathrm{react}}_{j,r}!}
```

The network evolves

```math
\frac{dY_i}{dt} = \sum_r \left(\nu^{\mathrm{prod}}_{i,r} - \nu^{\mathrm{react}}_{i,r}\right)F_r
```

For a single-zone post-processing network, this is an ODE system, not a PDE, because there are no spatial derivatives yet.

## Current validation command

Run tests with:

```sh
julia --project=. test/runtests.jl
```
