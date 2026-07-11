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

- fixed-step explicit Euler integration
- fixed-step RK4 integration
- fixed-step backward Euler integration with finite-difference Newton iterations
- constant or time-dependent `rho(t)` and `T9(t)`
- adaptive explicit timestep control
- STARLIB factor-uncertainty sampling
- Monte Carlo uncertainty runs
- trajectory file reading and interpolation
- metadata-aware trajectory input with `AGEUNIT`, `TUNIT`, and `RHOUNIT`
- initial abundance table parsing via `read_initial_abundances`
- flux diagnostics and integrated reaction flows
- diagnostic nuclear energy generation from Q-values
- total mass-fraction history and drift diagnostics
- abundance positivity diagnostics
- approximate weak charged-particle screening multiplier via `screening=:weak`
- baryon-number and charge-conservation checks
- explicit STARLIB reverse-rate lookup with `find_reverse_rate`
- generated detailed-balance reverse rates for supported capture reactions
- unsupported STARLIB chapter-layout reporting
- STARLIB `26Al*` / `26Alm` isomer label support for the `al*6` entries
- multiproduct STARLIB chapter support for beta-delayed proton/alpha channels
- H-Ca nova network selection with active/inert/missing isotope diagnostics
- post-processing weak decay output via `iso_massfDECAY.DAT`
- Iliadis JCH1 comparison CSVs and terminal residual summaries
- integrated flux reports for the largest abundance residuals
- precomputed integer stoichiometry inside `ReactionNetwork`
- adaptive solver diagnostics for accepted/rejected steps, timestep range, and Newton iterations
- one-call single-zone post-processing from reaction labels and mass fractions
- full trajectory-indexed PPN abundance output as `iso_massfXXXXX.DAT` files plus a wide mass-fraction CSV
- threaded finite-difference Jacobian construction for backward Euler when Julia is started with multiple threads
- JINA REACLIB reading (`read_reaclib`), analytic fit evaluation (`reaclib_rate`), and evaluation onto the STARLIB temperature grid (`reaclib_rate_tables`)
- Iliadis-2002 baseline rate policy (`iliadis2002_rate_tables`): NACRE for A < 20, Iliadis et al. 2001 for A = 20-40, with weak rates, fallback labels, and a full provenance report
- per-step rate caching and an analytic Jacobian for backward Euler (`jacobian=:analytic`, the default; `:finite_difference` retained for validation)
- winvne nuclear data import (`read_winvne`): ground-state spins and partition functions applied to generated detailed-balance reverse rates
- Chugunov, DeWitt & Yakovlev (2007) screening (`screening=:chugunov`), valid from the weak- into the strong-screening regime, alongside the Salpeter `:weak` option
- one-call `run_ppn(trajectory, abundances; rates=..., screening=...)` post-processing API
- optional high-order stiff integration via the OrdinaryDiffEqBDF package extension (`solve_network_fbdf`, `run_ppn(...; method=:fbdf)`)
- a PrecompileTools workload so fresh sessions start with the solver hot path already compiled
- REPL docstrings for the public API (`?run_ppn`, `?iliadis2002_rate_tables`, ...)
- tests for parsing, interpolation, fluxes, RHS calculation, time evolution, solver equivalence, and user-facing workflows

## Reaction-rate data

All rate libraries live in `data/` (see `data/README.md`). The large files are
not tracked by Git; fetch the REACLIB files with:

```sh
./data/download_rates.sh
```

Two rate sources are supported:

- **STARLIB** (`data/starlib.dat`): tabulated rates with factor uncertainties,
  used for Monte Carlo uncertainty sampling.
- **JINA REACLIB**: analytic fits evaluated onto the STARLIB temperature grid.
  `data/reaclib_v1.0.dat` is the frozen ReaclibV1.0 snapshot,
  `data/reaclib_il01.dat` and `data/reaclib_nacr.dat` are the complete
  Iliadis 2001 and NACRE fit sets held by the JINA database.

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
                 rates=:iliadis2002, screening=:chugunov)

result.final_mass_fractions
result.solver_stats
result.rate_policy_report.counts
```

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

## Milestones

### Milestone 1: Iliadis baseline reproduction

Goal: make the single-zone PPN output comparable to the Iliadis JCH1 baseline
for the same initial composition, thermodynamic history, and decay convention.

Tasks:

- Reconstruct or replace the digitized `trajectory.input` with the actual
  one-zone thermodynamic history used for the Iliadis comparison.
- Track final-state and post-decay residuals against
  `outputs/iliadis2002_jch1/iso_massf00000.DAT`.
- Use the integrated flux report to identify whether each major residual is
  caused by trajectory mismatch, missing rate data, missing reaction channels,
  or solver/network behavior.
- Add a stable comparison artifact for notebooks: isotope, reference mass
  fraction, model mass fraction, ratio, log residual, group, and active/inert
  status.

### Milestone 2: H-Ca network completeness audit

Goal: verify that the active network can represent the same physical reaction
space as a 142-isotope nova post-processing network from hydrogen through
calcium.

Tasks:

- Keep the H-Ca selector close to the Iliadis isotope/reaction scope without
  silently adding irrelevant heavy-ion or neutron-only channels.
- Produce a tracked network audit table for selected, skipped, unsupported, and
  missing STARLIB rows.
- Separate skipped rows by reason: unsupported chapter layout, no active
  species path, disintegration bookkeeping, neutron-induced branch, heavy-ion
  branch, or missing rate data.
- Add tests that pin the expected parser behavior for multiproduct weak decays,
  pp-chain light-particle reactions, and common nova breakout branches.

### Milestone 3: Reaction-rate data provenance

Goal: make every important reaction rate traceable and replace weak STARLIB
coverage with better data where needed.

Tasks:

- Done: REACLIB import path (`read_reaclib`, `reaclib_rate_tables`) and the
  `iliadis2002_rate_tables` policy with per-reaction provenance reporting.
- Done: the PPN driver prints the source-label histogram of the active
  network and prefers explicit REACLIB reverse fits of the chosen label over
  generated detailed-balance rates in `--rates iliadis2002` mode.
- Done: the compilation coverage gap is closed with tabulated rates extracted
  from the papers themselves (`data/iliadis2001_rates.dat`, all 55 reactions;
  `data/nacre_rates.dat`, 72 of 86 with lower/upper limits as factor
  uncertainties), cross-validated against the independent JINA fits and
  applied as `il01tab`/`nacrtab` overrides by `iliadis2002_rate_tables()`.
  See `data/README.md` for the validation methodology.
- Open: identify the h-burning and He-burning reactions whose STARLIB rates
  are missing or clearly not suitable for the Iliadis baseline.

### Milestone 4: Post-decay policy

Goal: make post-calculation decay physically explicit rather than hard-coded to
an Iliadis-specific mode.

Tasks:

- Keep `--decay-time S` as the user-facing control.
- Write post-decay output only to `iso_massfDECAY.DAT`, leaving the final
  network state untouched.
- Add a decay report showing active weak branches, half-lives/rates, parents,
  daughters, and mass moved during the post-decay interval.
- Decide how to handle unresolved long-lived isomers and branches that require
  data not present in the selected rate tables.

### Milestone 5: Solver robustness and performance

Goal: keep the expanded network stable enough for routine local nova runs while
preserving transparent diagnostics.

Tasks:

- Done: analytic Jacobian assembly with per-step rate caching replaces the
  finite-difference default (`examples/benchmark_solver.jl` measured a
  ~360x speedup on the 144-species nova network; `jacobian=:finite_difference`
  remains for validation).
- Done: optional FBDF stiff integration as a package extension. Install the
  solver once with `import Pkg; Pkg.add("OrdinaryDiffEqBDF")`, then
  `using OrdinaryDiffEqBDF` enables `solve_network_fbdf` and
  `run_ppn(...; method=:fbdf)`. Plain ReacNetJl carries no heavy solver
  dependencies.
- Move to sparse factorization once networks exceed a few hundred species.
- Add optional flux-history and energy-history CSV outputs beside the
  mass-fraction CSV.
- Track Newton failure modes by timestep, temperature, density, and dominant
  reaction flux.
- Evaluate whether SciML/DifferentialEquations.jl is worth adding once the
  network physics and data provenance are stable.

### Milestone 6: Beyond single-zone PPN

Goal: prepare for MPPN/TPPN without mixing it into the current validation
target.

Tasks:

- Keep single-zone PPN trajectory-prescribed and energy-diagnostic only.
- Design the MPPN data model: zones, zone masses, thermodynamic histories, and
  mixing coefficients.
- Define what outputs should remain compatible with the current
  `iso_massfXXXXX.DAT` format.

## Learning path

Useful topics to study while building:

1. Mass fraction `X_i` vs abundance `Y_i`
2. Thermonuclear reaction-rate units
3. STARLIB chapter/reaction types
4. Stoichiometry and reaction-network ODEs
5. Explicit vs implicit ODE methods
6. Stiffness in nuclear reaction networks
7. Lognormal uncertainty sampling
8. Monte Carlo uncertainty propagation

## Current validation command

Run tests with:

```sh
julia --project=. test/runtests.jl
```
