# ReacNetJl

ReacNetJl is an early-stage Julia project for building single-zone and eventually multi-zone post-processing nuclear reaction networks for thermonuclear nova environments.

The long-term goal is to provide a user-friendly, fast, and physically accurate tool for nuclear astrophysics reaction-rate uncertainty studies using experimental STARLIB data.

## Current status

The project currently supports the first minimal single-zone network workflow:

```text
STARLIB table -> Reaction -> ReactionNetwork -> dY/dt -> time evolution
```

Implemented so far:

- STARLIB-style rate table reading
- reaction label parsing, e.g. `18F(p,α)15O`
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
- constant or time-dependent `rho(t)` and `T9(t)`
- basic rate multiplier support for future uncertainty studies
- tests for parsing, interpolation, fluxes, RHS calculation, and time evolution

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

## Roadmap milestones

### Milestone 1: Species registry and abundance helpers

Goal: make initial conditions physically grounded and easy to use.

Add:

- `species_from_name(name)`
- a species registry for common particles and nuclei
- `Z` and `A` inference from names like `f18`, `o15`, `he4`, `p`
- `abundances_from_mass_fractions(network, Xdict)`
- `mass_fractions_from_abundances(network, Y)`
- optional mass-fraction normalization checks

Why this matters:

Users usually specify initial composition as mass fractions. The solver evolves abundances. This milestone gives a safe bridge between user input and solver input.

---

### Milestone 2: Network construction from reaction labels

Goal: make network creation easy from STARLIB data.

Add:

- `reaction_from_label(tables, label; source=nothing)`
- `network_from_labels(tables, labels; species=nothing, source=nothing)`
- automatic species inference from selected reactions
- helpful errors when a STARLIB reaction is missing or ambiguous

Example target API:

```julia
tables = read_starlib()

network = network_from_labels(
    tables,
    ["18F(p,α)15O", "18F(p,γ)19Ne"],
)
```

---

### Milestone 3: First real single-zone example

Goal: demonstrate the current network on a physically relevant nova reaction.

Add an example script:

```text
examples/single_zone_18f.jl
```

The example should:

- read STARLIB
- build a small `18F` destruction network
- define nova-like fixed conditions, e.g. `T9 = 0.2`, `rho = 1000.0`
- define initial mass fractions
- convert mass fractions to abundances
- evolve the network
- print final mass fractions
- optionally save a CSV-like output file

Possible first network:

```text
18F(p,α)15O
18F(p,γ)19Ne
```

---

### Milestone 4: Adaptive timestep control

Goal: make explicit integration safer without adding dependencies yet.

Add a simple adaptive timestep option that limits the maximum fractional abundance change per step:

```math
\max_i \left|\frac{\Delta Y_i}{Y_i}\right| < \epsilon
```

Possible API:

```julia
solve_network_adaptive(
    network,
    Y0,
    tspan,
    dt_initial,
    rho,
    T9;
    max_fractional_change=0.01,
)
```

This is not a replacement for a stiff solver, but it will help avoid unstable fixed-step runs while the project remains dependency-light.

---

### Milestone 5: STARLIB uncertainty sampling

Goal: use STARLIB factor uncertainties correctly.

STARLIB rates are commonly sampled as lognormal variations:

```math
R_{\mathrm{sampled}}(T) = R_{\mathrm{recommended}}(T) f_{\mathrm{unc}}(T)^p
```

where

```math
p \sim \mathcal{N}(0,1)
```

Add:

- interpolation of `factor_uncertainty`
- `sampled_interpolate_rate(table, T9, p)`
- reaction-level uncertainty parameter storage
- solver support for `rate_p_values`

Possible API:

```julia
p_values = randn(length(network.reactions))

times, history = solve_network(
    network,
    Y0,
    tspan,
    dt,
    rho,
    T9;
    rate_p_values=p_values,
)
```

---

### Milestone 6: Monte Carlo driver

Goal: run repeated networks for uncertainty evaluation.

Add:

- `run_monte_carlo`
- random sampling of STARLIB rate uncertainties
- collection of final abundances
- simple percentile summaries
- reproducible random seeds

Example target API:

```julia
results = run_monte_carlo(
    network,
    Y0,
    tspan,
    dt,
    rho,
    T9;
    nruns=1000,
    seed=1234,
)
```

Output should eventually support:

- final abundance distributions
- median values
- 16th/84th percentile intervals
- sensitivity diagnostics

---

### Milestone 7: Conservation and physics checks

Goal: catch physically invalid reactions or bookkeeping mistakes.

Add:

- baryon number conservation checks
- charge conservation checks
- warnings for unsupported STARLIB chapters
- optional network validation report

For reaction

```text
18F(p,α)15O
```

check:

```math
18 + 1 = 4 + 15
```

and

```math
9 + 1 = 2 + 8
```

---

### Milestone 8: Energy generation diagnostics

Goal: compute nuclear energy generation from reaction Q-values.

Add:

```math
\epsilon \propto \sum_r Q_r F_r
```

This will initially be diagnostic only, not coupled back to temperature.

---

### Milestone 9: Trajectory post-processing

Goal: evolve a single-zone network over hydrodynamic nova trajectories.

Add:

- trajectory reader for `(time, T9, rho)` data
- interpolation of trajectory profiles
- solve network over trajectory conditions

Target flow:

```text
trajectory file -> T9(t), rho(t) -> solve_network
```

---

### Milestone 10: Stiff solver research and implementation

Goal: move toward production-quality integration.

Explicit Euler and RK4 are useful for learning and small tests, but real nuclear networks are often stiff. Later options:

- implement a basic implicit backward Euler method
- implement Newton iterations and finite-difference Jacobians
- exploit sparse reaction-network structure
- evaluate whether using SciML/DifferentialEquations.jl is worthwhile despite adding a dependency

This should come after the network API and physical bookkeeping are stable.

## Recently completed next steps

The first eight roadmap implementation steps are now complete:

1. Added `species_from_name` and a built-in element/species registry.
2. Extended `ReactionNetwork` with `species_info`, so each network species has accessible `A` and `Z`.
3. Added mass-fraction to abundance conversion for an entire network.
4. Added abundance to mass-fraction conversion for output.
5. Added tests for species parsing and abundance conversion.
6. Added `reaction_from_label` and `network_from_labels`.
7. Created `examples/single_zone_18f.jl`.
8. Ran the example using real STARLIB data.

## Recently completed diagnostics work

The flux-analysis and graph-diagnostics implementation plan is now complete:

1. Added `reaction_string(reaction)`.
2. Added `reaction_fluxes(network, Y, rho, T9)`.
3. Added `reaction_flux_history(network, history, times, rho, T9)`.
4. Added `integrated_fluxes(times, flux_history)`.
5. Added `species_flux_balance(network, Y, rho, T9)`.
6. Added `reaction_edges(network)` for graph-like external plotting data.
7. Updated `examples/single_zone_18f.jl` to print integrated reaction fluxes.
8. Added `examples/oxygen_fluorine_mini_network.jl`.

## Recently completed validation work

The basic physics/bookkeeping validation layer is now complete:

1. Added `reaction_conservation(reaction)`.
2. Added `network_validation_report(network; throw_on_error=false)`.
3. Validation checks reaction baryon-number conservation.
4. Validation checks reaction charge conservation.
5. Validation checks consistency of `species`, `species_info`, and `species_index`.
6. Validation checks that reaction species are present in the network.
7. Examples now call `network_validation_report(network; throw_on_error=true)` before solving.

## Recently completed uncertainty work

The first STARLIB uncertainty propagation layer is now complete:

1. Added `interpolate_factor_uncertainty(table, T9)`.
2. Added `sampled_interpolate_rate(table, T9, p)`.
3. Connected sampled STARLIB rates to `reaction_flux` using `rate_p_value`.
4. Connected sampled STARLIB rates to `reaction_fluxes`, `network_rhs`, and `solve_network` using `rate_p_values`.
5. Added `run_monte_carlo(...)` for repeated single-zone uncertainty runs.
6. Added `examples/monte_carlo_18f.jl`.

## Immediate next steps

Recommended next implementation order:

1. Add adaptive timestep control for safer explicit integration.
2. Improve weak decay / beta-decay label support and STARLIB chapter handling.
3. Add optional CSV output for examples and Monte Carlo summaries.
4. Add trajectory post-processing from `(time, T9, rho)` profiles.
5. Research/implement a basic implicit method for stiff networks.
6. Expand toward a larger nova ppn network once decay and timestep handling are safer.

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
