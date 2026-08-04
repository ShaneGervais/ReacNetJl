#=
Run a whole-network STARLIB/NACRE lognormal Monte Carlo sensitivity study:
baseline first, then `--trials` independent MC trials, each sampling every
reaction that carries real rate-uncertainty data from `p ~ Normal(0, 1)`
(reactions without uncertainty data stay at their nominal/recommended rate --
flagged once here, up front, not per trial).

## Why lognormal-only, not per-parameter PDFs

The physically complete version of this study samples individual nuclear
input parameters -- resonance energies (Gaussian), resonance strengths /
partial widths / non-resonant S-factors (lognormal), upper limits of Γ and
γω (Porter-Thomas), interference (binary) -- and propagates them through the
rate formula, the way STARLIB's own rates were originally generated
(Longland et al. 2010, RatesMC). That requires STARLIB's per-resonance INPUT
parameter tables. This repo (and everywhere else searched on this machine)
only has STARLIB's aggregate OUTPUT: a (T9, recommended rate, lognormal
factor uncertainty) grid per reaction -- the already-collapsed result of
that Monte Carlo, not its inputs. Sampling `rate * factor_uncertainty^p`,
`p ~ Normal(0,1)`, from that aggregate factor uncertainty is the correct and
only rigorous method the data we actually have supports; it's also exactly
how STARLIB's documentation says the factor uncertainty is meant to be used.
See `ReacNetJl.sample_rate_p_values_all` / `has_rate_uncertainty`.

## Output layout

    <output_base>/baseline/            the unfactored run
    <output_base>/mc/trial_0001/       trial 1 (seed = seed_base + 1)
    <output_base>/mc/trial_0002/       trial 2 (seed = seed_base + 2)
    ...

Each trial is its own `run_ppn.jl --mc-all` subprocess (isolation + its own
run.log, which also records that trial's sampled-vs-unsampled reaction list),
bounded to `--jobs` concurrent processes at a time.

Usage:
    julia --project=. SensitivityStudy/run_monte_carlo.jl <nova_case> <output_base> \
        [--option 1|2|3] [--trials N] [--jobs N] \
        [--screening weak|chugunov|none] [--seed-base N]

Defaults: option=3, trials=100, jobs=4, screening=chugunov, seed-base=0

Example:
    julia --project=. SensitivityStudy/run_monte_carlo.jl \
        ne_nova_1.15_12_X_weiss_mixed SensitivityStudy/mc_runs --option 3 --trials 200 --jobs 8
=#

using ReacNetJl

const SCRIPT_DIR = @__DIR__
const PROJECT_ROOT = dirname(SCRIPT_DIR)

include(joinpath(SCRIPT_DIR, "rate_options.jl"))

function parse_args(args)
    positional = String[]
    option = 3
    trials = 100
    jobs = 4
    screening_arg = "chugunov"
    seed_base = 0
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--option"
            option = parse(Int, args[i+1]); i += 2
        elseif arg == "--trials"
            trials = parse(Int, args[i+1]); i += 2
        elseif arg == "--jobs"
            jobs = parse(Int, args[i+1]); i += 2
        elseif arg == "--screening"
            screening_arg = args[i+1]; i += 2
        elseif arg == "--seed-base"
            seed_base = parse(Int, args[i+1]); i += 2
        else
            push!(positional, arg); i += 1
        end
    end
    return positional, option, trials, jobs, screening_arg, seed_base
end

const USAGE = """
Usage: julia --project=. SensitivityStudy/run_monte_carlo.jl <nova_case> <output_base> \\
    [--option 1|2|3] [--trials N] [--jobs N] \\
    [--screening weak|chugunov|none] [--seed-base N]
"""

positional, option, ntrials, njobs, screening_arg, seed_base = parse_args(ARGS)
length(positional) >= 2 || error(USAGE)
nova_case, output_base = positional[1], positional[2]

case_dir = joinpath(SCRIPT_DIR, "nova_cases", nova_case)
trajectory_path = joinpath(case_dir, "trajectory.input")
abundance_path = if isfile(joinpath(case_dir, "initial_abundance_jch1.dat"))
    joinpath(case_dir, "initial_abundance_jch1.dat")
elseif isfile(joinpath(case_dir, "initial_abundance.dat"))
    joinpath(case_dir, "initial_abundance.dat")
else
    error("no initial_abundance_jch1.dat or initial_abundance.dat in $case_dir")
end
isfile(trajectory_path) || error("no trajectory.input in $case_dir")

println("nova case       = ", nova_case)
println("output base     = ", output_base)
println("option          = ", option, "  (", OPTION_DESCRIPTIONS[option], ")")
println("trials          = ", ntrials)
println("parallel jobs   = ", njobs)
println("seed base       = ", seed_base)
println()

# ---------- baseline, in-process (also gives us the network to report uncertainty coverage) ----------

println("running baseline...")
baseline_dir = joinpath(output_base, "baseline")
screening = screening_arg == "none" ? nothing : Symbol(screening_arg)
baseline_result = if option == 1
    run_ppn(trajectory_path, abundance_path; rates=:iliadis2002, screening=screening, output_dir=baseline_dir)
else
    run_ppn(trajectory_path, abundance_path; tables=rate_tables_for_option(option), screening=screening, output_dir=baseline_dir)
end
println("baseline done -> ", baseline_dir)
println("  validated reactions = ", baseline_result.validation.num_reactions)
println("  network species     = ", length(baseline_result.network.species))

# Same final_state.csv snapshot run_ppn.jl's CLI writes, so decay_ppn.jl (or
# an analysis notebook) can use the baseline the same way as any other run.
let profiles = trajectory_profiles(baseline_result.trajectory)
    final_time = last(baseline_result.times)
    final_T9 = profiles.T9(final_time)
    final_all = merge(baseline_result.inert_mass_fractions, baseline_result.final_mass_fractions)
    open(joinpath(baseline_dir, "final_state.csv"), "w") do io
        println(io, "# final_time_s,", final_time)
        println(io, "# final_T9,", final_T9)
        println(io, "# option,", option)
        println(io, "species,mass_fraction")
        for (name, value) in sort(collect(final_all); by=first)
            println(io, name, ",", value)
        end
    end
end
println()

# ---------- report which reactions have real uncertainty data, once, up front ----------

coverage = sample_rate_p_values_all(baseline_result.network)
n_sampled = length(coverage.sampled_reactions)
n_unsampled = length(coverage.unsampled_reactions)
println(n_sampled, " / ", n_sampled + n_unsampled, " reactions carry real rate-uncertainty data and will be sampled each trial.")
println("WARNING: ", n_unsampled, " reaction(s) have no uncertainty data and will stay at their nominal rate every trial:")
for label in coverage.unsampled_reactions
    println("  - ", label)
end
println()

# ---------- dispatch the trials, one run_ppn.jl subprocess per trial ----------

run_ppn_jl = joinpath(SCRIPT_DIR, "run_ppn.jl")
mc_dir = joinpath(output_base, "mc")

println(ntrials, " trial(s) to launch across ", njobs, " parallel job(s)...")
println()

function run_trial(trial)
    seed = seed_base + trial
    trial_dir = joinpath(mc_dir, "trial_" * lpad(trial, 4, '0'))
    mkpath(trial_dir)
    log_path = joinpath(trial_dir, "run.log")
    cmd = `julia --project=$PROJECT_ROOT $run_ppn_jl $trajectory_path $abundance_path $trial_dir --option $option --mc-all --seed $seed --screening $screening_arg`
    ok = open(log_path, "w") do io
        try
            run(pipeline(cmd; stdout=io, stderr=io))
            true
        catch
            false
        end
    end
    return (trial=trial, ok=ok, dir=trial_dir)
end

results = asyncmap(run_trial, 1:ntrials; ntasks=max(1, njobs))

failed = [r for r in results if !r.ok]
println()
println("Monte Carlo sweep complete: ", ntrials - length(failed), " / ", ntrials, " trials succeeded")
if !isempty(failed)
    println("FAILED trials (see run.log in each trial dir for details):")
    for r in failed
        println("  - trial ", r.trial, " -> ", r.dir)
    end
end
