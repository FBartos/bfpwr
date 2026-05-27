# Simulation Fixture Corpus Plan

This document describes the planned GitHub-only simulation fixture corpus for
`bfpwr`. The corpus is meant to validate fixed and sequential Bayes factor
power/sample-size implementations against reusable Monte Carlo trajectories. It
is not meant to be part of CRAN checks.

## Goals

- Generate reusable simulated statistic trajectories under a deliberate set of
  design priors and adversarial settings.
- Store enough trajectory information to query fixed designs, sequential
  designs, power curves, and sample-size behavior after the simulations are run.
- Recompute Bayes factor trajectories from stored statistics using the package
  BF functions, so many analysis priors can reuse the same simulated data.
- Keep ordinary package checks fast and skip all corpus-dependent tests when the
  files are not present, especially on CRAN.
- Make the corpus reproducible, chunkable on a compute server, and easy to
  extend with new design or analysis cases.

## Core Design

The corpus separates data-generating design cases from analysis cases.

Design cases simulate data/statistic trajectories under one design prior:

```text
design_case_id
test_family
design_prior
look_grid
nsim
seed/chunk seeds
```

Analysis cases point to one design case and apply one analysis prior and set of
Bayes factor thresholds:

```text
analysis_case_id
design_case_id
bf_function
analysis_prior
thresholds
```

This gives two benefits. Different design priors are represented by genuinely
different simulated datasets, while many analysis priors can be evaluated
cheaply from the same stored statistic trajectories.

## Repository Layout

```text
simulations/
  PLAN.md
  README.md
  manifest.csv
  schema/
    design-case-v1.md
    analysis-case-v1.md
    trajectory-v1.md
    summary-v1.md
  registry/
    design-cases.R
    analysis-cases.R
  R/
    grids.R
    generators.R
    bayes-factors.R
    summarize.R
    query.R
    validate.R
  scripts/
    run_chunk.R
    merge_chunks.R
    materialize_analysis.R
    validate_case.R
    validate_corpus.R
  jobs/
    pbs-array.sh
  corpus/
    v1/
      designs/
        z/<design_case_id>/
          design.rds
          chunks/
            chunk-0001.rds
          chunk-index.rds
          chunk-index.csv
          sha256.txt
        t/<design_case_id>/
        binomial/<design_case_id>/
      analyses/
        z/<analysis_case_id>/
          analysis.rds
          outcomes.rds
          summary.rds
          summary.csv
          sha256.txt
        t/<analysis_case_id>/
        binomial/<analysis_case_id>/

package/
  inst/tinytest/
    test-simulation-fixtures-z.R
    test-simulation-fixtures-t.R
    test-simulation-fixtures-binomial.R
    test-simulation-fixtures-sequential.R
```

The full corpus lives outside `package/`. Package tests locate it through an
environment variable such as `BFPWR_SIM_CORPUS` and skip when it is absent.

## Look Grids

Use only two standard grids to keep generation and validation simple:

```r
short_grid <- seq(10, 500, by = 10)
long_grid  <- seq(10, 10000, by = 10)
```

`short_grid` is the default for broad coverage. `long_grid` is reserved for
selected adversarial cases: large-sample behavior, many-look sequential
designs, late stopping, and sample-size search validation.

All simulated look sizes should be multiples of 10. Fixed-design validation
queries the trajectory corpus by filtering to one look. Sample-size search
functions are primarily tested deterministically by checking that their returned
sample size gives the target power through the corresponding `p*()` function;
the simulation corpus can additionally validate empirical behavior at the
nearest grid look when useful.

## Simulation Model

For each replicate in a design case:

1. Draw `true_effect` from the design prior, or set it equal to the design mean
   for a point design prior.
2. Draw the full data sequence up to `max(look_grid)`.
3. Compute statistics from prefixes of that same sequence at each look.
4. Store statistic trajectories, not Bayes factor trajectories as the primary
   artifact.

This preserves the dependence structure across sequential looks. The true
effect must be constant across all looks within one replicate.

### z, Normal Moment, and Directional z Designs

Store:

```text
design_case_id
replicate_id
chunk_id
look
n
se
true_effect
estimate
```

For two-sample standardized mean-difference style cases, generate two full
groups up to `max_n` and compute the cumulative standardized mean difference at
each look. For one-sample or paired analogues, generate one full sequence of
observations or differences.

Analysis cases derive `log_bf01` using:

- `bf01()` for point-null z-test normal analysis priors
- `dirbf01()` for directional z-test priors
- `nmbf01()` for normal moment priors

### t Designs

Store:

```text
design_case_id
replicate_id
chunk_id
look
n
n1
n2
true_effect
t
estimate
```

Analysis cases derive `log_bf01` using `tbf01()`. Include `two.sample`,
`one.sample`, and `paired` designs, plus selected unequal `n1`/`n2` cases.

### Binomial Designs

Store:

```text
design_case_id
replicate_id
chunk_id
look
n
true_p
x
```

Draw one Bernoulli sequence up to `max_n` per replicate. Store cumulative
successes `x` at each look. Analysis cases derive `log_bf01` using
`binbf01()`.

## Design Priors

Use separate design cases for separate design priors.

For continuous-effect families:

- point design prior: `dpsd = 0`, fixed `true_effect = dpm`
- normal design prior: `true_effect ~ N(dpm, dpsd^2)`

For binomial families:

- point design prior: fixed `true_p = dp`
- beta/truncated beta design prior: draw `true_p` once per replicate, then draw
  the full Bernoulli sequence conditional on that value

The design metadata must contain the exact design-prior parameters so analytic
comparisons can call functions such as:

```r
pbf01(..., dpm = design$dpm, dpsd = design$dpsd)
pbf01seq(..., dpm = design$dpm, dpsd = design$dpsd)
ptbf01(..., dpm = design$dpm, dpsd = design$dpsd)
pbinbf01(..., dp = design$dp, da = design$da, db = design$db,
         dl = design$dl, du = design$du)
```

## Analysis Cases

Analysis cases should cover all implemented BF families:

- `bf01()`, `pbf01()`, `nbf01()`, `powerbf01()`
- `dirbf01()` through `pbf01seq(type = "directional")`
- `nmbf01()`, `pnmbf01()`, `nnmbf01()`, `powernmbf01()`
- `tbf01()`, `ptbf01()`, `ntbf01()`, `powertbf01()`
- `binbf01()`, `pbinbf01()`, `nbinbf01()`, `powerbinbf01()`
- `pbf01seq()` and `ptbf01seq()`

Thresholds should include common and adversarial values:

```r
k_h1 <- c(1/30, 1/10, 1/6, 1/3)
k_h0 <- c(3, 6, 10, 30)
```

Sequential cases should include asymmetric threshold pairs, for example
`k1 = 1/30, k0 = 6`, not only reciprocal thresholds.

## Initial Coverage Matrix

Start with a stratified set rather than a full factorial.

### z and Moment Cases

- `psd = 0` point alternative
- local and diffuse normal priors, e.g. `psd = c(0.1, 1, 2)`
- normal moment priors with selected positive `psd`
- `dpm` aligned with analysis prior, equal to the null, opposite sign, and
  misspecified in magnitude
- `dpsd = 0` and positive values such as `0.1`, `0.25`, and `1`
- `usd = c(1, sqrt(2), 2)`

### Directional z Cases

- positive, null, and wrong-direction design priors
- diffuse analysis priors
- cases where H0 or H1 boundaries are hard or impossible at early looks

### t Cases

- `type = c("two.sample", "one.sample", "paired")`
- `alternative = c("two.sided", "greater", "less")`
- default JZS prior and selected informed t priors
- point and uncertain design priors
- unequal `n1`/`n2` for two-sample designs
- known adaptive-root regression cases, including large fractional-n behavior
  represented through deterministic checks and nearby grid looks

### Binomial Cases

- `type = c("point", "direction")`
- `p0 = c(0.2, 0.5, 0.75)`
- flat, skewed, and concentrated beta analysis priors
- point and truncated beta design priors
- cases where power is stepwise and sample-size inversion can oscillate

### Sequential Cases

- one-stage equivalence to fixed-design functions
- comparison schedules with 2 through 20 looks
- stress schedules with 50, 100, and 200 looks, with 100/200 restricted to
  long-grid-compatible designs
- strict package-reference checks use `pbf01seq(strict = TRUE)`; approximate
  sequential references are not correctness oracles for fixture tests
- late-stopping and many-look adversarial cases

## Replications and Precision

Default to:

```r
nsim <- 10000
```

This gives worst-case Monte Carlo standard error around `0.005` for a
probability near `0.5`. Increase to `50000` or `100000` only for selected rare
tail cases or cases where the default precision is not enough to distinguish a
numerical regression from Monte Carlo noise.

Every summary should store Monte Carlo standard errors:

```text
pH1
pH0
pInc
EN
VarN
mcse_pH1
mcse_pH0
mcse_pInc
nsim
```

## Reproducibility and Chunking

Use deterministic chunk-level RNG.

Recommended metadata in `design.rds`:

```text
schema_version
design_case_id
test_family
look_grid_name
look_grid
nsim
chunk_size
rng_kind
master_seed
chunk_seeds
generator_version
bfpwr_version
git_sha
created_at
```

Prefer `L'Ecuyer-CMRG` for parallel reproducibility. Store chunk seeds rather
than one seed per replicate. Each chunk should be independently reproducible:

```text
chunk 1: replicate_id 1-1000, seed S1
chunk 2: replicate_id 1001-2000, seed S2
```

The chunk script writes:

```text
chunks/chunk-0001.rds
chunks/chunk-0002.rds
...
```

The merge step validates:

- all expected chunks exist
- expected replicate ranges are present
- no duplicate `(replicate_id, look)` keys
- every replicate has every look
- `true_effect` or `true_p` is constant across looks within a replicate
- row count equals `nsim * length(look_grid)`
- chunk hashes match the manifest

The chunk files are the canonical stored trajectory artifacts. The finalized
design writes `chunk-index.rds`, `chunk-index.csv`, and `sha256.txt`; a
monolithic `trajectories.rds` may be produced locally for debugging but should
not be part of the GitHub fixture corpus.

## Materialized Analysis Outputs

The primary stored artifact is the design trajectory. Analysis outputs can be
materialized from it when needed:

```text
analysis.rds
outcomes.rds
summary.rds
summary.csv
```

`outcomes.rds` should include:

```text
analysis_case_id
design_case_id
replicate_id
stop_look
stop_n
decision        # H1, H0, inconclusive
final_log_bf01
```

`summary.rds` and `summary.csv` should include:

```text
analysis_case_id
design_case_id
nsim
look_grid_name
k1
k0
pH1
pH0
pInc
EN
VarN
reference_pH1
reference_pH0
reference_pInc
reference_EN
reference_VarN
mcse_pH1
mcse_pH0
mcse_pInc
status
notes
```

`reference_*` values come from deterministic package functions. Optional full
BF trajectories should use `bf_chunks/chunk-*.rds` plus `bf-chunk-index.*`,
not a single monolithic BF trajectory file.

## Validation Rules

Use Monte Carlo standard-error aware tolerances for simulation comparisons:

```r
abs(observed - reference) <= max(abs_tol, z_tol * mc_se)
```

Recommended defaults:

```r
abs_tol <- 1e-8       # deterministic summaries
mc_abs_tol <- 0.002   # floor for MC comparisons
z_tol <- 4            # conservative for fixture tests
```

Use tighter tolerances for deterministic function identity tests:

- closed-form or exact deterministic comparisons: `1e-10` to `1e-8`
- root-finding and optimization checks: `1e-6` to `1e-4`
- t-test adaptive root regressions: keep branch-specific tolerances already
  used by the package tests

For probabilities near 0 or 1, compare the relevant tail directly and keep a
small absolute tolerance floor.

## Package Tests

The corpus-dependent tests should skip unless the corpus is explicitly
available. A helper can use this pattern:

```r
sim_root <- Sys.getenv("BFPWR_SIM_CORPUS")
run_sims <- identical(Sys.getenv("BFPWR_RUN_SIM_FIXTURES"), "true")

if (!run_sims || !nzchar(sim_root) || !dir.exists(sim_root)) {
    tinytest::exit_file("simulation corpus not available")
}
```

CRAN and ordinary package checks should not download or generate simulations.
GitHub Actions or local validation can set:

```text
BFPWR_RUN_SIM_FIXTURES=true
BFPWR_SIM_CORPUS=<repo>/simulations/corpus/v1
```

Tests should usually load summaries first. Full trajectories should be loaded
only for checks that need them, such as:

- recomputing sequential stopping outcomes from `log_bf01`
- querying fixed-design empirical power at selected looks
- checking trajectory invariants
- validating that materialized analysis outputs match design trajectories

## Server Workflow

The intended server workflow is:

```text
Rscript simulations/scripts/run_chunk.R \
  --design-case <design_case_id> \
  --chunk <chunk_id>

Rscript simulations/scripts/merge_chunks.R \
  --design-case <design_case_id>

Rscript simulations/scripts/materialize_analysis.R \
  --analysis-case <analysis_case_id>

Rscript simulations/scripts/validate_case.R \
  --analysis-case <analysis_case_id>
```

For array jobs, create one task per design-case/chunk pair. Analysis
materialization can run after all chunks for a design case have merged.

## Versioning

Use schema and corpus versions:

```text
schema_version: 1
corpus_version: v1
```

Changing column names or semantics requires a new schema version. Adding new
cases under the same schema only updates the manifest.

Case IDs should be stable, lowercase ASCII, and short. Detailed parameters live
in metadata, not only in filenames. Encode decimals with `p` and negatives with
`m`, for example:

```text
z-dnorm-m0p5-s0p1-gridlong-s20260524
z-normal-pm0-psd1-k0p1-10-on-z-dnorm-m0p5-s0p1-gridlong
t-twosample-greater-jzs-dpoint-0p5-gridshort-s20260524
binom-point-p0p5-beta1-1-dpoint0p6-gridshort-s20260524
```

## Implementation Steps

1. Add schema docs, grid helpers, and design/analysis registries.
2. Implement design trajectory generators for z, t, and binomial families.
3. Implement chunk execution and merge validation.
4. Implement BF materialization from stored statistics.
5. Implement summary generation and deterministic reference calculations.
6. Run developer smoke checks locally with small temporary cases.
7. Add gated `tinytest` files that validate summaries and selected trajectory
   queries when `BFPWR_SIM_CORPUS` is set.
8. Dry-run selected registry cases with `nsim`/`chunk_size` overrides before
   spending server time.
9. Inspect file size, runtime, MCSE, and failures; refine case registry.
10. Run the production v1 corpus and use it in GitHub-only validation.

## Open Decisions

- Exact first-pass case registry size and which cases receive `long_grid`.
- Whether generated chunked fixture files should be committed directly or
  attached as release artifacts.
- Whether raw generated data should ever be stored for a small debug subset.
- Preferred server job format for the production run.
