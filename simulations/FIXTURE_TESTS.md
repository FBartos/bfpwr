# Fixture Corpus Tests

The large simulation corpus is a GitHub-only test fixture. Package checks must
not require it by default. Fixture tests should therefore run only when the user
explicitly points the test process at a local corpus.

The current corpus is published as release assets at:

```text
https://github.com/FBartos/bfpwr/releases/tag/sim-corpus-v1
```

For package fixture tests, download only the fixture summaries:

```sh
Rscript simulations/scripts/download_corpus_release.R \
  --corpus-root simulations/corpus/v1 \
  --asset-set fixture-tests
```

For cached chunk validation, download the design and BF-prior chunks:

```sh
Rscript simulations/scripts/download_corpus_release.R \
  --corpus-root simulations/corpus/v1 \
  --asset-set chunk-validation
```

## Activation

Set these environment variables before running the fixture tests:

```sh
BFPWR_SIM_CORPUS=/path/to/bfpwr/simulations/corpus/v1
BFPWR_SIM_REPO=/path/to/bfpwr
```

If `BFPWR_SIM_CORPUS` is unset or points at an incomplete corpus, fixture tests
exit early with the release download command. `BFPWR_SIM_REPO` is optional when
the corpus is stored under the repository's `simulations/corpus/v1` path, but
setting it is clearer for CI.

## Test Layers

There are two separate validation layers.

Simulation-side validation checks the generated corpus itself:

1. **Chunk integrity**
   Verify that expected corpus chunk files exist and conform to the simulation
   schema. This catches broken downloads, incomplete finalization, and registry
   drift.

2. **Cached BF trajectory equivalence**
   Recompute selected BF trajectories from reusable design statistic chunks and
   compare them to cached BF chunks. This is the main regression check for
   `bf01()`, `dirbf01()`, `nmbf01()`, `tbf01()`, and `binbf01()`.

Package-side fixture tests consume the published fixture summaries:

1. **Fixture integrity**
   Verify that every file listed in a fixture summary manifest exists and
   matches its recorded byte size and SHA256 hash.

2. **Fixed-look MC versus package references**
   For selected fixed looks, compare empirical probabilities from fixture
   summaries against deterministic package reference functions, allowing MC
   error. The z-test checks use `pbf01()`, `pnmbf01()`, and one-look
   `pbf01seq(type = "directional")` references.

3. **Sequential decision summaries**
   Query cached BF trajectories for first crossing of `k1` and `k0`, and compare
   cumulative and final stopping summaries with `pbf01seq(strict = TRUE)` where
   z-test sequential references are computationally viable.

4. **Search-function checks**
   Use fixture-backed fixed-look summaries to exercise sample-size and
   threshold-search helpers at settings where closed forms are available.
   Sequential search targets are intentionally empty until the sequential search
   API is implemented.

## Initial Coverage

The active z fixed-summary package fixture tests are:

- `package/inst/tinytest/test-fixture-corpus-z-bf01-fixed-summary.R`
- `package/inst/tinytest/test-fixture-corpus-z-nmbf01-fixed-summary.R`
- `package/inst/tinytest/test-fixture-corpus-z-dirbf01-fixed-summary.R`

The active z sequential-summary package fixture tests are:

- `package/inst/tinytest/test-fixture-corpus-z-bf01-sequential-summary.R`
- `package/inst/tinytest/test-fixture-corpus-z-nmbf01-sequential-summary.R`
- `package/inst/tinytest/test-fixture-corpus-z-dirbf01-sequential-summary.R`

The active binomial fixed-summary package fixture tests are:

- `package/inst/tinytest/test-fixture-corpus-binomial-binbf01-point-fixed-summary.R`
- `package/inst/tinytest/test-fixture-corpus-binomial-binbf01-direction-fixed-summary.R`

The active binomial sequential-summary package fixture tests are:

- `package/inst/tinytest/test-fixture-corpus-binomial-binbf01-point-sequential-summary.R`
- `package/inst/tinytest/test-fixture-corpus-binomial-binbf01-direction-sequential-summary.R`

The active t fixed and sequential package fixture tests are:

- `package/inst/tinytest/test-fixture-corpus-t-tbf01-fixed-summary.R`
- `package/inst/tinytest/test-fixture-corpus-t-tbf01-sequential-summary.R`

The z fixed-summary tests cover deterministic fixed-design summary rows across:

- BF thresholds 3, 10, and 30
- H1 and H0 evidence tails
- short and long look grids
- null, alternative, and adversarial designs
- point and normal analysis priors

The binomial fixed tests cover point-null and directional `binbf01()`
summaries across:

- BF thresholds 3, 10, and 30
- H1 and H0 evidence tails
- short and long look grids
- null, alternative, and adversarial designs
- default beta, matched-null, concentrated beta, truncated-directional, and
  large-shape analysis-prior settings

The binomial fixed tests compare selected rows against `pbinbf01()` with
Monte Carlo-aware tolerances. The binomial sequential fixtures are materialized
for future package-function validation, but there is currently no
package-level sequential binomial oracle. The active package tests therefore
only validate manifest integrity, deterministic probability/count identities,
registered schedule coverage, and selected summary-row presence. Raw
cached-trajectory recomputation belongs in simulation-side validation rather
than package fixture tests until a sequential binomial package function exists.

The t fixed tests compare stratified rows against `ptbf01()` across thresholds
3/10/30, H1/H0 tails, short and long grids, null/alternative/adversarial
designs, one-sample/paired/two-sample/unequal two-sample paths, Cauchy and
Student-t priors, shifted nulls, and one-sided/two-sided alternatives. The t
sequential tests validate manifest and deterministic summary integrity for all
registered schedules, including the 50/100/200-look stress schedules, and
compare selected one-sided 20-look summaries against
`ptbf01seq(strict = TRUE)`. Shifted-null unequal two-sample and paired
sequential t settings are currently kept as structural checks because the
strict package oracle shows known approximation gaps at early looks.

The sequential z tests verify all registered schedule IDs are present. They
compare stratified 5-look cumulative/final summaries against
`pbf01seq(strict = TRUE)` for `bf01()`, `nmbf01()`, and `dirbf01()`. The
`bf01()` and `dirbf01()` tests also include explicit schedule-shape reference
checks at 2, 20, and 50 looks, plus curated strict long-run checks at 100 and
200 looks. Exact moment-prior sequential references are much more expensive, so
`nmbf01()` package fixture tests stop at the 5-look exact-reference layer while
the 50/100/200-look moment summaries remain covered by manifest and
deterministic-summary validation.

The raw cached `z-bf01` trajectory equivalence check lives in the simulation
tooling instead:

```sh
Rscript simulations/scripts/validate_z_bf_chunks.R \
  --corpus-root simulations/corpus/v1
```

By default it covers one representative cached chunk for each z-test BF family:

- `bf01()` with a point analysis prior under a standard null design
- `nmbf01()` with a moment prior under a standard null design
- `dirbf01()` with a directional prior under a standard null design

The simulation-side check validates chunk presence, schema validity, key
stability, and cached `log_bf01` equality against current package
recomputation for chunk 1. It can also validate all cached BF priors selected by
a fixture set, for example:

```sh
Rscript simulations/scripts/validate_z_bf_chunks.R \
  --corpus-root simulations/corpus/v1 \
  --fixture-set z-nmbf01-fixed-core-v1 \
  --chunk-id 1 \
  --bf-type moment

Rscript simulations/scripts/validate_z_bf_chunks.R \
  --corpus-root simulations/corpus/v1 \
  --fixture-set z-dirbf01-fixed-core-v1 \
  --chunk-id 1 \
  --bf-type directional
```

The corresponding binomial cached-BF check is:

```sh
Rscript simulations/scripts/validate_binomial_bf_chunks.R \
  --corpus-root simulations/corpus/v1

Rscript simulations/scripts/validate_binomial_bf_chunks.R \
  --corpus-root simulations/corpus/v1 \
  --fixture-set binom-binbf01-point-fixed-core-v1 \
  --chunk-id 1 \
  --bf-type point

Rscript simulations/scripts/validate_binomial_bf_chunks.R \
  --corpus-root simulations/corpus/v1 \
  --fixture-set binom-binbf01-direction-fixed-core-v1 \
  --chunk-id 1 \
  --bf-type direction
```

The binomial raw-cache validator uses a default `log_bf01` tolerance of
`1e-10`, because large-shape beta stress cases can differ across R/platform
builds by a few `1e-12` from floating-point roundoff in large beta-function
calculations.

## Expansion Order

Add new fixture tests in this order:

1. `z-bf01` fixed-look MC versus `pbf01()`.
2. `z-dirbf01` cached trajectories.
3. `z-nmbf01` cached trajectories.
4. `tbf01` cached trajectories after the long-running rerun and finalization
   finish.

Keep each test file focused on one BF family/type so failures identify the
affected implementation quickly.

## Derived Fixture Summaries

The corpus also supports precomputed fixture summaries under
`simulations/corpus/v1/fixtures/<family>/<mode>/<fixture_set_id>/`.

Fixed summaries contain empirical BF tail probabilities, symmetric decision
probabilities, and first-grid sample-size search results. Sequential summaries
contain cumulative stopping probabilities, final stopping/sample-size summaries,
and, once enabled, first-look search results for registered schedules. The
current z sequential summaries write a typed zero-row `search-summary.rds`.
The z sequential schedules use a comparison grid with 2 through 20 looks and
three stress schedules with 50, 100, and 200 looks. The 100- and 200-look
stress schedules are long-grid only.

The raw BF chunks remain the canonical data. Derived summaries are regenerated
with:

```sh
Rscript simulations/scripts/materialize_fixed_fixture.R \
  --corpus-root simulations/corpus/v1 \
  --fixture-set z-bf01-fixed-core-v1

Rscript simulations/scripts/materialize_fixed_fixture.R \
  --corpus-root simulations/corpus/v1 \
  --fixture-set z-nmbf01-fixed-core-v1

Rscript simulations/scripts/materialize_fixed_fixture.R \
  --corpus-root simulations/corpus/v1 \
  --fixture-set z-dirbf01-fixed-core-v1

Rscript simulations/scripts/materialize_sequential_fixture.R \
  --corpus-root simulations/corpus/v1 \
  --fixture-set z-bf01-sequential-core-v1

Rscript simulations/scripts/materialize_sequential_fixture.R \
  --corpus-root simulations/corpus/v1 \
  --fixture-set z-nmbf01-sequential-core-v1

Rscript simulations/scripts/materialize_sequential_fixture.R \
  --corpus-root simulations/corpus/v1 \
  --fixture-set z-dirbf01-sequential-core-v1

Rscript simulations/scripts/materialize_fixed_fixture.R \
  --corpus-root simulations/corpus/v1 \
  --fixture-set binom-binbf01-point-fixed-core-v1

Rscript simulations/scripts/materialize_fixed_fixture.R \
  --corpus-root simulations/corpus/v1 \
  --fixture-set binom-binbf01-direction-fixed-core-v1

Rscript simulations/scripts/materialize_sequential_fixture.R \
  --corpus-root simulations/corpus/v1 \
  --fixture-set binom-binbf01-point-sequential-core-v1

Rscript simulations/scripts/materialize_sequential_fixture.R \
  --corpus-root simulations/corpus/v1 \
  --fixture-set binom-binbf01-direction-sequential-core-v1

Rscript simulations/scripts/materialize_fixed_fixture.R \
  --corpus-root simulations/corpus/v1 \
  --fixture-set t-tbf01-fixed-core-v1 \
  --all-registry-priors

Rscript simulations/scripts/materialize_sequential_fixture.R \
  --corpus-root simulations/corpus/v1 \
  --fixture-set t-tbf01-sequential-core-v1 \
  --all-registry-priors
```

The initial z fixture specs use symmetric evidence thresholds 3, 10, and 30.
Fixed-design search targets include 10%, 30%, 80%, 90%, and 95% for both H1 and
H0 evidence. Sequential-design search targets are intentionally empty for now.

## Exhaustive Pre-PR Validation

The pre-PR fixture suite is intentionally stronger than the package tinytests.
It validates every materialized fixture summary row, hard-fails deterministic
invariants, and uses distributional MC diagnostics for agreement with package
reference functions:

```sh
Rscript simulations/scripts/validate_fixture_suite.R \
  --corpus-root simulations/corpus/v1 \
  --fixture-set all \
  --level exhaustive \
  --output-dir simulations/corpus/v1/fixture-validation
```

Use `--require-all` when every registered fixture set is expected to be
materialized. Without it, `--fixture-set all` validates the materialized subset
and warns about missing fixture directories.

Deterministic failures include manifest byte/hash mismatches, duplicate keys,
invalid probabilities/counts/MCSEs, impossible log-BF status counts, broken
fixed or sequential search summaries, and summary tables that disagree with
their `overview.csv`.

Package-reference agreement is evaluated on all fixed-look rows where a
deterministic reference function is available. Sequential z summaries are
checked against a curated set of strict `pbf01seq()` references covering
thresholds 3/10/30 on short schedules and BF 10 on 50/100/200-look stress
schedules. Moment-prior sequential references are intentionally limited to 2-
and 5-look schedules, and uncertain normal-prior references are kept to short
schedules, because exact many-look references are slow. Sequential t summaries
are checked against selected strict `ptbf01seq()` references at 20 looks for
one-sided greater/less cases covering ordinary alternatives, one-sample local
effects, and wrong-direction adversarial designs; additional shifted-null,
paired, unequal, and 50/100/200-look t schedules remain covered structurally.

Use `--sequential-reference-profile none` to skip the strict sequential package
oracle while retaining deterministic fixture validation.

The validator reports Holm-corrected exact binomial row tests as diagnostics,
but hard-fails package-reference agreement only for severe binomial prediction
interval misses. This avoids treating expected many-comparison MC fluctuations
as implementation failures while still flagging practically incompatible rows.
It also reports aggregate residual diagnostics, including mean/sd of MC
z-residuals, coverage within 1/2/3 MCSE, outlier counts, and grouped residual
summaries.
