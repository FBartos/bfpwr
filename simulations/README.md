# bfpwr Simulation Fixtures

This directory contains the architecture for the GitHub-only simulation fixture
corpus described in `PLAN.md`.

For a reviewer-facing map of the report artifact, release corpus, and
validation entry points, see `OVERVIEW.md`.

The corpus is split into two layers:

- **Design cases** generate reusable data/statistic trajectories under one
  design prior.
- **Analysis cases** apply one analysis prior and BF threshold pair to an
  existing design case.

No generated fixture corpus is committed by default. The helper functions and
scripts can be smoke-tested against a temporary corpus directory. When the
corpus is generated, design trajectory chunks are the canonical stored
artifacts; finalized designs write `chunk-index.rds`, `chunk-index.csv`, and
`sha256.txt` rather than a monolithic `trajectories.rds`.

The published development corpus is stored as GitHub Release assets on the
fork:

```text
https://github.com/FBartos/bfpwr/releases/tag/sim-corpus-v1
```

Download only the package-test fixture summaries with:

```sh
Rscript simulations/scripts/download_corpus_release.R \
  --corpus-root simulations/corpus/v1 \
  --asset-set fixture-tests
```

Download the full corpus for chunk validation and development with:

```sh
Rscript simulations/scripts/download_corpus_release.R \
  --corpus-root simulations/corpus/v1 \
  --asset-set all
```

Production cases are maintained in `simulations/registry/` using four tiers:
`short`, `long`, `adversarial`, and `regression`. The R registry is the source
of truth; CSV manifests in `simulations/registry/manifests/` are generated for
review.

## Smoke Test

From the repository root:

```sh
Rscript -e "source('simulations/scripts/common.R'); bfpwr_sim_source_package(); smoke <- bfpwr_sim_smoke_run(tempfile('bfpwr-sim-')); print(smoke$summary)"
```

To materialize derived BF trajectories for an analysis case, pass
`--include-bf-trajectories` to `scripts/materialize_analysis.R`. BF trajectories
are written as `bf_chunks/chunk-*.rds` plus `bf-chunk-index.*`; the primary
design artifact remains the reusable statistic trajectory chunks.

## Corpus Tests

Package tests that depend on the full corpus should be gated by environment
variables:

```text
BFPWR_RUN_SIM_FIXTURES=true
BFPWR_SIM_CORPUS=<repo>/simulations/corpus/v1
BFPWR_SIM_REPO=<repo>
```

When those variables are not set, the corpus-dependent tests skip with download
instructions.

Binomial analyses currently have exact fixed-look references through
`pbinbf01()`, but no sequential deterministic reference in the package. Their
sequential summaries therefore record `status = "no_sequential_reference"` and
should be validated through fixed-look queries unless a sequential binomial
function is added later.

## Slow BF Prior Chunks

Some `tbf01()` BF-prior chunks are deterministic but much slower than ordinary
chunks. Use the split workflow for those tasks instead of giving one original
1000-replicate chunk a long walltime. The default split is 25 replicates so
subjobs stay comfortably below the 2-hour PBS worker limit even for the slow
one-sided `t` cases:

```sh
Rscript simulations/scripts/flag_slow_bf_prior_plan.R \
  --plan-index-file simulations/corpus/v1/rerun-bf-prior-plan-tfix-union.csv

Rscript simulations/scripts/split_bf_prior_plan.R \
  --plan-index-file simulations/corpus/v1/rerun-bf-prior-plan-tfix-union.csv \
  --corpus-root simulations/corpus/v1 \
  --skip-present \
  --only-slow \
  --split-size 25 \
  --output-file simulations/corpus/v1/rerun-bf-prior-plan-tfix-union-split.csv \
  --merge-output-file simulations/corpus/v1/rerun-bf-prior-plan-tfix-union-split-merge.csv
```

The split plan writes deterministic partials under
`bf_split_chunks/chunk-XXXX/part-YYYY.rds`. The merge plan combines complete
partials back into the normal `bf_chunks/chunk-XXXX.rds` files, so the existing
finalization script and corpus summaries still apply.
`--skip-present` uses a cheap non-empty-file check to avoid rereading hundreds
of large completed chunks while building recovery plans. Add
`--validate-present` only when explicitly auditing existing chunk contents.

The current slow-task rule flags `t` BF prior rows when any of these hold:

- one-sided `t` BF grid has at least 200,000 BF evaluations;
- one-sided prior direction is opposite the design-prior mean;
- long-tier one-sided `t` BF uses a shifted null.

The normal full-chunk runner now refuses to compute missing rows matching this
rule. Use the split workflow for those rows; override only with
`--allow-slow` or `BFPWR_ALLOW_SLOW_BF_PRIOR=1` when you intentionally want a
single full chunk.

On PBS, submit split rows first, then merge rows, then run the normal BF-prior
finalization jobs:

```sh
SPLIT_PLAN=/auto/brno2/home/fbartos/jobs/bfpwr/simulations/corpus/v1/rerun-bf-prior-plan-tfix-union-split.csv
qsub -J 1-$(($(wc -l < "$SPLIT_PLAN") - 1)) \
  -v BFPWR_SPLIT_PLAN_FILE="$SPLIT_PLAN" \
  simulations/jobs/pbs-bf-prior-split-plan-chunk.sh

MERGE_PLAN=/auto/brno2/home/fbartos/jobs/bfpwr/simulations/corpus/v1/rerun-bf-prior-plan-tfix-union-split-merge.csv
qsub -J 1-$(($(wc -l < "$MERGE_PLAN") - 1)) \
  -v BFPWR_SPLIT_PLAN_FILE="$SPLIT_PLAN",BFPWR_MERGE_PLAN_FILE="$MERGE_PLAN" \
  simulations/jobs/pbs-bf-prior-split-merge.sh
```
