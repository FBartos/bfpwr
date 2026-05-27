# Simulation Verification Overview

This directory contains the simulation-verification layer for `bfpwr`. It is
kept outside `package/vignettes` because the verification report depends on the
external simulation corpus and should not be part of ordinary package vignette
builds.

## Main Artifacts

- `fixture-verification.Rnw` - source for the simulation-verification report.
- `fixture-verification.pdf` - built report for reviewers and developers.
- `paper-verification.md` - fixed manuscript and sequential-paper value audit.
- `FIXTURE_TESTS.md` - details on the fixture-corpus tinytests and validator
  layers.
- `README.md` - simulation architecture and operational workflow.

## Corpus

The raw simulation corpus is not committed to git. It is published as GitHub
Release assets:

```text
https://github.com/FBartos/bfpwr/releases/tag/sim-corpus-v1
```

Download the fixture summaries needed by package tests:

```sh
Rscript simulations/scripts/download_corpus_release.R \
  --corpus-root simulations/corpus/v1 \
  --asset-set fixture-tests
```

Download the full corpus for chunk-level validation and development:

```sh
Rscript simulations/scripts/download_corpus_release.R \
  --corpus-root simulations/corpus/v1 \
  --asset-set all
```

Then set:

```text
BFPWR_SIM_CORPUS=<repo>/simulations/corpus/v1
BFPWR_SIM_REPO=<repo>
```

## Report Rendering

The report is rendered from the repository root:

```sh
Rscript simulations/scripts/render_fixture_verification_background.R
```

The render script checks for the local corpus first and prints release download
instructions if required files are missing.

## Validation Entry Points

Use ordinary package tests without the corpus env vars; fixture tests skip with
download instructions. Use these commands after downloading the corpus:

```sh
Rscript simulations/scripts/validate_fixture_suite.R \
  --corpus-root simulations/corpus/v1 \
  --fixture-set all \
  --level integrity

Rscript simulations/scripts/validate_z_bf_chunks.R \
  --corpus-root simulations/corpus/v1 \
  --fixture-set z-bf01-fixed-core-v1 \
  --chunk-id 1

Rscript simulations/scripts/validate_binomial_bf_chunks.R \
  --corpus-root simulations/corpus/v1 \
  --fixture-set binom-binbf01-point-fixed-core-v1 \
  --chunk-id 1
```

The package tinytests use the same corpus location through `BFPWR_SIM_CORPUS`
and `BFPWR_SIM_REPO`.
