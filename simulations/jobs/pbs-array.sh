#!/bin/sh
# Placeholder PBS array launcher.
#
# Expected environment:
#   BFPWR_SIM_CORPUS  Path to simulations/corpus/v1
#   BFPWR_DESIGN_CASE Design case id
#   PBS_ARRAY_INDEX   Chunk id
#
# Example:
#   qsub -J 1-10 simulations/jobs/pbs-array.sh

set -eu

if [ -z "${BFPWR_SIM_CORPUS:-}" ]; then
  echo "BFPWR_SIM_CORPUS is not set" >&2
  exit 1
fi

if [ -z "${BFPWR_DESIGN_CASE:-}" ]; then
  echo "BFPWR_DESIGN_CASE is not set" >&2
  exit 1
fi

Rscript simulations/scripts/run_chunk.R \
  --corpus-root "$BFPWR_SIM_CORPUS" \
  --design-case "$BFPWR_DESIGN_CASE" \
  --chunk "${PBS_ARRAY_INDEX:-1}"
