#!/bin/sh
#PBS -N bfpwr-bf-finalize
#PBS -l select=1:ncpus=1:mem=4gb:scratch_local=2gb
#PBS -l walltime=02:00:00
#PBS -j oe

set -eu

cleanup_scratch() {
  if command -v clean_scratch >/dev/null 2>&1; then
    clean_scratch
  fi
}
trap cleanup_scratch TERM EXIT

if [ -z "${BFPWR_PROJECT_DIR:-}" ]; then
  BFPWR_PROJECT_DIR="/auto/brno2/home/fbartos/jobs/bfpwr"
fi

if [ -z "${BFPWR_SIM_CORPUS:-}" ]; then
  BFPWR_SIM_CORPUS="$BFPWR_PROJECT_DIR/simulations/corpus/v1"
fi

PLAN_INDEX="${BFPWR_PLAN_INDEX:-${PBS_ARRAY_INDEX:-1}}"
BFPWR_R_MODULE="${BFPWR_R_MODULE:-r/4.5.1-gcc-10.2.1-zmneq6c}"
BFPWR_R_LIBS_USER="${BFPWR_R_LIBS_USER:-/auto/brno2/home/fbartos/Rpackages-45}"
export BFPWR_R_MODULE BFPWR_R_LIBS_USER

MODULE_INIT="/cvmfs/software.metacentrum.cz/modulefiles/5.3.1/libexec/modulecmd.tcl"
MODULE_TCL="/cvmfs/software.metacentrum.cz/modulefiles/5.3.1/bin/tclsh"

eval "$("$MODULE_TCL" "$MODULE_INIT" sh autoinit)"
module load "$BFPWR_R_MODULE"
mkdir -p "$BFPWR_R_LIBS_USER"
export R_LIBS_USER="$BFPWR_R_LIBS_USER"

RUN_DIR="${SCRATCHDIR:-/tmp}/${USER:-fbartos}/bfpwr-bf-finalize-${PBS_JOBID:-manual}-${PLAN_INDEX}"
mkdir -p "$RUN_DIR"
cp -a "$BFPWR_PROJECT_DIR/package" "$RUN_DIR/package"
mkdir -p "$RUN_DIR/simulations"
cp -a "$BFPWR_PROJECT_DIR/simulations/R" "$RUN_DIR/simulations/R"
cp -a "$BFPWR_PROJECT_DIR/simulations/scripts" "$RUN_DIR/simulations/scripts"
cp -a "$BFPWR_PROJECT_DIR/simulations/registry" "$RUN_DIR/simulations/registry"

cd "$RUN_DIR"
if [ -n "${BFPWR_PLAN_INDEX_FILE:-}" ]; then
  Rscript simulations/scripts/finalize_bf_prior_plan.R \
    --corpus-root "$BFPWR_SIM_CORPUS" \
    --plan-index "$PLAN_INDEX" \
    --plan-index-file "$BFPWR_PLAN_INDEX_FILE"
else
  Rscript simulations/scripts/finalize_bf_prior_plan.R \
    --corpus-root "$BFPWR_SIM_CORPUS" \
    --plan-index "$PLAN_INDEX"
fi

echo "completed BF prior finalization index $PLAN_INDEX"
