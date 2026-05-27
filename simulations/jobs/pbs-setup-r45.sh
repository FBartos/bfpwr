#!/bin/sh
#PBS -N bfpwr-setup-r45
#PBS -l select=1:ncpus=1:mem=4gb:scratch_local=1gb
#PBS -l walltime=01:00:00
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

BFPWR_R_MODULE="${BFPWR_R_MODULE:-r/4.5.1-gcc-10.2.1-zmneq6c}"
BFPWR_R_LIBS_USER="${BFPWR_R_LIBS_USER:-/auto/brno2/home/fbartos/Rpackages-45}"
export BFPWR_PROJECT_DIR BFPWR_R_MODULE BFPWR_R_LIBS_USER

MODULE_INIT="/cvmfs/software.metacentrum.cz/modulefiles/5.3.1/libexec/modulecmd.tcl"
MODULE_TCL="/cvmfs/software.metacentrum.cz/modulefiles/5.3.1/bin/tclsh"

eval "$("$MODULE_TCL" "$MODULE_INIT" sh autoinit)"
module load "$BFPWR_R_MODULE"
mkdir -p "$BFPWR_R_LIBS_USER"
export R_LIBS_USER="$BFPWR_R_LIBS_USER"

cd "$BFPWR_PROJECT_DIR"

Rscript simulations/scripts/setup_cluster_r45.R \
  --lib "$BFPWR_R_LIBS_USER" \
  --repo-root "$BFPWR_PROJECT_DIR"

Rscript simulations/scripts/record_environment.R \
  --output-file "$BFPWR_PROJECT_DIR/simulations/corpus/v1/environment/setup-r45-${PBS_JOBID:-manual}.rds"
