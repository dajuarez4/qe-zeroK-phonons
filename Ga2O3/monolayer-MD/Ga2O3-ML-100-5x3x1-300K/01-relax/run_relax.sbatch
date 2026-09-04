#!/bin/bash
#SBATCH -J ML100-relax
#SBATCH -A DMR26002
#SBATCH -p vm-small
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -c 4
#SBATCH -t 48:00:00
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err

set -eo pipefail
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=spread

cd "$SLURM_SUBMIT_DIR"
mkdir -p tmp
module purge
module load intel/19.1.1 impi/19.0.9

# QE 7.5 installation already verified in the user's Lonestar6 environment.
QE_BIN="${QE_BIN:-${SCRATCH}/apps/qe-7.5-atomic-fixed-ntyp128/bin/pw.x}"
if [[ ! -x "$QE_BIN" ]]; then
  echo "ERROR: QE executable not found or not executable: $QE_BIN" >&2
  exit 1
fi

ibrun "$QE_BIN" -in relax.in > relax.out
grep -q "JOB DONE." relax.out
grep -q "End of BFGS Geometry Optimization" relax.out
