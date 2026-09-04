#!/bin/bash
#SBATCH -J qe-vmsmall-test
#SBATCH -A DMR26002
#SBATCH -p vm-small
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -c 4
#SBATCH -t 00:05:00
#SBATCH -o qe-vmsmall-test.%j.out
#SBATCH -e qe-vmsmall-test.%j.err

set -eo pipefail
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=spread

echo "Started: $(date)"
echo "Host:    $(hostname)"
echo "Job ID:  ${SLURM_JOB_ID}"
echo "Tasks:   ${SLURM_NTASKS}"

module purge
module load intel/19.1.1 impi/19.0.9

QE_BIN="${QE_BIN:-${SCRATCH}/apps/qe-7.5-atomic-fixed-ntyp128/bin/pw.x}"
if [[ ! -x "$QE_BIN" ]]; then
  echo "ERROR: QE executable not found or not executable: $QE_BIN" >&2
  exit 1
fi

echo "QE executable: $QE_BIN"
echo "Testing QE through the Lonestar6 MPI launcher..."
ibrun "$QE_BIN" -help

echo "PASS: private QE 7.5 and ibrun are available."
echo "Finished: $(date)"
