#!/bin/bash
#SBATCH -J ML001-smoke
#SBATCH -A DMR26002
#SBATCH -p vm-small
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -c 4
#SBATCH -t 00:15:00
#SBATCH -o smoke-test.%j.out
#SBATCH -e smoke-test.%j.err

set -eo pipefail
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=spread

cd "$SLURM_SUBMIT_DIR"
mkdir -p tmp-smoke
module purge
module load intel/19.1.1 impi/19.0.9

QE_BIN="${QE_BIN:-${SCRATCH}/apps/qe-7.5-atomic-fixed-ntyp128/bin/pw.x}"
if [[ ! -x "$QE_BIN" ]]; then
  echo "ERROR: QE executable not found or not executable: $QE_BIN" >&2
  exit 1
fi

echo "QE executable: $QE_BIN"
ibrun "$QE_BIN" -in smoke_test.in > smoke_test.qe.out
grep -q "JOB DONE." smoke_test.qe.out
grep -q "convergence has been achieved" smoke_test.qe.out
echo "PASS: 75-atom 001 SCF smoke test completed successfully."
