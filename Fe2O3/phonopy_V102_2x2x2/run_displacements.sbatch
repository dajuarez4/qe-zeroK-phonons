#!/bin/bash
#SBATCH --job-name=Fe2O3-V102-FD
#SBATCH --account=jakar_partner_physj
#SBATCH --qos=jakar_medium_physj
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=38
#SBATCH --exclusive
#SBATCH --array=0-5%2
#SBATCH --output=%A_%a.log

set -euo pipefail

pwd
hostname
date
ulimit -s unlimited

# Intel's initialization script references optional variables that may be
# unset. Temporarily disable Bash nounset, then restore strict checking.
set +u
source /opt/intel/oneapi/setvars.sh -ofi_internal=1 --force
set -u
module purge
module load gnu12
module load openmpi4/4.1.5
module load tbb
module load compiler-rt
module load mkl

# Avoid threaded MKL inside each MPI rank. QE's distributed ScaLAPACK
# diagonalizer crashes with this QE 7.0/OpenMPI/MKL combination, so the
# pw.x launch below uses serial subspace diagonalization (-ndiag 1).
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

DISPLACEMENTS=(disp-001 disp-002 disp-003 disp-004 disp-005 disp-006)
POINT=${DISPLACEMENTS[$SLURM_ARRAY_TASK_ID]}

cd "$SLURM_SUBMIT_DIR/displacements/$POINT"

mpirun --mca pml ob1 --mca btl ^openib,uct \
  -np "$SLURM_NTASKS" /shared/quantum-espresso/bin/pw.x \
  -npool 2 -ndiag 1 \
  -in alpha_Fe2O3.fd.scf.in > alpha_Fe2O3.fd.scf.out

grep -q 'JOB DONE' alpha_Fe2O3.fd.scf.out
date
