#!/bin/bash
#SBATCH --job-name=Al69-FD
#SBATCH --account=jakar_general
#SBATCH --qos=jakar_medium_general
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=38
#SBATCH --exclusive
#SBATCH --output=%A_%a.log

set -euo pipefail

if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
  echo "Submit this as the array printed by generate_displacements.sh."
  exit 1
fi

pwd
hostname
date
ulimit -s unlimited

set +u
source /opt/intel/oneapi/setvars.sh -ofi_internal=1 --force
set -u
module purge
module load gnu12
module load openmpi4/4.1.5
module load tbb
module load compiler-rt
module load mkl

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

mapfile -t points < <(find "$SLURM_SUBMIT_DIR/displacements" -mindepth 1 -maxdepth 1 -type d -name 'disp-*' | sort)
if [ "$SLURM_ARRAY_TASK_ID" -ge "${#points[@]}" ]; then
  echo "Array index $SLURM_ARRAY_TASK_ID exceeds ${#points[@]} displacement folders."
  exit 1
fi
cd "${points[$SLURM_ARRAY_TASK_ID]}"

mpirun --mca pml ob1 --mca btl ^openib,uct \
  -np "$SLURM_NTASKS" /shared/quantum-espresso/bin/pw.x \
  -npool 2 -ndiag 1 \
  -in Al69_Ga2O3.fd.scf.in > Al69_Ga2O3.fd.scf.out

grep -q 'JOB DONE' Al69_Ga2O3.fd.scf.out
date
