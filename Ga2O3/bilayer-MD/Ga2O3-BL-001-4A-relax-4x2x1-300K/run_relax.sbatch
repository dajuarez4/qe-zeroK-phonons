#!/bin/bash
#SBATCH --job-name=Ga2O3-001-4A-relax
#SBATCH --account=jakar_partner_physj
#SBATCH --qos=jakar_medium_physj
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=38
#SBATCH --exclusive
#SBATCH --output=%j.log

set -eo pipefail

cd "$SLURM_SUBMIT_DIR"
mkdir -p tmp/relax

pwd
hostname
date
ulimit -s unlimited

# Do not use `set -u`: Intel setvars references optional unset variables.
source /opt/intel/oneapi/setvars.sh -ofi_internal=1 --force
module purge
module load gnu12
module load openmpi4/4.1.5
module load tbb
module load compiler-rt
module load mkl

mpirun --mca pml ob1 --mca btl '^openib,uct' -np "$SLURM_NTASKS" \
  /shared/quantum-espresso/bin/pw.x \
  -in Ga2O3-BL-001.relax.in > Ga2O3-BL-001.relax.out

grep -q "JOB DONE." Ga2O3-BL-001.relax.out
grep -q "End of BFGS Geometry Optimization" Ga2O3-BL-001.relax.out
if grep -q "maximum number of steps has been reached" Ga2O3-BL-001.relax.out; then
  echo "ERROR: relaxation reached nstep without force convergence" >&2
  exit 1
fi

echo "Relaxation converged. Add Ga2O3-BL-001.relax.out to the preparation folder."
echo "Do not start MD until the q-grid/TDEP thermal configuration is generated."
date
