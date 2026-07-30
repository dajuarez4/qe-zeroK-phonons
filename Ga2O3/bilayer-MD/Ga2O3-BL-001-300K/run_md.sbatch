#!/bin/bash
#SBATCH --job-name=Ga2O3-001-MD
#SBATCH --account=jakar_partner_physj
#SBATCH --qos=jakar_medium_physj
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=38
#SBATCH --exclusive
#SBATCH --output=%j.log

# Intel oneAPI's setvars.sh references optional variables that may be unset,
# so nounset (-u) cannot be enabled before sourcing the cluster environment.
set -eo pipefail

cd "$SLURM_SUBMIT_DIR"
mkdir -p tmp

pwd
hostname
date
ulimit -s unlimited

source /opt/intel/oneapi/setvars.sh -ofi_internal=1 --force
module purge
module load gnu12
module load openmpi4/4.1.5
module load tbb
module load compiler-rt
module load mkl

mpirun --mca pml ob1 --mca btl ^openib,uct \
  -np "$SLURM_NTASKS" /shared/quantum-espresso/bin/pw.x \
  -in Ga2O3-BL-001.md.in > Ga2O3-BL-001.md.out

date
