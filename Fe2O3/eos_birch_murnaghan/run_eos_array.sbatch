#!/bin/bash
#SBATCH --job-name=Fe2O3EOS
#SBATCH --account=jakar_general
#SBATCH --qos=jakar_medium_general
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=38
#SBATCH --exclusive
#SBATCH --array=0-8%3
#SBATCH --output=%A_%a.log

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

POINTS=(V092 V094 V096 V098 V100 V102 V104 V106 V108)
POINT=${POINTS[$SLURM_ARRAY_TASK_ID]}

cd "$SLURM_SUBMIT_DIR/$POINT" || exit 1

mpirun --mca pml ob1 --mca btl ^openib,uct \
  -np $SLURM_NTASKS /shared/quantum-espresso/bin/pw.x \
  -in alpha_Fe2O3.eos.in > alpha_Fe2O3.eos.out

date
