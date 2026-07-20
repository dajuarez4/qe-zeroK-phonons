#!/bin/bash
#SBATCH --job-name=AlSymPh
#SBATCH --account=jakar_partner_physj
#SBATCH --qos=jakar_medium_physj
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=38
#SBATCH --exclusive
#SBATCH --output=%j.log

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

QE=/shared/quantum-espresso/bin
MPI_ARGS="--mca pml ob1 --mca btl ^openib,uct"

mpirun $MPI_ARGS -np $SLURM_NTASKS ${QE}/ph.x \
  -in alpha_Al2O3.ph_grid_2x2x2.in > alpha_Al2O3.ph_grid_2x2x2.out
mpirun $MPI_ARGS -np 1 ${QE}/q2r.x \
  -in alpha_Al2O3.q2r.in > alpha_Al2O3.q2r.out

for index in 1 2 3 4; do
  mpirun $MPI_ARGS -np 1 ${QE}/matdyn.x \
    -in "alpha_Al2O3.path${index}.matdyn.in" \
    > "alpha_Al2O3.path${index}.matdyn.out"
done

mpirun $MPI_ARGS -np 1 ${QE}/matdyn.x \
  -in alpha_Al2O3.phdos.in > alpha_Al2O3.phdos.out

python3 plot_alpha_al2o3_phonons.py

date
