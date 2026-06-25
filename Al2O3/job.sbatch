#!/bin/bash
#SBATCH --job-name=test
#SBATCH --account=jakar_partner_physj
#SBATCH --qos=jakar_medium_physj
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --output=%j.log
#SBATCH --error=%j.err

pwd
hostname
date
ulimit -s unlimited

module purge
module load gnu12/12.2.0
module load tbb/2021.13
module load mkl/2024.2

export OMP_NUM_THREADS=1

# OpenMPI root directory
export MPI_ROOT=/opt/ohpc/pub/mpi/openmpi4-gnu12/4.1.5
export PATH=${MPI_ROOT}/bin:$PATH
export LD_LIBRARY_PATH=${MPI_ROOT}/lib:$LD_LIBRARY_PATH

# Intel MKL libraries
export LD_LIBRARY_PATH=/opt/intel/oneapi/mkl/2024.2/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/opt/intel/oneapi/2024.2/lib:$LD_LIBRARY_PATH

QE=/shared/quantum-espresso/bin/ph.x
INPUT=Al2O3.phG.in
OUTPUT=Al2O3.phG.out

echo "Using mpirun:"
which mpirun
mpirun --version

echo "Checking missing libraries:"
ldd ${QE} | grep "not found" || echo "No missing libraries"

echo "Running ph.x..."
mpirun -np ${SLURM_NTASKS} ${QE} -in ${INPUT} > ${OUTPUT}

date
