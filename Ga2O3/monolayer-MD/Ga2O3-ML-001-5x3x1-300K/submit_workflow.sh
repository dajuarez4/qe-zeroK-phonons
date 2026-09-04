#!/bin/bash
set -euo pipefail
root="$(cd "$(dirname "$0")" && pwd)"
relax_submit=$(cd "$root/01-relax" && sbatch --parsable run_relax.sbatch)
# On federated Slurm systems --parsable may return JOBID;CLUSTER.
relax_id="${relax_submit%%;*}"
if [[ ! "$relax_id" =~ ^[0-9]+$ ]]; then
    echo "ERROR: could not extract a numeric relaxation job ID from: $relax_submit" >&2
    exit 1
fi
md_submit=$(cd "$root/02-md-300K" && sbatch --parsable --dependency="afterok:${relax_id}" run_md.sbatch)
md_id="${md_submit%%;*}"
echo "relax job: $relax_id"
echo "MD job:    $md_id (afterok:$relax_id)"
