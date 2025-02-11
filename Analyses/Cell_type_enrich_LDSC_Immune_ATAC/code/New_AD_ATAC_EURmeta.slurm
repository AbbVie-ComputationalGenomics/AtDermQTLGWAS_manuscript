#!/bin/sh
#SBATCH --time=90:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH --job-name="AD_ATAC_EURmeta"
#SBATCH --output=AD_multi_tissue.log
#SBATCH --mail-type=END

echo "SLURM_JOB_ID="$SLURM_JOB_ID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR

echo "working directory = "$SLURM_SUBMIT_DIR

#Variables to be set
dir=/projects/abv/users/olivamx2/Projects/AD/reproduce_fanying/sLDSC/src/ldsc

cd ${dir}

module load miniconda3

conda activate ldsc

./ldsc.py \
--h2-cts /data/AD_EUR.sumstats.gz \
--ref-ld-chr /data/1000G_Phase3_baseline_v1.2/baseline. \
--out LDSC_Immune_ATAC \
--ref-ld-chr-cts /data/Calderon_Corces_optimal_overlap_ATAC.ldcts \
--w-ld-chr /data/weights_hm3_no_hla/weights.

#
