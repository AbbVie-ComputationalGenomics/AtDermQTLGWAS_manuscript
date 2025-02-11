#!/bin/sh
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH --job-name="Miami_plot"
#SBATCH --output=miamiplot.log
#SBATCH --mail-type=END

cd /projects/abv/GRC/rileybm/projects/AD_metagwas_paper/Figures/fig.GWAS_miami_plot_comparisons

# Load software
module load r/4.2.0

# Run R script
Rscript --vanilla code/Miami_combo_script.R