#!/bin/bash

#SBATCH -J format_ImmuNexUT_eQTLs
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --mem=20GB
#SBATCH --time=04:00:00
#SBATCH -p magnus
#SBATCH -o /scratch/users/olivamx2/AD/logs/format_ImmuNexUT_eQTLs.%x.%j.o 
#SBATCH -e /scratch/users/olivamx2/AD/logs/format_ImmuNexUT_eQTLs.%x.%j.e

module load htslib

dir=DIR
cell=CELL

cat <(echo "molecular_trait_id	chromosome	position	ref	alt	variant	ma_samples	maf	pvalue	beta	se	type	ac	an	r2	molecular_trait_object_id	gene_id	median_tpm	rsid" > $dir/data/ImmuNexUT_eQTLs/$cell.all.tsv; tail -n+2  $dir/data/ImmuNexUT_eQTLs/$cell'_nominal.txt' | awk '{gsub(/chr/,"",$9);print $1"\t"$9"\t"$10"\t"$7"\t"$8"\tchr"$9"_"$10"_"$7"_"$8"\tNA\tNA\t"$12"\t"$13"\tNA\tNA\tNA\tNA\tNA\tNA\t"$1"\tNA\t"$6 }' | sort --parallel=16 -k2,2n -k3,3n ) > $dir/data/ImmuNexUT_eQTLs/$cell.all.tsv
bgzip -f $dir/data/ImmuNexUT_eQTLs/$cell.all.tsv
tabix -c m -f -s 2 -b 3 -e 3  $dir/data/ImmuNexUT_eQTLs/$cell.all.tsv.gz

