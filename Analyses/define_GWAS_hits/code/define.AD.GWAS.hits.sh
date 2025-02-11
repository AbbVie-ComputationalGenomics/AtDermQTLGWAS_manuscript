module load htslib
module load python/python-3.8.2
module load R/R-3.6.1

pheno=$1 # AD_GWAS_meta_FE_DF9_hg38_forMeri
window=$2 # clumping window in Mbs
thres=$3 # GWAS hit threshold

wdir=$PWD
scratch_dir=/ui/abv/olivamx2/scratch-global/AD
AD_GWAS_dir=$wdir/data
  
# Define clumps for each AD GWAS
bash code/plinkClump.2.sh $pheno $window $thres

# Merge clumps into larger overlapping windows
echo AD_GWAS GWAS_hit_index chr start end num_bp num_snps lead_snp_hg19 lead_snp lead_snp_p span > $AD_GWAS_dir/$pheno.AD.GWAS.hits.pvalue"$thres".window"$window"Mb.coords.txt
index_hit=0
for chr in $(seq 1 22); do
	echo $chr
	if [[ -f $AD_GWAS_dir/$pheno.chr"$chr".UKBunrelatedLD.pvalue"$thres".window"$window"Mb.clumped ]]; then
		$wdir/code/parse_clumps.R $window $thres $pheno $chr
		num_loci=$(ls $scratch_dir/$pheno.chr"$chr".UKBunrelatedLD.pvalue"$thres".window"$window"Mb.clumped.loci*.txt | wc -l)
		for locus in $(seq 1 $num_loci); do
			index_hit=$(($index_hit+1))
			num_snps=$(tail -n+2 $scratch_dir/$pheno.chr"$chr".UKBunrelatedLD.pvalue"$thres".window"$window"Mb.clumped.loci"$locus".txt | wc -l)
			start=$(tail -n+2 $scratch_dir/$pheno.chr"$chr".UKBunrelatedLD.pvalue"$thres".window"$window"Mb.clumped.windows.txt | head -n $locus | tail -n 1 | cut -d' ' -f 3)
			end=$(tail -n+2 $scratch_dir/$pheno.chr"$chr".UKBunrelatedLD.pvalue"$thres".window"$window"Mb.clumped.windows.txt | head -n $locus | tail -n 1 | cut -d' ' -f 4)
			num_bp=$((end-start))
			lead_snp_and_pvalue=$(tabix $scratch_dir/$pheno.txt.gz $chr":"$start"-"$end | sort -gk7 |  cut -f 5,6,7 | head -n 1)
			echo $pheno $index_hit $chr $start $end $num_bp $num_snps $lead_snp_and_pvalue >> $AD_GWAS_dir/$pheno.AD.GWAS.hits.pvalue"$thres".window"$window"Mb.coords.txt
			echo $pheno $index_hit $chr $start $end $num_bp $num_snps $lead_snp_and_pvalue 
		done
	fi
done
num_windows=$(grep -v chr $scratch_dir/$pheno.chr*.UKBunrelatedLD.pvalue"$thres".window"$window"Mb.clumped.windows.txt | wc -l)
echo $pheno has $num_windows GWAS hits
