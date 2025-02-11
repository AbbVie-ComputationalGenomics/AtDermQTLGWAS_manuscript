gcta=$HOME/software/gcta
scratch_dir=/ui/abv/olivamx2/scratch-global/AD
wdir=$PWD

pheno=$1 # AD meta-GWAS e.g. AD_GWAS_meta_FE_DF9_hg38_forMeri  
thres=$2 # e.g. 5e-08 (GWAS hit min pval)
window=$3 # e.g. 2 (Mb)
thres_cond=$4 # e.g. 5e-08 (GWAS hit conditional min pval)

# The GCTA data for this is set up so that the *.ma files contain summary-level information per associated locus
# For more information on GCTA formatting, see here: http://cnsgenomics.com/software/gcta/GCTA_UserManual_v1.24.pdf
# Significance was set at 5e-8

for chr in $(seq 1 22); do
bfile=$wdir/data/ukb_imp_chr"$chr"_v3.15000_ukb_unrelated_samples.filtered_variants.info_score_dot3.maf_dotzerozeroone
	echo $pheno $chr $scratch_dir/$pheno.chr"$chr".UKBunrelatedLD.pvalue"$thres".window"$window"Mb.clumped
	if [[ -f $scratch_dir/$pheno.chr"$chr".UKBunrelatedLD.pvalue"$thres".window"$window"Mb.clumped ]]; then 
		cojo_file=$scratch_dir/$pheno.chr"$chr".UKBunrelatedLD.GCTA.txt
		hits=$(ls $scratch_dir/$pheno.chr"$chr".UKBunrelatedLD.pvalue$thres.window"$window"Mb.clumped.loci* | wc -l)
		for locus in $(seq 1 $hits); do
		    snp_file=$scratch_dir/$pheno.chr"$chr".UKBunrelatedLD.pvalue$thres.window"$window"Mb.clumped.loci$locus.txt
		    # use cojo-slct to identify independently-associated SNPs that reside in the same genomic locus
		    gcta --bfile $bfile \
			--extract $snp_file \
			--chr $chr \
			--cojo-file $cojo_file \
			--cojo-p $thres_cond \
			--cojo-slct \
			--thread-num 4 \
			--out $wdir/data/$pheno.chr"$chr".UKBunrelatedLD.pvalue$thres.window"$window"Mb.clumped.loci$locus

			num_indep=$(tail -n+2 $wdir/data/$pheno.chr"$chr".UKBunrelatedLD.pvalue$thres.window"$window"Mb.clumped.loci$locus.jma.cojo | wc -l);
			if [[ $num_indep -gt 1 ]]; then
				for i in $(tail -n+2 $wdir/data/$pheno.chr"$chr".UKBunrelatedLD.pvalue$thres.window"$window"Mb.clumped.loci$locus.jma.cojo | cut -f 2); do
					grep -v -w $i $wdir/data/$pheno.chr"$chr".UKBunrelatedLD.pvalue$thres.window"$window"Mb.clumped.loci$locus.jma.cojo | tail -n+2  | cut -f 2 > $wdir/data/$pheno.chr"$chr".UKBunrelatedLD.pvalue$thres.window"$window"Mb.clumped.loci$locus.jma.cojo.$i;
			    # use cojo-slct to calculate stats conditioning on independently-associated SNPs that reside in the same genomic locus
				    	gcta --bfile $bfile \
					    --extract $snp_file \
					    --chr $chr \
					    --cojo-file $cojo_file \
					    --cojo-cond $wdir/data/$pheno.chr"$chr".UKBunrelatedLD.pvalue$thres.window"$window"Mb.clumped.loci$locus.jma.cojo.$i \
					    --thread-num 4 \
					    --out $wdir/data/$pheno.chr"$chr".UKBunrelatedLD.pvalue$thres.window"$window"Mb.clumped.loci$locus.$i
				done
			fi
		done
	fi
done
