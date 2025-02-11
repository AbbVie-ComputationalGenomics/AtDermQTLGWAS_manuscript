#!/bin/sh

pheno=$1 # AD meta-GWAS e.g. AD_GWAS_meta_FE_DF9_hg38_forMeri
window=$2 # window. e.g. 2 (2000kb, 2Mb)
thres=$3 # GWAS hit p-value threshols e.g. 5e-08
scratch_dir=/ui/abv/olivamx2/scratch-global/AD
wdir=$PWD
ukb_dir=$wdir/data
AD_GWAS_dir=$wdir/data

module load htslib
# Prepare AD GWAS
echo "Prepare $pheno GWAS data"
echo "1) Sort $pheno GWAS data"
ln -s /projects/abv/GRC/rileybm/projects/Atopic_Dermatitis/GWAMA_FE_6_cohort_meta_FG_R9/$pheno.txt.gz $AD_GWAS_dir
#cp file of sample sizes to data/AD_GWAS_FE_sample_count.v2.txt

# Format we want:
# The columns rsID and n_cases were added by 28th Apr 2022. n_cases is an estimation in many cases.
# CHR	POS	REF	ALT	MARKER	PVAL	BETA	SE	Z	ALT_AF	n_studies	n_samples	hg19_CHR	hg19_POS	rsID	n_cases

# Format of file provided by Bridget:

# CHR POS REF ALT MARKER rsID PVAL BETA SE Z ALT_AF n_studies n_samples effects hg19_CHR hg19_POS
#1 10177 A AC 1:10177:A:AC rs367896724 0.344458 -0.019829 0.020975 -0.945377 0.40052 3 232105 ?-??++ 1 10177
echo 'CHR	POS	REF	ALT	MARKER	rsID	PVAL	BETA	SE	Z	ALT_AF	n_studies	n_samples	effects	n_cases' > $AD_GWAS_dir/$pheno.txt
awk 'FNR==NR{a[$1]=$3;next}{if (a[$13]) {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,a[$13]} else {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,"NA"}}' <(tail -n+2 $AD_GWAS_dir/AD_GWAS_FE_sample_count.v2.txt) <(zcat $AD_GWAS_dir/$pheno.txt.gz | tail -n+2) | tr ' ' '\t' | sort --parallel 8 -nk1,1 -nk2,2n >> $AD_GWAS_dir/$pheno.txt
ln -s $AD_GWAS_dir/$pheno.txt $scratch_dir

echo "2) Compress $pheno GWAS data"
bgzip -f $scratch_dir/$pheno.txt
echo "3) tabix $pheno GWAS data" 
tabix -c C -f -s 1 -b 2 -e 2 $scratch_dir/$pheno.txt.gz
echo "Done"

# make sure plink really knows how many threads it can use! this might be system-specific variable
export OMP_NUM_THREADS=8

# No clumps for chr13, 15, 18 and 21, apparently
for chr in {1..22}; do
echo chr $chr
# GWAS data SHOULD BE formatted in the same way as is expected by GCTA: SNP A1 A2 freq b se p N
gwas_gcta_file=$scratch_dir/$pheno.chr$chr.UKBunrelatedLD.GCTA.txt

echo "SNP A1 A2 freq b se p N" > $gwas_gcta_file
# A1 corresponds to effect allele, agnostic to genome build
tabix $scratch_dir/$pheno.txt.gz $chr | \
awk '{print $6"\t"$4"\t"$3"\t"$11"\t"$8"\t"$9"\t"$7"\t"$13}' | \
tr "\t" " " >> $gwas_gcta_file

#CHR POS REF ALT MARKER PVAL BETA SE Z ALT_AF n_studies n_samples hg19_CHR hg19_POS rsID
#1 10177 A AC 1:10177:A:AC 0.344458 -0.019829 0.020975 -0.945377 0.40052 3 232105 1 10177 rs367896724

# run the clumping setting genome-wide significance at p = 5e-8 and pruning down to an r2 = 0.05, and allowing p-values down to 0.05
# see the Plink 1.9 documentation for full details

module load plink/plink
clumpfile=$AD_GWAS_dir/$pheno.chr$chr.UKBunrelatedLD.pvalue"$thres".window"$window"Mb
echo $clumpfile

tmpdir=/ui/abv/olivamx2/scratch-global/UKB 
plink --bfile $tmpdir/ukb_imp_chr"$chr"_v3.15000_ukb_unrelated_samples.filtered_variants.info_score_dot3.maf_dotzerozeroone \
    --clump $gwas_gcta_file \
    --clump-field p \
    --clump-snp-field SNP \
    --clump-p1 $thres \
    --clump-kb "$window"000 \
    --clump-r2 0.05 \
    --clump-p2 0.05 \
    --out $clumpfile \
    --memory 16000 \
    --threads 1

done
