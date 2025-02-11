module load htslib
module load python/python-3.8.2
module load R/R-3.6.1

pheno=$1 # AD_GWAS_12meta_FE_hg38_forMeri
window=$2 # clumping window in Mbs (e.g. 1)
thres=$3 # GWAS hit threshold (e.g. 5e-08)

tmpdir=/ui/abv/olivamx2/scratch-global/UKB # NOTE: customize and create the scratchdir accordingly: e.g. mkdir /ui/abv/olivamx2/scratch-global/UKB. Modify subsequent scripts specifying the tmpdir of choice
wdir=$PWD

# Step 1
# Prepare UKB data for different analysis, i.e. retrieve genotypes, select sample sets based on relatedness, ancestry, etc.
# NOTE: kinship.txt has been copied directly, as gosset is out of comission. The symbolic links rely on unchanged original data alocation
bash $wdir/code/prepare_UKB_data.sh

scratch_dir=/ui/abv/olivamx2/scratch-global/AD 

# Step 2
# Select UKB variants and samples to feed LD info to the clumping process. Agnostic to the AD GWAS, but needed for clumping and other steps.
# NOTE: create the scratchdir accordingly
for chr in $(seq 1 22); do sed "s/\$chr/$chr/g" code/plinkClump.template.qsub > tmp/plinkClump.template.chr$chr.qsub; sbatch tmp/plinkClump.template.chr$chr.qsub;done
#sed "s/\$casesfile/$AD_GWAS_12meta_nsamples_ncases/g" code/format.GWAS.sh | sed "s/\$pheno/$pheno/g"> tmp/format.GWAS.window"$window".thres"$thres".qsub; sbatch tmp/format.GWAS.window"$window".thres"$thres".qsub;

# Step 3
# Define clumps for each AD GWAS
# Run only when code/plinkClump.template.qsub finishes
for chr in $(seq 1 22); do sed "s/\$chr/$chr/g" code/plinkClump.2.template.qsub | sed "s/\$window/$window/g" |  sed "s/\$thres/$thres/g" | sed "s/\$pheno/$pheno/g" > tmp/plinkClump.2.template.chr$chr.window"$window".thres"$thres".qsub; sbatch tmp/plinkClump.2.template.chr$chr.window"$window".thres"$thres".qsub;done

# Merge clumps into larger overlapping windows
# Run only when code/plinkClump.2.template.qsub finishes
echo AD_GWAS GWAS_hit_index chr start end num_bp num_snps lead_snp_hg19 lead_snp lead_snp_p> data/$pheno.AD.GWAS.hits.pvalue"$thres".window"$window"Mb.coords.txt

index_hit=0
for chr in $(seq 1 22); do
        echo $chr
        if [[ -f $scratch_dir/$pheno.chr"$chr".UKBunrelatedLD.pvalue"$thres".window"$window"Mb.clumped ]]; then
                ./code/parse_clumps.R $window $thres $pheno $chr
                num_loci=$(ls $scratch_dir/$pheno.chr"$chr".UKBunrelatedLD.pvalue"$thres".window"$window"Mb.clumped.loci*.txt | wc -l)
                for locus in $(seq 1 $num_loci); do
                        index_hit=$(($index_hit+1))
                        num_snps=$(tail -n+2 $scratch_dir/$pheno.chr"$chr".UKBunrelatedLD.pvalue"$thres".window"$window"Mb.clumped.loci"$locus".txt | wc -l)
                        start=$(tail -n+2 $scratch_dir/$pheno.chr"$chr".UKBunrelatedLD.pvalue"$thres".window"$window"Mb.clumped.windows.txt | head -n $locus | tail -n 1 | cut -d' ' -f 3)
                        end=$(tail -n+2 $scratch_dir/$pheno.chr"$chr".UKBunrelatedLD.pvalue"$thres".window"$window"Mb.clumped.windows.txt | head -n $locus | tail -n 1 | cut -d' ' -f 4)
                        num_bp=$((end-start))
                        lead_snp_and_pvalue=$(tabix $scratch_dir/$pheno.txt.gz $chr":"$start"-"$end | sort -gk7 |  cut -f 5,6,7 | head -n 1)
                        echo $pheno $index_hit $chr $start $end $num_bp $num_snps $lead_snp_and_pvalue >> data/$pheno.AD.GWAS.hits.pvalue"$thres".window"$window"Mb.coords.txt
                        echo $pheno $index_hit $chr $start $end $num_bp $num_snps $lead_snp_and_pvalue 
                done
        fi
done
num_windows=$(grep -v chr $scratch_dir/$pheno.chr*.UKBunrelatedLD.pvalue"$thres".window"$window"Mb.clumped.windows.txt | wc -l)
echo $pheno has $num_windows GWAS hits

## Disregard code below. Kept for historical reasons



module load R/R-3.6.3

tmpdir=/ui/abv/olivamx2/scratch-global/AD # NOTE: customize and create the scratchdir accordingly: e.g. mkdir /ui/abv/olivamx2/scratch-global/AD. Modify subsequent scripts specifying the tmpdir of choice
pheno=AD_GWAS_meta_FE_DF9_hg38_forMeri
window=$1 # 1 (Mb)
thres=$2 # 5e-08
thres_cojo=$3 # 5e-08

# Define AD GWAS clumps, given UKB LD data
bash code/define.AD.GWAS.hits.sh $pheno $window $thres

# Identify independent signals per AD GWAS loci
bash code/gctaSelect.sh $pheno $thres $window $thres_cojo
