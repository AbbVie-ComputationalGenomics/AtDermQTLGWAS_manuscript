pheno=AD_GWAS_12meta_FE_hg38_forMeri
window=$1 # 1 (Mb)
thres=$2 # 5e-08
thres_cojo=$3 # 5e-08

# Define AD GWAS clumps, given UKB LD data, for trans-ancestry and ancestry-stratified AD GWASs
bash code/define.AD.GWAS.hits.sh $pheno $window $thres # located in define_GWAS_hits folder
for ancestry in $(echo ASN AFR AMR EUR); do bash code/define.AD.GWAS.hits.ancestry.sh $pheno $window $thres $ancestry; done

# Identify independent signals per AD GWAS loci, for trans-ancestry and ancestry-stratified AD GWASs (ASN AFR AMR EUR)
for chr in $(seq 1 22); do sed "s/\$chr/$chr/g" code/gctaSelect.template.qsub | sed "s/\$window/$window/g" |  sed "s/\$thres/$thres/g" | sed "s/\$pheno/$pheno/g" | sed "s/\$thres_cojo/$thres_cojo/g" > tmp/gctaSelect.template.chr$chr.window"$window".thres"$thres".thres_cojo"$thres_cojo"qsub; bash tmp/gctaSelect.template.chr$chr.window"$window".thres"$thres".thres_cojo"$thres_cojo"qsub;done
for ancestry in $(echo ASN AFR AMR EUR); do for chr in $(seq 1 22); do sed "s/\$chr/$chr/g" code/gctaSelect.ancestry.template.qsub | sed "s/\$ancestry/$ancestry/g" | sed "s/\$window/$window/g" |  sed "s/\$thres/$thres/g" | sed "s/\$pheno/$pheno/g" | sed "s/\$thres_cojo/$thres_cojo/g" > tmp/gctaSelect.$ancestry.template.chr$chr.window"$window".thres"$thres".thres_cojo"$thres_cojo"qsub; sbatch tmp/gctaSelect.$ancestry.template.chr$chr.window"$window".thres"$thres".thres_cojo"$thres_cojo"qsub;done;done

# Retrieve annotations of plink and GCTA conditional signal, for trans-ancestry and ancestry-stratified AD GWASs that have GWAS hit clumps (ASN EUR)
paste  <(tail -n+2 data/AD_GWAS_12meta_FE_hg38_forMeri.AD.GWAS.hits.pvalue5e-08.window1Mb.coords.txt) <(bedtools merge -c 4 -o collapse,count -i <(tail -n+2 /ui/abv/olivamx2/scratch-global/AD/AD_GWAS_12meta_*1Mb*raw* | grep rs |  sort -nk1,1 -k3,3n | awk '{print $1"\t"$3"\t"$4"\t"$2}') | cut -f 4,5) <(for chr in $(seq 1 22); do bash code/gctaSelect.parse.sh AD_GWAS_12meta_FE_hg38_forMeri 5e-08 1 5e-08 $chr all; done)| tr ' ' '\t'
paste  <(tail -n+2 data/AD_GWAS_EUR_FE_hg38_forMeri.AD.GWAS.hits.pvalue5e-08.window1Mb.coords.txt) <(bedtools merge -c 4 -o collapse,count -i <(tail -n+2 /ui/abv/olivamx2/scratch-global/AD/AD_GWAS_EUR_*1Mb*raw* | grep rs |  sort -nk1,1 -k3,3n | awk '{print $1"\t"$3"\t"$4"\t"$2}') | cut -f 4,5) <(for chr in $(seq 1 22); do bash code/gctaSelect.parse.sh AD_GWAS_EUR_FE_hg38_forMeri 5e-08 1 5e-08 $chr all; done)| tr ' ' '\t'
paste  <(tail -n+2 data/AD_GWAS_ASN_FE_hg38_forMeri.AD.GWAS.hits.pvalue5e-08.window1Mb.coords.txt) <(bedtools merge -c 4 -o collapse,count -i <(tail -n+2 /ui/abv/olivamx2/scratch-global/AD/AD_GWAS_ASN_*1Mb*raw* | grep rs |  sort -nk1,1 -k3,3n | awk '{print $1"\t"$3"\t"$4"\t"$2}') | cut -f 4,5) <(for chr in $(seq 1 22); do bash code/gctaSelect.parse.sh AD_GWAS_ASN_FE_hg38_forMeri 5e-08 1 5e-08 $chr all; done)| tr ' ' '\t'

# Define non-redundant set of AD GWAS loci combining trans-ancestry and ancestry-stratified AD GWAS loci
bash code/define.AD.GWAS.hits.ancestry.merged.sh $pheno $window $thres $thres_cojo

# Merge GWAS+QTL signal and perform colocalization
bash code/merge_gwas_hits_QTLs.and.coloc.sh pQTLs pvalue
bash code/merge_gwas_hits_QTLs.and.coloc.sh eQTLs beta
bash code/merge_gwas_hits_QTLs.and.coloc.sh ImmuNexUT_eQTLs pvalue
bash code/merge_gwas_hits_QTLs.and.coloc.sh mQTLs beta
bash code/merge_gwas_hits_QTLs.and.coloc.sh sQTLs/exon beta
bash code/merge_gwas_hits_QTLs.and.coloc.sh sQTLs/txrev beta
bash code/merge_gwas_hits_QTLs.and.coloc.sh sQTLs/tx beta

# Summarize results
bash summarize.results.sh
