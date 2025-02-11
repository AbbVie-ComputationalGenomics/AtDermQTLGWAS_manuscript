#Activate LDSC and calculate heritability for the AD EUR 4-way meta-analysis

module load miniconda3

conda activate ldsc

cd /data/
wget https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/eur_w_ld_chr.tar.bz2
bunzip2 eur_w_ld_chr.tar.bz2
tar -xvf eur_w_ld_chr.tar 

./ldsc.py \
--h2 /data/AD_EUR.sumstats.gz \
--ref-ld-chr /data/eur_w_ld_chr/ \
--w-ld-chr /data/eur_w_ld_chr/ \
--two-step 1000 \
--samp-prev 0.095 \
--pop-prev 0.15 \
--out /LDSC_Immune_ATAC/AD_EUR_meta_LDSC



#pop-prev 4 cohort EUR meta = 42,963 cases / 451,435 total = 0.095

#Total Liability scale h2: 0.0967 (0.0107)
#Lambda GC: 1.1683
#Mean Chi^2: 1.2925
#Intercept: 1.0496 (0.0123)
#Ratio: 0.1696 (0.0421)