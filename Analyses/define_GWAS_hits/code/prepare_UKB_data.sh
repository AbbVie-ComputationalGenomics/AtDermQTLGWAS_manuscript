wdir=$PWD  # dir=define_GWAS_hits
## UKB data
UKB=$wdir/data
mkdir -p $UKB
cd $UKB

## Gosset

# From Gosset:
#scp /mnt/nvme/users/grundaj/UKBioBank/21574/kinship.txt olivamx2@mghn1.ncsa.illinois.edu:$PWD

## Magnus

# Genotypes
ln -s /projects/abv-ukb/UKBiobank-uncurated/Imputed_GTs/ukb_imp_chr*_v3.bgen* .
# Genotypes info ( "SNP", "rsid", "POS", "Allele1", "Allele2", "MAF", "MinorAllele", "INFO", where INFO is imputation quality. )
ln -s /projects/abv-ukb/UKBiobank-uncurated/Imputed_GTs/ukb_mfi_chr*_v3.txt .
# Sample id
ln -s /projects/abv-ukb/UKBiobank-uncurated/Imputed_GTs/ukb26041_imp_chr1_v3_s487395.sample .
# Subset 15000 unrelated samples
get_seeded_random()
{
  seed="$1"
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt \
    </dev/zero 2>/dev/null
}

#sort -k1,1 ukb26041_imp_chr1_v3_s487395.sample > ukb26041_imp_chr1_v3_s487395.sample.sorted.txt
sort -k1,1 ukb26041_imp_chr1_v3_s487395.sample > ukb26041_imp_chr1_v3_s487395.sample.sorted.txt
join ukb26041_imp_chr1_v3_s487395.sample.sorted.txt <(grep '0$' kinship.txt | shuf --random-source=<(get_seeded_random 42) -n 15000 | cut -f 1 | sort -k1,1) | awk '{print $1"\t"$2}'> 15000_ukb_unrelated_samples.txt

#Subset ancestry-specific individuals
ln -s /projects/grc/finemapping_pipeline/ref_panel_gds_files/ancestry_keys/ukbb/*_id_ukbb.txt . 

for ancestry in $(echo AFR AMR CSA EAS EUR MID); do
	echo $ancestry
	join ukb26041_imp_chr1_v3_s487395.sample.sorted.txt <(grep '0$' kinship.txt | grep -w -F -f "$ancestry"_id_ukbb.txt | shuf --random-source=<(get_seeded_random 42) -n 2500 | cut -f 1 | sort -k1,1) | awk '{print $1"\t"$2}'> 2500_ukb_unrelated_samples.$ancestry.txt
done
cd ..

