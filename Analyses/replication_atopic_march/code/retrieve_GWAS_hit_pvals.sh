# Scan Atopic march phenos in EBI GWAS Catalog reported GWAS hits. Pick suggestive (5e-06) and GWS (5e-08) hits. For many studies, nothing is reported at a suggestive level (even ifit can be found by scanning full summary stats)

for trait_class in $(tail -n+2 GWAS_CAT_traits.tsv | cut -f 1 | uniq); do # for AD/eczema | allergy | asthma (trait class)
> $trait_class.alternative.lite.hits.txt
  for hit in $(cat GWAS_hits.tsv | tr "\t" "!"); do # for each AD GWAS region (AD GWAS hit lead snp +/- 500kb)
    chr=$(echo $hit | cut -d'!' -f 1)
    start=$(echo $hit | cut -d'!' -f 2)
    stop=$(echo $hit | cut -d'!' -f 3)
    lead_snp_rsid=$(echo $hit | cut -d'!' -f 4)
    lead_snp_coords=$(echo $hit | cut -d'!' -f 5)
    lead_snp_p=$(echo $hit | cut -d'!' -f 6)
    ancestry=$(echo $hit | cut -d'!' -f 7)
    nearest_gene=$(echo $hit | cut -d'!' -f 8)
    #ontology_terms=$(grep $trait_class GWAS_CAT_traits.tsv | cut -f 3 | tr "\n" "|" | sed "s/|$//g") # grep mapping ontology terms
    #echo $ontology_terms
    # for that region, grep P<1e-06 GWAS hits mapping to trait class corresponding ontology terms
    echo $(tabix data/alternative.lite.gz $chr":"$start"-"$stop | grep -f <(grep $trait_class GWAS_CAT_traits.tsv | cut -f 3) | awk -F "\t" '{if($28 <= 5e-06 && 5e-08 < $28) {print "Reported:"$37":Suggestive:"$22"_"$28} if(5e-08 >= $28) {print "Reported:"$37":GWS:"$22"_"$28}}' | tr "\n" "|") >> $trait_class.alternative.lite.hits.txt 
    #awk -v ontoterms=$ontology_terms '{gsub(".*\\/","",$36);if($36 ~ ontoterms) {print $0}}' #DISCARDED
  done > $trait_class.alternative.lite.hits.txt
done

# grep -f <(tail -n+2 GWAS_CAT_traits.tsv | cut -f 3 | uniq) studies | grep -vf <(cut -f 3 GWAS_CAT_studies.tsv) | cut -f 1,3,6,7,13,14,15 | grep 2024 > 2024_studies.txt # DISREGARDED: keep largest studies (below)

# Use largest asthma study (GCST010043) to retrieve the smallest pvalue variant per AD GWAS hit and assign GWAS significance category
for region in $(awk '{print $1":"$2"-"$3}' GWAS_hits.hg19.tsv | head -n 101); do   line=$(tabix data/HanY_prePMID_asthma_Meta-analysis_UKBB_TAGC.txt.gz  $region | sort -gk7,7 | head -n 1 | tr "\t" "|");   snp=$(echo $line | cut -d'|' -f1);   p=$(echo $line | cut -d'|' -f7);   category=$(echo $p | awk '{category="NotSignif";if($1 <= 5e-06 && 5e-08 < $1) {category="Suggestive"} if(5e-08 >= $1) {category="GWS"}; print category}');   echo "FullSumStats:GCST010043:"$category":"$snp"_"$p; done > asthma.fullsummarystats.hits.txt 

# Use largest AD studies (GCST90244788 EUR and GCST90244787 MULTI)
for region in $(awk '{print $1":"$2"-"$3}' GWAS_hits.hg19.tsv | head -n 101); do     line=$(cat <(tabix data/GCST90244787_buildGRCh37.tsv.gz  $region | cut -f 1,2; tabix data/GCST90244788_buildGRCh37.tsv.gz  $region | cut -f 1,2) | sort -gk2,2 | head -n 1 |  tr "\t" "|");     snp=$(echo $line | cut -d'|' -f1);     p=$(echo $line | cut -d'|' -f2);   category=$(echo $p | awk '{category="NotSignif";if($1 <= 5e-06 && 5e-08 < $1) {category="Suggestive"} if(5e-08 >= $1) {category="GWS"}; print category}');   echo "FullSumStats:GCST90244787/GCST90244788:"$category":"$snp"_"$p; done > ADeczema.fullsummarystats.hits.txt

# Use largest allergic study (GCST005038) 
for region in $(awk '{print $1":"$2"-"$3}' GWAS_hits.hg19.tsv | head -n 101); do   line=$(tabix data/SHARE-without23andMe.LDSCORE-GC.SE-META.v0.gz  $region | sort -gk8,8 | head -n 1 | tr "\t" "|");   snp=$(echo $line | cut -d'|' -f1);   p=$(echo $line | cut -d'|' -f8);   category=$(echo $p | awk '{category="NotSignif";if($1 <= 5e-06 && 5e-08 < $1) {category="Suggestive"} if(5e-08 >= $1) {category="GWS"}; print category}');   echo "FullSumStats:GCST005038:"$category":"$snp"_"$p; done > allergy.fullsummarystats.hits.txt

paste <(paste ADeczema.alternative.lite.hits.txt ADeczema.fullsummarystats.hits.txt | tr "\t" "|") <(paste allergy.alternative.lite.hits.txt allergy.fullsummarystats.hits.txt | tr "\t" "|") <(paste asthma.alternative.lite.hits.txt asthma.fullsummarystats.hits.txt | tr "\t" "|") > SupplementaryTable11.txt
