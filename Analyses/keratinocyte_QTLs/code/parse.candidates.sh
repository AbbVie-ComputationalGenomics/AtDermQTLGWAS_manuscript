echo 'Gene	Stimulus	molecular_trait_id	chromosome	position	ref	alt	variant	ma_samples	maf	pvalue	beta	se	type	ac	an	r2	molecular_trait_object_id	gene_id	median_tpm	rsid' | sed 's/\t/,/g' > data/candidates.stats.csv

for candidate in $(tail -n+2 data/candidates.tsv | sed 's/\t/|/g' | grep -v '\*'); do
  gene=$(echo $candidate | awk -F '|' '{print $1}')
  clump=$(echo $candidate | awk -F '|' '{print $2}')
  rs=$(echo $candidate | awk -F '|' '{print $5}')
  echo $gene $clump $rs
  for file in $(ls data/keratinocyt*gz); do 
	  echo $gene $(echo $file | sed 's/data.*tinocytes_//g' | sed 's/.all.tsv.gz//g') $(tabix $file $clump | grep ":$gene:" | grep -w $rs) | tr ' ' ',' >> data/candidates.stats.csv
  done
done
