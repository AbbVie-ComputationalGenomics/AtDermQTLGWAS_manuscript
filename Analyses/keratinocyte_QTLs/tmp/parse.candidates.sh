for candidate in $(cat tmp/candidates.tsv); do
  gene=$(echo $candidate | awk -F '.' '{print $1}')
  chr=$(echo $candidate | awk -F '.' '{print $2}')
  start=$(echo $candidate | awk -F '.' '{print $3}')
  end=$(echo $candidate | awk -F '.' '{print $4}')
  rs=$(echo $candidate | awk -F '.' '{print $5}')
  echo $gene $chr':'$start'-'$end $rs
  for file in $(ls data/keratinocyt*gz); do 
  echo $gene $rs $file $(tabix $file $chr':'$start'-'$end | grep ":$gene:" | grep -w $rs)
  done
done