for tissue in $(echo "arm,leg,axilla,back,face,scalp,acral" | tr "," " "); do
  #echo $tissue |  xargs -I {} sbatch --partition=grc  --nodes=1 --cpus-per-task=1 --mem-per-cpu=50g --time=05:00:00 --wrap="./code/perform.Keratinocyte.DE.analysis.R {}"
  #echo $tissue |  xargs -I {} sbatch --partition=grc  --nodes=1 --cpus-per-task=1 --mem-per-cpu=50g --time=05:00:00 --wrap="./code/perform.Keratinocyte_Subtype.DE.analysis.R {}"
  for keratinocyte_subtype in $(echo Basal_Keratinocytes Cycling_Keratinocytes Differentiated_Keratinocytes Keratinized_Keratinocytes); do
	  echo "$tissue $keratinocyte_subtype" |  xargs -I {} sbatch --partition=grc  --nodes=1 --cpus-per-task=1 --mem-per-cpu=50g --time=05:00:00 --wrap="./code/perform.Keratinocyte_Subtype.DE.analysis.R {}"
  done
done

