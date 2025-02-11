for tissue in $(echo "arm,leg,axilla,back,face,scalp,acral" | tr "," " "); do
  echo $tissue |  xargs -I {} sbatch --partition=grc  --nodes=1 --cpus-per-task=1 --mem-per-cpu=50g --time=05:00:00 --wrap="./code/read.umichigan.seurat.stratify.per.bodysite.R {}"
done
