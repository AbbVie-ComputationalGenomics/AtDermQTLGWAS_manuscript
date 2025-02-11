wdir=$PWD
scratchdir=/scratch/users/olivamx2/AD
mkdir -p $scratchdir/logs

## ImmuNexUT ##
# https://humandbs.biosciencedbc.jp/en/hum0214-v6
mkdir -p $wdir/data/ImmuNexUT_eQTLs
bash code/download.ImmuNexUT_eQTLs.sh

# sample sizes were provided by an ImmuNexUT contact. Cp sample_sizes.txt to $wdir/data/ImmuNexUT_eQTLs/sample_sizes.txt
for cell in $(cut -f 1 $wdir/data/ImmuNexUT_eQTLs/sample_sizes.txt); do
	echo $cell;
	sed "s/CELL/$cell/g" $wdir/code/format.ImmuNexUT_eQTLs.sh | sed "s/DIR/$wdir/g"> $scratchdir/logs/$cell.slurm;
	sbatch $scratchdir/logs/$cell.slurm;
done

## eQTL Catalog ##
# URL
