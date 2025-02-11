wdir=$PWD
scratchdir=/scratch/users/olivamx2/AD
pushd $scratchdir

wget -r -nH --cut-dirs=2 --no-parent --reject="index.html*" https://ddbj.nig.ac.jp/public/ddbj_database/gea/experiment/E-GEAD-000/E-GEAD-397/ # count data. Unfiltered?
wget -r -nH --cut-dirs=2 --no-parent --reject="index.html*" https://ddbj.nig.ac.jp/public/ddbj_database/gea/experiment/E-GEAD-000/E-GEAD-398/ # eQTLs conditional FDR<0.05
wget -r -nH --cut-dirs=2 --no-parent --reject="index.html*" https://ddbj.nig.ac.jp/public/ddbj_database/gea/experiment/E-GEAD-000/E-GEAD-420/ # eQTLs primary full stats

# Work with eQTLs primary full stats only.
unzip gea/experiment/E-GEAD-000/E-GEAD-420/E-GEAD-420.processed.zip 
bgzip -d nominal_eQTL_2.tar.gz
ln -s $scratchdir/*_nominal.txt $wdir/data/ImmuNexUT_eQTLs 

popd
