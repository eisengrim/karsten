NCFILE="/output/subdomain_GP1_0001.nc"
MATFILE="/EcoII/acadia_uni/workspace/observed/GP/Drifter/GP_all.mat"
OUTFILE="/array/home/119865c/test/"

for i in `ls /EcoII/acadia_uni/workspace/simulated/FVCOM/dngridCSR/drifter_runs/BFRIC_0.015/GP/`:
do
    PATH=${f}${NCFILE}
    python -c "import drifter_statistics; calculate_stats($NCFILE, $MATFILE, 'GP', '2017-01-06', plot=True, outpath=$OUTPATH, multi=True)"
done

NCFILE="/array/home/119865c/workspace/acadia_bof_v2_3d_20170106.nc"
MATFILE="/array/home/119865c/workspace/MP_BMP_20170106.nc"
python -c "import drifter_statistics; calculate_stats($NCFILE, $MATFILE, 'MP', '2017-01-06', plot=True, outpath=$OUTPATH, multi=True)"
