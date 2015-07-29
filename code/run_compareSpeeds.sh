
for i in $(ls /EcoII/acadia_uni/workspace/simulated/FVCOM/dngridCSR/drifter_runs/BFRIC_0.015/GP/):
do
    python compareSpeeds.py -l GP --dir $i --verbose
done

for x in $(ls /EcoII/acadia_uni/workspace/simulated/FVCOM/dngridCSR/drifter_runs/BFRIC_0.015/DG/):
do
    python compareSpeeds.py -l DG --dir $x --verbose
done

for k in $(ls /EcoII/acadia_uni/workspace/simulated/FVCOM/dngridCSR/drifter_runs/BFRIC_0.015/PP/):
do
    python compareSpeeds.py -l PP --dir $k --verbose
done

echo 'PROGRAM COMPLETE, MY LORD.'
