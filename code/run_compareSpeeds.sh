
for i in `ls /EcoII/acadia_uni/workspace/simulated/FVCOM/dngridCSR/drifter_runs/BFRIC_0.015/GP/`:
do
    python compareSpeeds.py -l GP --dir $i --verbose
    python compareSpeeds.py -l GP --dir $i --verbose --ratio 0.9348 
done

for x in `ls /EcoII/acadia_uni/workspace/simulated/FVCOM/dngridCSR/drifter_runs/BFRIC_0.015/DG/`:
do
    python compareSpeeds.py -l DG --dir $x --verbose
    python compareSpeeds.py -l DG --dir $x --verbose --ratio 1.227
done

for k in `ls /EcoII/acadia_uni/workspace/simulated/FVCOM/dngridCSR/drifter_runs/BFRIC_0.015/PP/`:
do
    python compareSpeeds.py -l PP --dir $k --verbose
    python compareSpeeds.py -l PP --dir $k --verbose --ratio 0.9377
done


echo 'PROGRAM COMPLETE, MY LORD.'
