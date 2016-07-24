## MESH is name of mesh file. FIELD is the field config file.
## Ensure interpolated velocity files for u, v, and w are written from
## runInterpAvg.py. Then run parseNETCDF4.py for the actual files. 
## Field data for foods may be generated using generateFood.py.
## Use the tag -? on any function for help and info.
## To see the run time of multi-threaded processes, use the command:
##       ps axo pid,etime,args
## Gradients must be computed for z, depth and all fields. 
~/HP-IBM/prepmesh -vb MESH.ini
~/HP-IBM/convert_to_dat -vf MESH.ini
~/HP-IBM/calcdepth -vi 0 -o 5 MESH.ini
~/HP-IBM/fielddata -vo MESH.ini FIELD.ini
~/HP-IBM/precomputegrads -vi 5 MESH.ini
~/HP-IBM/floaters -vj 4 MESH.ini
~/HP-IBM/vtkexport -vt 480 MESH.ini
~/HP-IBM/vtkparticles -vf MESH.ini
