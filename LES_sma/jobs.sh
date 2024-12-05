#!/bin/sh 

### General options 
### -- set the job Name -- 
#BSUB -J sma
### -- ask for number of cores (default: 1) -- 
#BSUB -n 32
### -- specify that the cores must be on the same host -- 
#BSUB -R "span[hosts=1]"
### -- set walltime limit: hh:mm -- 
#BSUB -W 12:00 
### -- specify that we need 2GB of memory per core/slot -- 
#BSUB -R "rusage[mem=2GB]"
### -- specify that we want the job to get killed if it exceeds 3 GB per core/slot -- 
#BSUB -M 3GB
### -- set the email address -- 
# please uncomment the following line and put in your e-mail address,
# if you want to receive e-mail notifications on a non-default address
##BSUB -u s232439@dtu.dk
### -- send notification at start -- 
# #BSUB -B 
### -- send notification at completion-- 
#BSUB -N 
### -- Specify the output and error file. %J is the job-id -- 
### -- -o and -e mean append, -oo and -eo mean overwrite -- 
#BSUB -o Output_%J.out 
#BSUB -e Error_%J.err 

# -- load the OpenFOAM module --
module load OpenFoam/v2312/gcc-12.3.0-binutils-2.40-openmpi-4.1.6


# add the path to your case directory:
cd /work3/s232439/s232439-v2312/run/assignment_02/LES_sma

# maybe not working
[ ! -d "0" ] && mkdir 0
# Restore initial conditions
if [ -d "0.orig" ]; then
    cp -r 0.orig/* 0/
else
    echo "Error: Directory '0.orig' not found. Cannot restore initial conditions."
    exit 1
fi

### mesh
blockMesh 
decomposePar 
mpirun -np 32 pimpleFoam -parallel
reconstructPar


### postprocess
touch f.foam

