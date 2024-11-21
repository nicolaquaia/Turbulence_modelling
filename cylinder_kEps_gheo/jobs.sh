#!/bin/sh 
### General options 
### -- set the job Name -- 
#BSUB -J assignment2
### -- ask for number of cores (default: 1) -- 
#BSUB -n 8 
### -- specify that the cores must be on the same host -- 
#BSUB -R "span[hosts=1]"
### -- set walltime limit: hh:mm -- 
#BSUB -W 02:00 
### -- specify that we need 2GB of memory per core/slot -- 
#BSUB -R "rusage[mem=2GB]"
### -- specify that we want the job to get killed if it exceeds 3 GB per core/slot -- 
#BSUB -M 3GB
### -- set the email address -- 
# please uncomment the following line and put in your e-mail address,
# if you want to receive e-mail notifications on a non-default address
##BSUB -u s233119@dtu.dk
### -- send notification at start -- 
#BSUB -B 
### -- send notification at completion-- 
#BSUB -N 
### -- Specify the output and error file. %J is the job-id -- 
### -- -o and -e mean append, -oo and -eo mean overwrite -- 
#BSUB -o Output_%J.out 
#BSUB -e Error_%J.err 

# -- load the OpenFOAM module --
module load OpenFoam/v2312/gcc-12.3.0-binutils-2.40-openmpi-4.1.6


# add the path to your case directory:
cd /work3/s233119/s233119-v2312/run/cylinder_kEps

canCompile || exit 0    # Dynamic code

# fix: file in folder 0.orig and not 0
mkdir 0
restore0Dir

### mesh
blockMesh 
mirrorMesh
# fix: mirrowed mesh in 0.0001 and not constant
rm -rf constant/polyMesh
cp -r 0.0001/polyMesh constant

### solver
pimpleFoam


### postprocess
touch f.foam

# add openfoam commands:
# blockMesh
# rm -r 0
# cp -r 0.orig 0
# decomposePar  -force 


# mpirun pimpleFoam -parallel > Output.log


# Extra tips:
# To submit the jobscript, type: bsub < jobscript , then you get a job ID to track the status, kill the job etc.
# To track the job (pending, running, done, exiting): bstat (no ID needed)
# To kill a job: bkill jobID
# 
# Read more about submitting jobs to DTU HPC: https://www.hpc.dtu.dk/?page_id=2868#Parallel
# This file is modified for 46115 Turbulence Modeling - E24
