#!/bin/bash
#PBS -q gpuq
#PBS -l nodes=1:ppn=6
#PBS -l walltime=12:00:00
#PBS -V
#PBS -j oe
#PBS -N bridging
######################################
outputfile=__OUTFILE__
bridgingdir=~/miktoarm/SCFT_Mikto/bridging/bridging-scripts
phase=__PHASE__
######################################

cd $PBS_O_WORKDIR
outdir=${PBS_O_WORKDIR}
rundir=${outdir}
username=`whoami`

############# TO USE LOCAL SCRATCH FOR INTERMEDIATE IO, UNCOMMENT THE FOLLOWING
#if [ ! -d /scratch_local/${username} ]; then
#  rundir=/scratch_local/${username}/${PBS_JOBID}
#  mkdir -p $rundir
#  cp ${PSB_O_WORKDIR}/* $rundir
#  cd $rundir
#fi
#####################################################

# Fetch the device ID for the GPU that has been assigned to the job
GPUDEV=`cat $PBS_GPUFILE | awk '{print $1}'`
if [ -z $GPUDEV ]; then
  echo "ERROR finding $PBS_GPUFILE; using default GPU deviceid=0"
  GPUDEV=0
fi

## Prepare the run by substituting the CUDA select device line
## Check whether the line exists first
#grep "CUDA_[Ss]elect[Dd]evice" ${inputfile} > /dev/null
#if [ $? -ne 0 ]; then
#  echo "CUDA_SelectDevice line not found in $inputfile"
#  exit 1
#fi
#sed -i "s/\(CUDA_[Ss]elect[Dd]evice\).*/\1 = ${GPUDEV}/g" ${inputfile}

# Run the job
#${polyftsdir}/PolyFTSGPU.x ${inputfile} > ${outdir}/${outputfile}
EXEC=$bridgingdir/bridging-${phase}.sh
if [ -e $EXEC ]; then
    $EXEC ../ "all" > ${outdir}/${outputfile}
else
    echo "Error $EXEC does not exist!"
fi


# Copy back results
if [ "$rundir" != "$outdir" ]; then
  mv * ${outdir}
fi

# Force good exit code here - e.g., for job dependency
exit 0

