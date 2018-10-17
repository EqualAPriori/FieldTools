#!/bin/bash
#PBS -q batch
#PBS -l nodes=1:ppn=6
#PBS -l walltime=24:00:00
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

cat $PBS_NODEFILE > nodes

# Run the job
#EXEC=$bridgingdir/bridging-${phase}.sh
EXEC=$bridgingdir/bridging-wrapper-FK.sh
if [ -e $EXEC ]; then
    #$EXEC ${phase} ../ "all" > ${outdir}/${outputfile}
    $EXEC ${phase} ../ > ${outdir}/${outputfile}
else
    echo "Error $EXEC does not exist!"
fi

# Copy back results
if [ "$rundir" != "$outdir" ]; then
  mv * ${outdir}
fi

# Force good exit code here - e.g., for job dependency
exit 0
