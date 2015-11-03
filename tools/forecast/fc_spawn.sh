#!/bin/sh


# For the most current run in the TOP directory,
# see if we can start the next cycle
# we can start a new cycle if
#  * the previous run has been archived
#  * or the current run is past 24:00 hours


TOP=/home/jattema/WRF/WRFV3
ARCDIR=/projects/0/sitc/archive2

cd $TOP
RUNDIR=`ls -d $TOP/2015-??-?? | sort | tail -1`

# check if this is acutally a proper directory
if [ ! -d "$RUNDIR" ]; then
    exit
fi

DATESTART=`basename $RUNDIR`
YYYY=`date -d $DATESTART +%Y`
MM=`date -d $DATESTART +%m`
DD=`date -d $DATESTART +%d`

# Check location..
if [ -f "$ARCDIR/${YYYY}/${MM}/${DD}/wrfout_d01_${YYYY}-${MM}-${DD}_00:00:00.nc" ]; then
    # .. in archive
    NCFILE="$ARCDIR/${YYYY}/${MM}/${DD}/wrfout_d01_${YYYY}-${MM}-${DD}_00:00:00.nc"
    echo "Found in archive: $NCFILE"
elif [ -f "$RUNDIR/wrfout_d01_${YYYY}-${MM}-${DD}_00:00:00" ]; then
    # .. in rundir
    NCFILE="$RUNDIR/wrfout_d01_${YYYY}-${MM}-${DD}_00:00:00"
    echo "Found in rundir: $NCFILE"
else
    # not found!
    echo "Cannot find matching netcdf file"
    exit
fi

STEP=`ncdump -h "$NCFILE" | grep UNLIMITED | tr -d "[a-zA-Z=;()/]"`
echo "at step $STEP"
if (( $STEP > 8 )); then
    NEWRUN=`date -d "$DATESTART 24 hours"  +%F`
    echo "Ready for the next cycle, next STARTDATE is $NEWRUN"

    JOB=$TOP/run.${NEWRUN}.job
    cat $TOP/tools/forecast/run.template | sed "s=%RUNDIR%=${TOP}/${NEWRUN}=" | sed "s=%STARTDATE%=$NEWRUN=">"$JOB"
    sbatch "$JOB"
fi
