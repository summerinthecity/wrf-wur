#!/bin/bash

#################################################################
# Config

# all errors are fatal
set -e

# CDO and NCO paths
if [ `hostname | grep -i cartesius` ]; then
    module load cdo
    module load nco
fi
if [ `hostname | grep maq06` ]; then
    export PATH=/home/escience/local/bin:$PATH
fi

# main data file
strt=~/SST/watertemp.nc

#-----------------------------------------------------------------


#################################################################
# Command line arguments

function usage {
    echo "Interpolate Rijkswaterstaat observations to a WRF grid"
    echo ""
    echo "Usage: $0 date griddes output"
    echo ""
    echo " - griddes is a NetCDF file containing the required grid, or a grid"
    echo "   description. Grid descriptions can be made using 'cdo griddes wrfout.nc',"
    echo "   but be sure to select only the relevant grid from that output!"
    echo ""
    echo " - date can be anything like YYYY-MM-DD or 'yesterday'"
    echo ""
    exit -1
}

# date for SST
WHEN=$1

# griddescriptor
GRID=$2

# output
OUTPUT=$3

if [ x$WHEN = x ]; then
   usage
fi

if [ x$GRID = x ]; then
   usage
fi

if [ x$OUTPUT = x ]; then
   usage
fi

#-----------------------------------------------------------------




#################################################################
# Interpolate

# tempfiles
tempa=${strt}.nobnds
tempb=${strt}.daymax
tempc=${strt}.tint
tempd=${strt}.remap
tempe=${strt}.last
tempf=${strt}.kelvin

function clean {
    rm -f ${tempa} ${tempb} ${tempc} ${tempd} ${tempe} ${tempf}
}



# start clean
clean

# remove time_bnds to prevent confusion in cdo
ncks -x -v time_bnds -O -o ${tempa} ${strt}

# find day maxima
cdo daymax ${tempa} ${tempb}

# interpolate over time
cdo inttime,2014-06-01,00:00,1days ${tempb} ${tempc}

# inverse-distance weighted nearest neighbour
cdo remapdis,${GRID} ${tempc} ${tempd} 

# take last time step
ncks -d time,`date -d $WHEN +'%F'` -O -o ${tempe} ${tempd}

# another fill missing pass on the spatial
# cdo fillmiss,4 ${tempe} sst_d0${i}.nc
cdo remapdis,${GRID} ${tempe} ${tempf}

cdo addc,273.15 ${tempf} $OUTPUT

# clean up
clean

#-----------------------------------------------------------------
