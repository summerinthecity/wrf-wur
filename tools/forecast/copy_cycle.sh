#!/bin/bash

set -e

###################################
#       C  o  n  f  i  g 
###################################

CYCLEFIELDS="TSLB,SMOIS,SH2O,SMCREL,CANWAT,TSK"
TMPFILE=temp.`basename $0`.$$

MANUAL="
`basename $0` copies surface fields ($CYCLEFIELDS) from a WRF output to a WRF input file.
To be used with for a warm start, or as a part of a forecast cycle.

NOTE: 
  This script is unaware of nesting, so run explicitly for each domain.
  No interpolation or re-gridding performed. Make sure the grids match.

Usage:
`basename $0` wrfoutput time wrfinput

  wrfoutput  An outputfile from a WRF run
  time       Time index of field to copy
  wrfinput   wrfinput_d?? file
"


# check arguments
#####################################

if [[ $# != 3 ]]; then
    echo "Got $# arguments"
    echo $MANUAL
    exit -1
fi

if [[ ! -f "$1" ]]; then
    echo "Cannot find input file $1. Aborting"
    exit -1
fi
if [[ ! -f "$3" ]]; then
    echo "Cannot find output file $1. Aborting"
    exit -1
fi


cp "$3" "${TMPFILE}"
ncks -C -A -o "${TMPFILE}" -d Time,$2 -v ${CYCLEFIELDS} "$1"
mv "${TMPFILE}" "$3"

# float TSLB  (Time, soil_layers_stag, south_north, west_east) ; TSLB:description = "SOIL TEMPERATURE" ;
# float SMOIS (Time, soil_layers_stag, south_north, west_east) ; SMOIS:description = "SOIL MOISTURE" ;
# float SH2O  (Time, soil_layers_stag, south_north, west_east) ; SH2O:description = "SOIL LIQUID WATER" ;
# float SMCREL(Time, soil_layers_stag, south_north, west_east) ; SMCREL:description = "RELATIVE SOIL MOISTURE" ;
# float CANWAT(Time,                   south_north, west_east) ; CANWAT:description = "CANOPY WATER" ;
# float TSK   (Time,                   south_north, west_east) ; TSK:description = "SURFACE SKIN TEMPERATURE" ;

# DONT copy these
# ###############

# float TMN       (Time,                   south_north, west_east) ; TMN:description = "SOIL TEMPERATURE AT LOWER BOUNDARY" ;
#    should come from climatological geogrid files (as is WRF's default)

# float TBL_URB_IN(Time, soil_layers_stag, south_north, west_east) ; TBL_URB_IN:description = "WALL LAYER TEMPERATURE" ;
# float TGL_URB_IN(Time, soil_layers_stag, south_north, west_east) ; TGL_URB_IN:description = "ROAD LAYER TEMPERATURE" ;
# float TRL_URB_IN(Time, soil_layers_stag, south_north, west_east) ; TRL_URB_IN:description = "ROOF LAYER TEMPERATURE" ;
#    done in separate script: copy_urb_init
# sst (
#    done in separate script: copy_urb_init
