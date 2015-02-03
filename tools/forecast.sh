#!/bin/bash

###########################3
# Forecast config:

# location of configuration file
CONFIG=/home/jiska/forecast.config
ARCHIVE=/home/jiska/archive            # Must exist


# location of external tools
NCDUMP=ncdump
NC3TONC4=nc3tonc4



##################################
# FIXME read from config file
NDOMS=5
DATESTART="2012-09-06"
##################################


MANUAL="
Control WRF forecast runs.

   $0 command option

prepare:
  all
  date         Set the datetime for the run
  boundaries   Download boundaries from NCEP and run ungrib, metgrid
  urban        Copy urban fields from previous run
  sst          Set lake temperatures

run:
  all
  real         Run real.exe
  wrf          Queue the run for execution

zip:
  all 
  ts           Zip TS files
  netcdf       Convert netCDF output to netCDF-4 format with compression
  log          Zip WRF log files

archive:
  all
  ts           Copy TS files to archive dir
  netcdf       Copy netCDF to archive dir
  log          Copy log file to archive dir

clean:
  all 
  input        Remove WRF boundaries
  output       Remove all WRF output files

status         Print forecast status
"


function log {
    message=$1

    stamp=`date -u +"[%s] %F %T"`
    printf '%s : %s\n' "$stamp" "$message" >> "$CONFIG"
}

function help {
    printf '%s\n' "$MANUAL"
    exit -1
}

function status {
    echo "Status!"
    log "yes..."
}

# Remove all input files from a WRF run,
# after a 'clean input', a 'prepare all' must be run
#
# Assumes:
#    CWD is the 'wrf' run directory
function clean_input {

    # starting state and boundaries
    # -----------------------------

    rm -f wrfinput*
    rm -r wrfbdy*    
}

#
# Remove all output files from a WRF run,
# after a 'clean output' the WRF forecast can be resubmitted
#
# Assumes:
#    CWD is the 'wrf' run directory
function clean_output {
    FORCE=$1

    # log files
    # -------------------------

    rm -f rsl.out.*
    rm -f rs.error.*
    rm -f logs.zip

    # wrfout
    # ------------------------

    rm -f wrfout*

    # time series
    # -----------

    if [ -f tslist ]; then
        STATIONS=`cat tslist | awk '{print $2}'`

        for d in `seq -f '%02.0f' 1 $NDOMS`; do
            for s in $STATIONS; do
                rm -f ${s}.d${d}.*
            done
        done
    else
        printf "$0 [$LINENO]: Can't open tslist file, time series not removed\n"
    fi
}

# Modifies the namelist.input to contain the required date
# Thed date can be given as function argument, or is set to one day after the last forecast
#
# Assumes:
#    CWD is the 'wrf' run directory
# Required env:
#    NDOMS
function prepare_date {
    if [[ -z "$NDOMS" ]]; then
        printf "$0 [$LINENO]: NDOMS not set, aborting\n"
        exit 1;
    fi

    DATESTART=$1
    if [ ! "$DATESTART" ]; then
        printf "$0 [$LINENO]: No valid date given, aborting\n"
        exit 1;
    fi

    YEAR=`echo $DATESTART | sed 's/-/ /g' | awk '{print $1}'`
    MONTH=`echo $DATESTART | sed 's/-/ /g' | awk '{print $2}'`
    DAY=`echo $DATESTART | sed 's/-/ /g' | awk '{print $3}'`

    # create a new namelist
    cat namelist.input | \
    sed -e "s/^.*START_YEAR.*$/START_YEAR = $NDOMS*$YEAR/"    \
        -e "s/^.*START_MONTH.*$/START_MONTH = $NDOMS*$MONTH/" \
        -e "s/^.*START_DAY.*$/START_DAY = $NDOMS*$DAY/"    >  \
    namelist.new || (
        printf "$0 [$LINENO]: Unable to create namelist, aborting\n"
        exit 1
    )

    mv namelist.new namelist.input || (
        printf "$0 [$LINENO]: Unable to write namelist, aborting\n"
        exit 1
    )
}

# Zip WRF timeseries files
# Assumes:
#    CWD is the wrf 'run' directory    
#    tslist is in the CWD
# Required env:
#    NDOMS
function zip_ts {
    if [[ -z "$NDOMS" ]]; then
        printf "$0 [$LINENO]: NDOMS not set, aborting\n"
        exit 1;
    fi

    if [ -f tslist ]; then
        STATIONS=`cat tslist | awk '{print $2}'`
    else
        printf "$0 [$LINENO]: Can't open tslist file, aborting\n"
        exit 1
    fi

    for d in `seq -f '%02.0f' 1 $NDOMS`; do
        for s in $STATIONS; do

            # if there are ts files for this station and domain to process,
            if [ `ls ${s}.d${d}.* > /dev/null` ]; then

                ZIPPED=${s}.d${d}.zip
                WORKING=${s}.d${d}.zip.working

                # and there is no zip file yet
                if [ ! -f "$ZIPPED" ]; then
                    # ... remove old files,
                    rm -rf "$WORKING"

                    # ... zip,
                    zip "$WORKING" ${s}.d${d}.* || (
                        printf "$0 [$LINENO]: Can't zip ts files, aborting\n"
                        exit 1
                    )

                    # ... and copy archive
                    mv "$WORKING" "$ZIPPED"
                fi
            fi
            
        done
    done
}

# Zip WRF log files
# Assumes:
#    CWD is the wrf 'run' directory    
function zip_log {
    ZIPPED="logs.zip"
    WORKING="logs.zip.working"
    FILES="rsl.out.0000 tslist namelist.input"

    # check if the log files exist
    for f in $FILES; do
        if [ ! -f $f ]; then
            printf "$0 [$LINENO]: Can't find %s, aborting\n", "$f"
            exit 1
        fi
    done

    # If there is no zip file yet
    if [ ! -f "$ZIPPED" ]; then

        # ... remove old files,
        rm -rf "$WORKING"

        # ... zip,
        zip "$WORKING" $FILES || (
            printf "$0 [$LINENO]: Can't zip log files %s, aborting\n", "$FILES"
            exit 1
        )

        # ... and copy archive
        mv "$WORKING" "$ZIPPED"
    fi
}

# Zip WRF netCDF output
# Assumes:
#    CWD is the wrf 'run' directory    
# Required env:
#    NDOMS, DATESTART
function zip_netcdf {
    if [[ -z "$NDOMS" || -z "$DATESTART" ]]; then
        printf "$0 [$LINENO]: NDOMS or DATESTART not set, aborting\n"
        exit 1;
    fi

    for i in `seq -f '%02.0f' 1 $NDOMS`; do
        NCDF3="wrfout_d${i}_${DATESTART}_00:00:00"
        TMPFL="wrfout_d${i}_${DATESTART}.working"
        NCDF4="wrfout_d${i}_${DATESTART}_00:00:00.nc" 

        # if the output file exists
        if [ -f "${NCDF3}" ]; then

            # .. but the netCDF4 version does not exist
            if [ ! -f "${NCDF4}" ]; then

                # convert the file from netCDF3 to netCDF4
                $NC3TONC4 -o "${NCDF3}" "${TMPFL}" && mv -f "${TMPFL}" "${NCDF4}" 
            fi

            # test if all went ok
            $NCDUMP -h "${NCDF3}" > "${NCDF3}.dump" 
            $NCDUMP -h "${NCDF4}" > "${NCDF4}.dump" 
            cmp "${NCDF3}.dump" "${NCDF4}.dump" || (
                printf "$0 [$LINENO]: Error in converting ${NCDF3}"
                exit 1 
            )

            # clean up
            rm -f "${NCDF3}" "${NCDF3}.dump" "${NCDF4}.dump" "${TMPFL}"
        fi
    done
}

# extract options and their arguments into variables.
case "$1" in
    prepare)
        case "$2" in
            "date")
                echo "Setting date"
                prepare_date $3
            ;;
        esac ;;
    run) 
        # TODO
        echo TODO
    ;;
    zip)
        case "$2" in
            "all")
                echo "Zipping netcdf"
                zip_netcdf
                echo "Zipping TS"
                zip_ts
                echo "Zipping log"
                zip_log
            ;;
            "netcdf") 
                echo "Zipping netcdf"
                zip_netcdf
            ;;
            "ts") 
                echo "Zipping TS"
                zip_ts
            ;;
            "log") 
                echo "Zipping log"
                zip_log
            ;;
            *) echo "Internal error!" ; exit 1 ;;
        esac ;;
    clean)
        case "$2" in
            "all")
                echo "Cleaning input"
                clean_input
                echo "Cleaning output"
                cleen_output
            ;;
            "input")
                echo "Cleaning input"
                clean_input
            ;;
            "output")
                echo "Cleaning output"
                cleen_output
            ;;
            esac
    ;;
    help)
        help
    ;;
    status)
        status
        exit
    ;;
    archive)
        case "$2" in
        *) echo TOOD
        esac ;;
    *) echo "Internal error!" ; exit 1 ;;
esac 

log "$1 $2 $3"

## dictionaries in bash
#declare -A keyValues
#keyValues["key1"]=val1
#keyValues+=( ["key2"]=val2 ["key3"]=val3 )
#for key in ${!keyValues[@]}
#do
#    echo ${key} ${keyValues[${key}]}
#done
#
#echo "Value of key1: "${keyValues[key1]}


# float TSLB(Time, soil_layers_stag, south_north, west_east) ; TSLB:description = "SOIL TEMPERATURE" ;
# float SMOIS(Time, soil_layers_stag, south_north, west_east) ; SMOIS:description = "SOIL MOISTURE" ;
# float SH2O(Time, soil_layers_stag, south_north, west_east) ; SH2O:description = "SOIL LIQUID WATER" ;
# float SMCREL(Time, soil_layers_stag, south_north, west_east) ; SMCREL:description = "RELATIVE SOIL MOISTURE" ;
# float CANWAT(Time, south_north, west_east) ; CANWAT:description = "CANOPY WATER" ;

# int FNDSOILW(Time) ; FNDSOILW:description = "SOILW_LOGICAL" ;

# float TSK(Time, south_north, west_east) ; TSK:description = "SURFACE SKIN TEMPERATURE" ;
# float TMN(Time, south_north, west_east) ; TMN:description = "SOIL TEMPERATURE AT LOWER BOUNDARY" ;
# float TBL_URB_IN(Time, soil_layers_stag, south_north, west_east) ; TBL_URB_IN:description = "WALL LAYER TEMPERATURE" ;
# float TGL_URB_IN(Time, soil_layers_stag, south_north, west_east) ; TGL_URB_IN:description = "ROAD LAYER TEMPERATURE" ;
# float TRL_URB_IN(Time, soil_layers_stag, south_north, west_east) ; TRL_URB_IN:description = "ROOF LAYER TEMPERATURE" ;

#	float TSLB(Time, soil_layers_stag, south_north, west_east) ;
#		TSLB:FieldType = 104 ;
#		TSLB:MemoryOrder = "XYZ" ;
#		TSLB:description = "SOIL TEMPERATURE" ;
#		TSLB:units = "K" ;
#		TSLB:stagger = "Z" ;
#		TSLB:coordinates = "XLONG XLAT" ;
#	float SMOIS(Time, soil_layers_stag, south_north, west_east) ;
#		SMOIS:FieldType = 104 ;
#		SMOIS:MemoryOrder = "XYZ" ;
#		SMOIS:description = "SOIL MOISTURE" ;
#		SMOIS:units = "m3 m-3" ;
#		SMOIS:stagger = "Z" ;
#		SMOIS:coordinates = "XLONG XLAT" ;
#	float SH2O(Time, soil_layers_stag, south_north, west_east) ;
#		SH2O:FieldType = 104 ;
#		SH2O:MemoryOrder = "XYZ" ;
#		SH2O:description = "SOIL LIQUID WATER" ;
#		SH2O:units = "m3 m-3" ;
#		SH2O:stagger = "Z" ;
#		SH2O:coordinates = "XLONG XLAT" ;
#	float SMCREL(Time, soil_layers_stag, south_north, west_east) ;
#		SMCREL:FieldType = 104 ;
#		SMCREL:MemoryOrder = "XYZ" ;
#		SMCREL:description = "RELATIVE SOIL MOISTURE" ;
#		SMCREL:units = "" ;
#		SMCREL:stagger = "Z" ;
#		SMCREL:coordinates = "XLONG XLAT" ;
#	float CANWAT(Time, south_north, west_east) ;
#		CANWAT:FieldType = 104 ;
#		CANWAT:MemoryOrder = "XY " ;
#		CANWAT:description = "CANOPY WATER" ;
#		CANWAT:units = "kg m-2" ;
#		CANWAT:stagger = "" ;
#		CANWAT:coordinates = "XLONG XLAT" ;
#
#	int FNDSOILW(Time) ;
#		FNDSOILW:FieldType = 106 ;
#		FNDSOILW:MemoryOrder = "0  " ;
#		FNDSOILW:description = "SOILW_LOGICAL" ;
#		FNDSOILW:units = "-" ;
#		FNDSOILW:stagger = "" ;
#
#	float TSK(Time, south_north, west_east) ;
#		TSK:FieldType = 104 ;
#		TSK:MemoryOrder = "XY " ;
#		TSK:description = "SURFACE SKIN TEMPERATURE" ;
#		TSK:units = "K" ;
#		TSK:stagger = "" ;
#		TSK:coordinates = "XLONG XLAT" ;
#	float TMN(Time, south_north, west_east) ;
#		TMN:FieldType = 104 ;
#		TMN:MemoryOrder = "XY " ;
#		TMN:description = "SOIL TEMPERATURE AT LOWER BOUNDARY" ;
#		TMN:units = "K" ;
#		TMN:stagger = "" ;
#		TMN:coordinates = "XLONG XLAT" ;
#	float TBL_URB_IN(Time, soil_layers_stag, south_north, west_east) ;
#		TBL_URB_IN:FieldType = 104 ;
#		TBL_URB_IN:MemoryOrder = "XYZ" ;
#		TBL_URB_IN:description = "WALL LAYER TEMPERATURE" ;
#		TBL_URB_IN:units = "K" ;
#		TBL_URB_IN:stagger = "Z" ;
#		TBL_URB_IN:coordinates = "XLONG XLAT" ;
#	float TGL_URB_IN(Time, soil_layers_stag, south_north, west_east) ;
#		TGL_URB_IN:FieldType = 104 ;
#		TGL_URB_IN:MemoryOrder = "XYZ" ;
#		TGL_URB_IN:description = "ROAD LAYER TEMPERATURE" ;
#		TGL_URB_IN:units = "K" ;
#		TGL_URB_IN:stagger = "Z" ;
#		TGL_URB_IN:coordinates = "XLONG XLAT" ;
#	float TRL_URB_IN(Time, soil_layers_stag, south_north, west_east) ;
#		TRL_URB_IN:FieldType = 104 ;
#		TRL_URB_IN:MemoryOrder = "XYZ" ;
#		TRL_URB_IN:description = "ROOF LAYER TEMPERATURE" ;
#		TRL_URB_IN:units = "K" ;
#		TRL_URB_IN:stagger = "Z" ;
#		TRL_URB_IN:coordinates = "XLONG XLAT" ;
#
