#!/bin/bash

# abort on any error (ie. non-zero exit status)
set -e

###########################3
# Forecast config:

CYCLESTEP="1 day"  # time between two forecast runs
CYCLELEN="2 days"  # length of a forecast run

# Index in netCDF file to use for copy_urban and copy_surface
CYCLEINDEX[01]=0
CYCLEINDEX[02]=0
CYCLEINDEX[03]=0
CYCLEINDEX[04]=0
CYCLEINDEX[05]=0

# location of configuration file
CONFIG=/home/jiska/forecast.config

# working directories
DATDIR=/home/jiska/WRF/tars/dat
WPSDIR=/home/jiska/WRF/tars/WPS
RUNDIR=/home/jiska/WRF/tars/WRFV3/run
ARCDIR=/home/jiska/archive

# location of external tools
NCDUMP=ncdump
NC3TONC4=nc3tonc4
TOOLS=/home/jiska/WRF/tars/WRFV3/tools
NAMELIST=$TOOLS/namelist.py
COPYURBAN=$TOOLS/copy_urb_init.sh


dateregex="([0-9][0-9][0-9][0-9]).([0-9][0-9]).([0-9][0-9])"


##################################
# FIXME read from config file
##################################


MANUAL="
Control WRF forecast runs.

   $0 command option

prepare:
  all
  date  <date> Set the datetime for the run, also accepts special date 'next', 
  boundaries   Download boundaries from NCEP and run ungrib, metgrid
  urban <date> Copy urban fields from the specified run, also accepts special date 'previous'
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

function forecastinit {
    NDOMS=`$NAMELIST --get domains:max_dom namelist.input`

    START_Y=`$NAMELIST --get time_control:start_year:0  namelist.input`
    START_M=`$NAMELIST --get time_control:start_month:0 namelist.input`
    START_D=`$NAMELIST --get time_control:start_day:0   namelist.input`
    DATESTART=`printf '%4i-%02i-%02i' ${START_Y} ${START_M} ${START_D}`

    if [ -f $RUNDIR/tslist ]; then
        STATIONS=`cat $RUNDIR/tslist | awk '{print $2}'`
    else
        printf "$0 [$LINENO]: Can't open tslist file!\n"
        STATIONS=""
    fi
}


function repeat {
    str=$1
    count=$2
   
    printf "%${count}s" | sed "s/ /${str},/g"
}


function log {
    message=$1

    stamp=`date -u +"[%s] %F %T"`
    printf '%s : %s\n' "$stamp" "$message" >> "$CONFIG"
}

# Argument:
#  station
#  domain
function ls_tslist {

    # suffixes for the ts files
    SX="(PH|QV|TH|TS|UU|VV)"

    RX=`printf '%s.d%02.0f.%s' $1 $2 $SX`

    ls -1 $RUNDIR | grep -P $RX || echo ""
}


function help {
    printf '%s\n' "$MANUAL"
}

function status {
    echo "Date start: $DATESTART"
    echo "NDOMS:      $NDOMS"

    # Check archive status
    archivedir $DATESTART ARCHIVE

    echo '------------------------------'
    echo 'Archive status'
    echo '------------------------------'
    echo -n "LOGS        "
    if [[ -f "$ARCHIVE/logs_${DATESTART}.zip" ]]; then
        echo done
    else
        if [[ -f "$RUNDIR/logs_${DATESTART}.zip" ]]; then
            echo "run archive log"
        else
            echo .
        fi
    fi

    # check if the wrfout*nc files are archived...
    echo -n "NETCDF:     " 
    if [[ `ls $ARCHIVE/wrfout*nc 2>/dev/null | wc -l ` == $NDOMS ]]; then
        echo done
    else
        # ... or if we can do archive netcdf
        if [[ `ls $RUNDIR/wrfout*nc` ]]; then
            echo "run archive netcdf"
        else
            echo .
        fi
    fi

    # As the number of time series files depends on if the station are within which nested sub domains
    # we cannot easily find out how many there should be
    echo -n "TS:         "
    if [ -f $RUNDIR/tslist ]; then
           
        # are there already some files archived?
        if [[ `ls $ARCHIVE/*.d??.zip 2>/dev/null` ]]; then

            # ...are there still files to archive?
            if [[ `ls $RUNDIR/*.d??.zip 2>/dev/null` ]]; then
                echo "run archive ts"
            else
                echo done
            fi
        else       
            # ...are there files to archive?
            if [[ `ls $RUNDIR/*.d??.zip 2>/dev/null` ]]; then
                echo "run archive ts"
            else
                echo .
            fi
        fi
    else
        echo "no tslist"
    fi

    echo '------------------------------'
    echo 'Zip status'
    echo '------------------------------'

    echo -n "LOGS:       "
    if [[ -f $RUNDIR/rsl.out.0000 || -f $RUNDIR/rsl.out ]]; then
         echo "run zip log"
    elif [[ -f "$RUNDIR/logs_${DATESTART}.zip" ]]; then
         echo done
    else
         echo .
    fi

    echo -n "NETCDF:     " 
    if [[ "$NDOMS" == `ls $RUNDIR/wrfout*nc 2>/dev/null | wc -l` ]]; then
        echo done
    else
        # check if the wrfout_d* files exist
        if [[ "$NDOMS" == `ls $RUNDIR/wrfout_d* 2> /dev/null | grep -v nc` ]]; then
            echo "run zip netcdf"
        else
            echo .
        fi
    fi

    echo -n "TS:         "
    FILES=""
    for s in $STATIONS; do
        for d in `seq 1 $NDOMS`; do
            FILES+=`ls_tslist $s $d`
        done
    done
    if [ ! -z "$FILES" ]; then
        echo "run zip ts"
    else
        echo .
    fi

    echo '------------------------------'
    echo 'Run status'
    echo '------------------------------'

    # output filename depends on parallel / serial run
    if [ -f rsl.out.0000 ]; then
        OUT="rsl.out.0000"
    elif [ -f rsl.out ]; then
        OUT="rsl.out"
    fi

    NEEDLE="SUCCES COMPLETE BOUNDARIES $DATESTART"
    BDY="not done"
    grep "${NEEDLE}" prepare_boundaries.log 1>/dev/null 2>&1 && BDY="done"

    WRF="not done"
    NEEDLE="SUCCESS COMPLETE WRF"
    grep "SUCCESS COMPLETE WRF" "$OUT" 1>/dev/null 2>&1 && WRF="done"


    echo -n "REAL:       "
    if [ "$BDY" == "done" ]; then
        echo "done"
    else
        echo "."
    fi

    echo -n "WRF:        "
    if [ "$WRF" == "done" ]; then
        echo "done"
    else 
        if [ "$BDY" == "done" ]; then
            echo "run wrf"
        else
            echo "."
        fi
    fi

    echo '------------------------------'
    echo 'Boundary status'
    echo '------------------------------'

echo
    
}

function splitdate {
    if [[ ! "$1" =~ $dateregex ]]; then
        printf "$0 [$LINENO]: No valid date given. Aborting.\n"
        exit 1
    fi
    if [[ $# != 4 ]]; then
        printf "$0 [$LINENO]: splitdate needs exactly 4 arguments. Aborting.\n"
        exit 1
    fi

    eval "$2='${BASH_REMATCH[1]}'"
    eval "$3='${BASH_REMATCH[2]}'"
    eval "$4='${BASH_REMATCH[3]}'"
}


# Create the archive dir for the given date; date as YYYY-MM-DD
# Dir is returned in the second argument $2, so use like this
# archivedir 2014-01-01 ARCDIR
function archivedir {

    # Construct archive directory name
    splitdate $1 YEAR MONTH DAY
    DEST=$ARCDIR/$YEAR/$MONTH/$DAY

    # Create it
    mkdir -p $DEST

    # return it in the second argument
    eval "$2='$DEST'"
}


# Move TS files to the archive
#
# Assumes:
#    tslist is in $RUNDIR
# Required env:
#    ARCDIR, RUNDIR, DATESTART, NDOMS
function archive_ts {
    if [[ -z "$NDOMS" || -z "$DATESTART" || ! -d "$ARCDIR" || ! -d "$RUNDIR" ]]; then
        printf "$0 [$LINENO]: NDOMS, DATESTART, ARCDIR, or RUNDIR not set. Aborting\n"
        exit -1
    fi;

    # Get the archive directory for this run
    archivedir $DATESTART DEST

    for d in `seq -f '%02.0f' 1 $NDOMS`; do
        for s in $STATIONS; do

            ZIPPED=${s}.d${d}.zip
            if [ -f $RUNDIR/$ZIPPED ]; then
                cp $RUNDIR/$ZIPPED $DEST
                rm -f $RUNDIR/$ZIPPED
            fi
        done
    done
}

# Move netCDF files to the archive
#
# Required env:
#    ARCDIR, RUNDIR, DATESTART, NDOMS
function archive_netcdf {
    if [[ -z "$NDOMS" || -z "$DATESTART" || ! -d "$ARCDIR" || ! -d "$RUNDIR" ]]; then
        printf "$0 [$LINENO]: NDOMS, DATESTART, ARCDIR, or RUNDIR not set. Aborting\n"
        exit -1
    fi;

    # Get the archive directory for this run
    archivedir $DATESTART DEST

    for d in `seq -f '%02.0f' 1 $NDOMS`; do
        NCDF4="wrfout_d${d}_${DATESTART}_00:00:00.nc" 
        cp $RUNDIR/$NCDF4 $DEST
        rm -f $RUNDIR/$NCDF4
    done
}

# Move log files to the archive
#
# Required env:
#    ARCDIR, RUNDIR, DATESTART, NDOMS
function archive_log {
    if [[ -z "$NDOMS" || -z "$DATESTART" || ! -d "$ARCDIR" || ! -d "$RUNDIR" ]]; then
        printf "$0 [$LINENO]: NDOMS, DATESTART, ARCDIR, or RUNDIR not set. Aborting\n"
        exit -1
    fi;

    # Get the archive directory for this run
    archivedir $DATESTART $DEST

    LOGS="logs_${DATESTART}.zip"

    if [ -f $RUNDIR/$LOGS ]; then
        cp $RUNDIR/$LOGS $DEST
        rm -f $RUNDIR/$LOGS
    else
        printf "$0 [$LINENO]: Cannot archive logs. Aborting\n"
        exit -1
    fi
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
    rm -f wrfbdy*    
    rm -f met_em.d??.????-??-??_??:??:??.nc

    # logs
    # ----

    rm -f prepare_boundaries.log
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

    rm -f rsl.out
    rm -f rsl.out.*
    rm -f rs.error.*
    rm -f $LOGS

    # wrfout
    # ------------------------

    rm -f wrfout*

    # time series
    # -----------

    if [ -f tslist ]; then

        for d in `seq -f '%02.0f' 1 $NDOMS`; do
            for s in $STATIONS; do
                rm -f ${s}.d${d}.*
            done
        done
    else
        printf "$0 [$LINENO]: Can't open tslist file, time series not removed\n"
    fi
}

# Modifies the namelist.input to for the required date
# Arguments:
#      $1   The date as YYYY-MM-DD or 'next'
#
# When 'next' is given, CYCLESTEP is added to the start date from the namelist.
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

    if [ ! "$1" ]; then
        printf "$0 [$LINENO]: No valid date or 'next' given, aborting\n"
        exit 1;
    fi

    if [ "next" == $1 ]; then
        DATESTART=`date --date "$DATESTART $CYCLESTEP" +%F`
    else
        DATESTART=$1
    fi

    # Set starting date

    splitdate $DATESTART YEAR MONTH DAY
    $NAMELIST --set time_control:start_year  `repeat $YEAR  $NDOMS` namelist.input
    $NAMELIST --set time_control:start_month `repeat $MONTH $NDOMS` namelist.input
    $NAMELIST --set time_control:start_day   `repeat $DAY   $NDOMS` namelist.input

    # Set ending date
    DATEEND=`date --date "$DATESTART $CYCLELEN" +%F`

    splitdate $DATEEND YEAR MONTH DAY
    $NAMELIST --set time_control:end_year  `repeat $YEAR  $NDOMS` namelist.input
    $NAMELIST --set time_control:end_month `repeat $MONTH $NDOMS` namelist.input
    $NAMELIST --set time_control:end_day   `repeat $DAY   $NDOMS` namelist.input
}


# Run the standard WRF commands to make input boundaries
# Required env:
#    WPSDIR, DATDIR, RUNDIR
function prepare_boundaries {
    OLDCWD=`pwd`

    if [[ ! -d "$WPSDIR"  || ! -d "$DATDIR" || ! -d "$RUNDIR" ]]; then
        printf "$0 [$LINENO]: One of WPSDIR, DATDIR, or RUNDIR not set. Aborting\n"
        exit -1
    fi;

    # clean start
    rm -f GRIBFILE.???
    rm -f FILE:????-??-??_??

    cd $WPSDIR
    rm -f prepare_boundaries.log
    ./link_grib.csh "$DATDIR/*" 2>&1 >> $RUNDIR/prepare_boundaries.log
    ./ungrib.exe                2>&1 >> $RUNDIR/prepare_boundaries.log
    ./metgrid.exe               2>&1 >> $RUNDIR/prepare_boundaries.log

    rm -f GRIBFILE.???
    rm -f FILE:????-??-??_??
    mv met_em.d??.????-??-??_??:??:??.nc $RUNDIR
    cd $OLDCWD
}

# Copy urban temperature fields TRL TBL TGL for cycling.
# When no date and index is given, try copying form DATESTART - CYCLESTEP
# Note that the index in the netcdf file (CYCLEINDEX) is hardcoded at the start
# of this file
#
# Arguments:
#    CYCLEDATE (optional)
# Required env:
#    RUNDIR, NDOMS
function prepare_urban {
    OLDCWD=`pwd`

    if [[ -z "$NDOMS" || ! -d "$RUNDIR" ]]; then
        printf "$0 [$LINENO]: NDOMS or RUNDIR not set. Aborting\n"
        exit -1
    fi
    if [ "previous" == "$1" ]; then
        CYCLEDATE=`date --date "$DATESTART $CYCLESTEP ago" +%F`
    else
        CYCLEDATE=$1
    fi
    archivedir $CYCLEDATE ARCDIR

    for d in `seq -f '%02.0f' 1 $NDOMS`; do
       $COPYURBAN ${ARCDIR}/wrfout_d${d}_${CYCLEDATE}_00:00:00.nc ${CYCLEINDEX[$d]} ${RUNDIR}/wrfinput_d${d}
    done
}

# Run real.exe
#
# Requires:
#    RUNDIR
#
function run_real {
    if [ ! -d "$RUNDIR" ]; then
        printf "$0 [$LINENO]: RUNDIR not set. Aborting\n"
        exit -1
    fi
    cd $RUNDIR 
    ./real.exe 2>&1 >> prepare_boundaries.log
}

# Run WRF
#
# Requires:
#    RUNDIR
#
function run_wrf {
    if [ ! -d "$RUNDIR" ]; then
        printf "$0 [$LINENO]: RUNDIR not set. Aborting\n"
        exit -1
    fi
    cd $RUNDIR 
    ./wrf.exe 2>&1
}


# Zip WRF timeseries files
#
#    station.d??.* -> station.d??.zip
#
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

    for d in $NDOMS; do
        for s in $STATIONS; do

            # if there are ts files for this station and domain to process,
            FILES=`ls_tslist $s $d`
            if [ ! -z "$FILES" ]; then

                ZIPPED=${s}.d${d}.zip
                WORKING=working.${s}.d${d}.zip

                # and there is no zip file yet
                if [ ! -f "$ZIPPED" ]; then

                    # ... remove old archive to be sure,
                    rm -f "$WORKING"

                    # ... zip,
                    zip "$WORKING" $FILES

                    # ... clean up
                    rm -f $FILES

                    # ... and copy archive
                    mv "$WORKING" "$ZIPPED"
                fi
            fi
            
        done
    done
}

# Zip WRF log files
#
#     $LOGS: tslist namelist.input prepare_boundaries.log rsl.out|rsl.out.0000  
#
# Assumes:
#    CWD is the wrf 'run' directory    
function zip_log {
    WORKING="${LOGS}.working"
    FILES="tslist namelist.input prepare_boundaries.log"

    # output filename depends on parallel / serial run
    if [ -f rsl.out.0000 ]; then
        FILES="rsl.out.0000 $FILES"
    elif [ -f rsl.out ]; then
        FILES="rsl.out $FILES"
    else
        printf "$0 [$LINENO]: Can't find rsl.out or rsl.out.0000, aborting\n"
        exit 1
    fi

    # check if the log files exist
    for f in $FILES; do
        if [ ! -f $f ]; then
            printf "$0 [$LINENO]: Can't find %s, aborting\n", "$f"
            exit 1
        fi
    done

    # If there is no zip file yet
    if [ ! -f "$LOGS" ]; then

        # ... remove old files,
        rm -rf "$WORKING"

        # ... zip,
        zip "$WORKING" $FILES || (
            printf "$0 [$LINENO]: Can't zip log files %s, aborting\n", "$FILES"
            exit 1
        )

        # ... and copy archive
        mv "$WORKING" "$LOGS"
    fi
}

# Zip WRF netCDF output by converting to netCDF4:
#
#      wrfout* -> wrfout*.nc
#
# Assumes:
#    CWD is the wrf 'run' directory    
# Required env:
#    NDOMS, DATESTART
function zip_netcdf {
    if [[ -z "$NDOMS" || -z "$DATESTART" ]]; then
        printf "$0 [$LINENO]: NDOMS or DATESTART not set, aborting\n"
        exit 1;
    fi

    for d in `seq -f '%02.0f' 1 $NDOMS`; do
        NCDF3="wrfout_d${d}_${DATESTART}_00:00:00"
        TMPFL="wrfout_d${d}_${DATESTART}.working"
        NCDF4="wrfout_d${d}_${DATESTART}_00:00:00.nc" 

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



#####################################################
#            M A I N
#####################################################

ls_tslist  veenk 10

# Immediately print help and exit
if [[ "$1" == "help" || "$1" == "-h" || "$1" == "-?" ]]; then
    help
    exit 0
fi


forecastinit

# extract options and their arguments into variables.
case "$1" in
    prepare)
        case "$2" in
            "date")       prepare_date $3 ;;
            "boundaries") prepare_boundaries ;;
            "urban")      prepare_urban $3 ;;
        esac
    ;;
    run) 
        case "$2" in
            "real") run_real ;;
            "wrf")  run_wrf  ;;
        esac
    ;;
    zip)
        case "$2" in
            "all")     zip_netcdf ; zip_ts : zip_log ;;
            "netcdf")  zip_netcdf ;;
            "ts")      zip_ts ;;
            "log")     zip_log ;;
            *)         echo "Internal error!" ; exit 1 ;;
        esac
    ;;
    clean)
        case "$2" in
            "all")     clean_input clean_output ;;
            "input")   clean_input ;;
            "output")  clean_output ;;
            *)         echo "Internal error!" ; exit 1 ;;
            esac
    ;;
    status) status ; exit
    ;;
    archive)
        case "$2" in
            "all")      archive_ts ; archive_netcdf ; archive_log ;;
            "ts")       archive_ts ;;
            "netcdf")   archive_netcdf ;;
            "log")      archive_log ;;
            *)         echo "Internal error!" ; exit 1 ;;
        esac ;;
    *)
        help
        printf "$0 [$LINENO]: Please indicate an action.\n"
        exit 1
    ;;
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
