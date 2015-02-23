#!/bin/bash

# abort on any error (ie. non-zero exit status)
set -e

# Uncomment following line for debug mode
#set -xv

###########################3
# Forecast config:

CYCLESTEP=24        # time between two forecast runs in hours
CYCLELEN=24         # length of a forecast run in hours
BOUNDARYINTERVAL=6  # time between boundaries, in hours

# Index in netCDF file to use for copy_cycle
CYCLEINDEX[01]=0
CYCLEINDEX[02]=0
CYCLEINDEX[03]=0
CYCLEINDEX[04]=0
CYCLEINDEX[05]=0

# location of configuration file
CONFIG=/home/jiska/forecast.config

# working directories
DATDIR=/data/hdd/ECMWF_OPAN
WPSDIR=/home/jiska/WRF/tars/WPS
RUNDIR=/home/jiska/WRF/tars/WRFV3/run
ARCDIR=/home/jiska/archive

# location of external tools
NCDUMP=ncdump
NC3TONC4=nc3tonc4
TOOLS=/home/jiska/WRF/tars/WRFV3/tools/forecast
NAMELIST=$TOOLS/namelist.py
COPYURBAN=$TOOLS/copy_urb_init.sh
COPYCYCLE=$TOOLS/copy_cycle.sh


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
  cycle <date> Copy cycle fields from the specified run, also accepts special date 'previous'
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

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#                    U  T  I  L  I  T  Y     F  U  N  C  T  I  O  N  S 
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

######################################################################
# Repeat a string str count times, using. Output is comma separated: `repeat 'x' 3` = "x,x,x,"
# Arguments:
#    str       string to repeat
#    count     number of times to repeat it 
######################################################################
function repeat {
    str=$1
    count=$2
   
    printf "%${count}s" | sed "s/ /${str},/g"
}

######################################################################
# Write a log message to the CONFIG file
# Arguments:
#    message   Mesage to write
######################################################################
function log {
    message=$1

    stamp=`date -u +"[%s] %F %T"`
    printf '%s : %s\n' "$stamp" "$message" >> "$CONFIG"
}

######################################################################
# Split a date given as YYYY.MM.DD into year, month, and day
# Arguments:
#    date     Date to parse
#    year     (return) name of output variable for the year
#    month    (return) name of output variable for the month
#    day      (return) name of output variable for the day
######################################################################
function splitdate {
    dateregex="([0-9][0-9][0-9][0-9]).([0-9][0-9]).([0-9][0-9])"

    if [[ ! "$1" =~ $dateregex ]]; then
        printf "$0 [$LINENO]: No valid date given. Aborting.\n"
        exit 1
    fi
    if [[ $# != 4 ]]; then
        printf "$0 [$LINENO]: splitdate needs exactly 4 arguments, not $#. Aborting.\n"
        exit 1
    fi

    eval "$2='${BASH_REMATCH[1]}'"
    eval "$3='${BASH_REMATCH[2]}'"
    eval "$4='${BASH_REMATCH[3]}'"
}

######################################################################
# List all available tslist files for a given station and domain
# Required env:
#    RUNDIR
# Argument:
#    station
#    domain
######################################################################
function ls_tslist {
    station=$1
    domain=$2

    # suffixes for the ts files
    SX="(PH|QV|TH|TS|UU|VV)"

    RX=`printf '%s.d%02.0f.%s' $1 $2 $SX`

    ls -1 $RUNDIR | grep -P $RX || echo ""
}

######################################################################
# Return a list of boundary filenames for the ECMWF operational analysis
# Required env:
#    BOUNDARYINTERVAL, CYCLELEN
# Arguments:
#    DATE      starting date in a format accepted by date
#    FILES     (return) variable that will be set to required files
######################################################################
function boundary_list_ecmwf_opan {
    if [[ -z "$BOUNDARYINTERVAL" || -z "$CYCLELEN" ]]; then
        printf "$0 [$LINENO]: BOUNDARYINTERVAL or CYCLELEN not set. Aborting\n"
        exit -1
    fi;

    start=$1
    out=$2

    LIST=""
    for i in `seq 0 $BOUNDARYINTERVAL $CYCLELEN`; do

        YMD=`date --date "$start $i hours" +%Y%m%d`
        HM=`date --date "$start $i hours" +%_k | tr -d ' '`
        if [[ "$HM" != "0" ]]; then
            HM=${HM}00
        fi
        LIST="$LIST $DATDIR/AN${YMD}${HM}sig"
    done

    eval "$out='${LIST}'"
}

######################################################################
# Initialize the forecast script
# Parses the namelist file and sets NDOMS, DATESTART, and STATIONS
# Required env:
#    NAMELIST, RUNDIR 
######################################################################
function forecastinit {
    if [[ -z "$NAMELIST" || -z "$RUNDIR" ]]; then
        printf "$0 [$LINENO]: NAMELIST or RUNDIR not set. Aborting\n"
        exit -1
    fi;

    NDOMS=`$NAMELIST --get domains:max_dom "$RUNDIR/namelist.input"`

    START_Y=`$NAMELIST --get time_control:start_year:0  "$RUNDIR/namelist.input"`
    START_M=`$NAMELIST --get time_control:start_month:0 "$RUNDIR/namelist.input"`
    START_D=`$NAMELIST --get time_control:start_day:0   "$RUNDIR/namelist.input"`
    DATESTART=`printf '%4i-%02i-%02i' ${START_Y} ${START_M} ${START_D}`

    if [ -f $RUNDIR/tslist ]; then
        STATIONS=`cat $RUNDIR/tslist | awk '{print $2}'`
    else
        printf "$0 [$LINENO]: Can't open tslist file!\n"
        STATIONS=""
    fi
}

######################################################################
# Check GRIB edition number of ECMWF files, 
# and convert to GRIB1 if necessary
# This will overwrite the original file
# Arguments:
#    FILES  List of files to process
######################################################################
function check_grib_version {
    if [ -z "$1" ]; then
        printf "$0 [$LINENO]: No files given, aborting\n"
        exit -1
    fi;

    F=$1
    while [ ! -z "$F" ]; do
        if [ ! -f "$F" ]; then
            printf "$0 [$LINENO]: File not found: $F, aborting\n"
            exit -1
        fi
        if [[ `grib_ls "$F" | grep  -P '^2 *ecmf'` ]]; then
            # Convert to GRIB1
            grib_set -f -s editionNumber=1 "${F}" "${F}.grib1" > /dev/null 2>&1
            mv -f "${F}.grib1" "${F}"
        fi
        shift
        F=$1 
    done
}


######################################################################
# Create the archive dir for the given date; date as YYYY-MM-DD
# Dir is returned in the second argument $2, so use like this
# archivedir 2014-01-01 ARCDIR
######################################################################
function archivedir {

    # Construct archive directory name
    splitdate $1 YEAR MONTH DAY
    DEST=$ARCDIR/$YEAR/$MONTH/$DAY

    # Create it
    mkdir -p $DEST

    # return it in the second argument
    eval "$2='$DEST'"
}




# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#                    F  o  r  e c  a  s  t     c  o  m  m  a  n  d  s 
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



######################################################################
# Print help and exit
######################################################################
function help {
    printf '%s\n' "$MANUAL"
}


######################################################################
# Determine forecast status
######################################################################
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
    if [[ `ls $ARCHIVE/wrfout*nc 2> /dev/null | wc -l` == $NDOMS ]]; then
        echo done
    else
        # ... or if we can do archive netcdf
        if [[ `ls $RUNDIR/wrfout*nc 2> /dev/null` ]]; then
            echo "run archive netcdf"
        else
            echo .
        fi
    fi

    # As the number of time series files depends on if the stations are within which nested sub domains
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
        if [[ `ls $RUNDIR/wrfout_d* 2> /dev/null | grep -v nc` ]]; then
            echo "run zip netcdf"
        else
            echo .
        fi
    fi

    echo -n "TS:         "
    foundts="no"
    foundtszip="no"
    for s in $STATIONS; do
        for d in `seq 1 $NDOMS`; do

            # Check for un-zipped files
            FILES=`ls_tslist $s $d`
            for F in $FILES; do
                if [ -f $F ]; then
                    foundts="yes"
                fi
            done

            # Check for zip archive
            ZIPPED=${s}.d${d}.zip
            if [ -f $ZIPPED ]; then
                foundtszip="yes"
            fi

        done
    done
    if [ $foundts == "yes" ]; then
        echo "run zip ts"
    else
        if [ $foundtszip == "yes" ]; then
            echo "done"
        else
            echo .
        fi
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

    NEEDLE="SUCCESS COMPLETE REAL_EM INIT"
    BDY="not done"
    grep "${NEEDLE}" prepare_boundaries.log 1>/dev/null 2>&1 && BDY="done"

    # Test for the wrfbdy_d?? and wrfinput_d??.
    # They should always be present for d01, but it depends on the nesting feedback setup if they are required for the other domains.
    if [[ ! -f wrfbdy_d01  || ! -f wrfinput_d01 ]]; then
        BDY="not done"
    fi

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


######################################################################
# Move TS files to the archive
#
# Assumes:
#    tslist is in $RUNDIR
# Required env:
#    ARCDIR, RUNDIR, DATESTART, NDOMS
######################################################################
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

######################################################################
# Move netCDF files to the archive
#
# Required env:
#    ARCDIR, RUNDIR, DATESTART, NDOMS
######################################################################
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

######################################################################
# Move log files to the archive
#
# Required env:
#    ARCDIR, RUNDIR, DATESTART, NDOMS
######################################################################
function archive_log {
    if [[ -z "$NDOMS" || -z "$DATESTART" || ! -d "$ARCDIR" || ! -d "$RUNDIR" ]]; then
        printf "$0 [$LINENO]: NDOMS, DATESTART, ARCDIR, or RUNDIR not set. Aborting\n"
        exit -1
    fi;

    # Get the archive directory for this run
    archivedir $DATESTART DEST

    LOGS="logs_${DATESTART}.zip"

    if [ -f $RUNDIR/$LOGS ]; then
        cp $RUNDIR/$LOGS $DEST
        rm -f $RUNDIR/$LOGS
    else
        printf "$0 [$LINENO]: Cannot archive logs. Aborting\n"
        exit -1
    fi
}

######################################################################
# Remove all input files from a WRF run,
# after a 'clean input', a 'prepare all' must be run
#
# Assumes:
#    CWD is the 'wrf' run directory
######################################################################
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

######################################################################
# Remove all output files from a WRF run,
# after a 'clean output' the WRF forecast can be resubmitted
#
# Assumes:
#    CWD is the 'wrf' run directory
######################################################################
function clean_output {
    FORCE=$1

    # log files
    # -------------------------

    rm -f rsl.out
    rm -f rsl.out.*
    rm -f rs.error.*
    rm -f logs_${DATESTART}.zip

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

######################################################################
# Modifies the namelist.input and namelist.wps to for the requested date
# Arguments:
#      $1   The date as YYYY-MM-DD or 'next'
#
# When 'next' is given, CYCLESTEP is added to the start date from the namelist.
#
# Required env:
#    NDOMS
######################################################################
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
        DATESTART=`date --date "$DATESTART $CYCLESTEP hours" +%F`
    else
        DATESTART=$1
    fi

    # namelist.input 
    # namelist.wps   (Format: 2006-08-16_12:00:00)
    # --------------------------------------------

    # Set starting date

    splitdate $DATESTART YEAR MONTH DAY
    $NAMELIST --set time_control:start_year  `repeat $YEAR  $NDOMS` $RUNDIR/namelist.input
    $NAMELIST --set time_control:start_month `repeat $MONTH $NDOMS` $RUNDIR/namelist.input
    $NAMELIST --set time_control:start_day   `repeat $DAY   $NDOMS` $RUNDIR/namelist.input

    $NAMELIST --set share:start_date `repeat ${YEAR}-${MONTH}-${DAY}_00:00:00 $NDOMS` $WPSDIR/namelist.wps

    # Set ending date
    DATEEND=`date --date "$DATESTART $CYCLELEN hours" +%F`

    splitdate $DATEEND YEAR MONTH DAY
    $NAMELIST --set time_control:end_year  `repeat $YEAR  $NDOMS` $RUNDIR/namelist.input
    $NAMELIST --set time_control:end_month `repeat $MONTH $NDOMS` $RUNDIR/namelist.input
    $NAMELIST --set time_control:end_day   `repeat $DAY   $NDOMS` $RUNDIR/namelist.input

    $NAMELIST --set share:end_date `repeat ${YEAR}-${MONTH}-${DAY}_00:00:00 $NDOMS` $WPSDIR/namelist.wps
}


######################################################################
# Run the standard WRF commands to make input boundaries
# Required env:
#    WPSDIR, DATDIR, RUNDIR
######################################################################
function prepare_boundaries {

    if [[ ! -d "$WPSDIR"  || ! -d "$DATDIR" || ! -d "$RUNDIR" || -z "$DATESTART" ]]; then
        printf "$0 [$LINENO]: One of WPSDIR, DATDIR, RUNDIR, or DATESTART not set. Aborting\n"
        exit -1
    fi;

    boundary_list_ecmwf_opan $DATESTART FILES

    # clean start
    rm -f $WPSDIR/GRIBFILE.???
    rm -f $WPSDIR/FILE:????-??-??_??
    rm -f $WPSDIR/PFILE:????-??-??_??
    rm -f $WPSDIR/PRES:????-??-??_??

    rm -f $RUNDIR/prepare_boundaries.log
    check_grib_version $FILES
    $WPSDIR/link_grib.csh $FILES    2>&1 >> $RUNDIR/prepare_boundaries.log
    $WPSDIR/ungrib.exe              2>&1 >> $RUNDIR/prepare_boundaries.log
    $WPSDIR/util/calc_ecmwf_p.exe   2>&1 >> $RUNDIR/prepare_boundaries.log
    if [[ `grep "^ERROR:" $RUNDIR/prepare_boundaries.log` ]]; then
        printf "$0 [$LINENO]: Error running calc_ecmwf.exe. Aborting\n"
        exit -1
    fi 
    $WPSDIR/metgrid.exe             2>&1 >> $RUNDIR/prepare_boundaries.log

    rm -f $WPSDIR/GRIBFILE.???
    rm -f $WPSDIR/FILE:????-??-??_??
    mv -f $WPSDIR/met_em.d??.????-??-??_??:??:??.nc $RUNDIR
}

######################################################################
# Copy urban temperature fields TRL TBL TGL, and some surface fields for cycling
# When no date and index is given, try copying form DATESTART - CYCLESTEP
# Note that the index in the netcdf file (CYCLEINDEX) is hardcoded at the start
# of this file
#
# Arguments:
#    CYCLEDATE (optional)
# Required env:
#    RUNDIR, NDOMS
######################################################################
function prepare_cycle {

    if [[ -z "$NDOMS" || ! -d "$RUNDIR" ]]; then
        printf "$0 [$LINENO]: NDOMS or RUNDIR not set. Aborting\n"
        exit -1
    fi
    if [ "previous" == "$1" ]; then
        CYCLEDATE=`date --date "$DATESTART $CYCLESTEP hours ago" +%F`
    else
        CYCLEDATE=$1
    fi
    archivedir $CYCLEDATE ARCDIR

    for d in `seq -f '%02.0f' 1 $NDOMS`; do
       $COPYURBAN ${ARCDIR}/wrfout_d${d}_${CYCLEDATE}_00:00:00.nc ${CYCLEINDEX[$d]} ${RUNDIR}/wrfinput_d${d}
       $COPYCYCLE ${ARCDIR}/wrfout_d${d}_${CYCLEDATE}_00:00:00.nc ${CYCLEINDEX[$d]} ${RUNDIR}/wrfinput_d${d}
    done
}

######################################################################
# Set sea-surface temperature
#
# Required env:
#    RUNDIR, NDOMS
######################################################################
function prepare_sst {
}

######################################################################
# Run real.exe
#
# Requires:
#    RUNDIR
######################################################################
function run_real {
    if [ ! -d "$RUNDIR" ]; then
        printf "$0 [$LINENO]: RUNDIR not set. Aborting\n"
        exit -1
    fi
    cd $RUNDIR 
    ./real.exe 2>&1 >> prepare_boundaries.log
}

######################################################################
# Run wrf.exe
#
# Requires:
#    RUNDIR
######################################################################
function run_wrf {
    if [ ! -d "$RUNDIR" ]; then
        printf "$0 [$LINENO]: RUNDIR not set. Aborting\n"
        exit -1
    fi
    cd $RUNDIR 
    ./wrf.exe 2>&1
}


######################################################################
# Zip WRF timeseries files
#
#    station.d??.* -> station.d??.zip
#
# Required env:
#    NDOMS, STATIONS, RUNDIR
######################################################################
function zip_ts {
    if [[ -z "$NDOMS" || -z "$RUNDIR" || -z "$STATIONS" ]]; then
        printf "$0 [$LINENO]: NDOMS, RUNDIR, or STATIONS not set, aborting\n"
        exit 1;
    fi

    OLDPWD=$PWD
    cd $RUNDIR

    for d in `seq -f '%02.0f' 1 $NDOMS`; do
        for s in $STATIONS; do

            # if there are ts files for this station and domain to process,
            FILES=`ls_tslist $s $d`
            if [ ! -z "$FILES" ]; then

                ZIPPED=$RUNDIR/${s}.d${d}.zip
                WORKING=$RUNDIR/working.${s}.d${d}.zip

                # and there is no zip file yet
                if [ ! -f "$ZIPPED" ]; then

                    # ... remove old archive to be sure,
                    rm -f "$WORKING"

                    # ... zip,
                    zip "$WORKING" $FILES > /dev/null

                    # ... clean up
                    rm -f $FILES

                    # ... and copy archive
                    mv "$WORKING" "$ZIPPED"
                fi
            fi
            
        done
    done

    cd $OLDPWD
}

######################################################################
# Zip WRF log files, and remove unzipped logs on success
#
#     $LOGS: tslist namelist.input prepare_boundaries.log rsl.out|rsl.out.0000  
#
# FIXME: deps
######################################################################
function zip_log {

    OLDPWD=$PWD
    cd $RUNDIR

    FILES="tslist namelist.input prepare_boundaries.log"
    CLEANUP="prepare_boundaries.log"

    # output filename depends on parallel / serial run
    if [ -f rsl.out.0000 ]; then
        RSL="rsl.out.0000"
    elif [ -f rsl.out ]; then
        RSL="rsl.out"
    else
        printf "$0 [$LINENO]: Can't find rsl.out or rsl.out.0000, aborting\n"
        exit 1
    fi
    FILES="$RSL $FILES"
    CLEANUP="$RSL $CLEANUP"

    # check if the log files exist
    for f in $FILES; do
        if [ ! -f $f ]; then
            printf "$0 [$LINENO]: Can't find %s, aborting\n", "$f"
            exit 1
        fi
    done

    # If there is no zip file yet
    LOGS="logs_${DATESTART}.zip"
    if [ ! -f "$LOGS" ]; then

        WORKING="${LOGS}.working"

        # ... remove old files,
        rm -rf "$WORKING"

        # ... zip,
        zip "$WORKING" $FILES > /dev/null || (
            printf "$0 [$LINENO]: Can't zip log files %s, aborting\n", "$FILES"
            exit 1
        )

        # ... and copy archive
        mv "$WORKING" "$LOGS"

        rm -f $CLEANUP
    fi

    cd $OLDPWD
}

######################################################################
# Zip WRF netCDF output by converting to netCDF4:
#
#      wrfout* -> wrfout*.nc
#
# Assumes:
#    CWD is the wrf 'run' directory    
# Required env:
#    NDOMS, DATESTART
######################################################################
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
                $NC3TONC4 -o "${NCDF3}" "${TMPFL}" && mv -f "${TMPFL}" "${NCDF4}" > /dev/null
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




# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#                                           M A I N
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

ls_tslist  veenk 1

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
            "cycle")      prepare_cycle $3 ;;
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
            "all")     zip_netcdf ; zip_ts ; zip_log ;;
            "netcdf")  zip_netcdf ;;
            "ts")      zip_ts ;;
            "log")     zip_log ;;
            *)         echo "Internal error!" ; exit 1 ;;
        esac
    ;;
    clean)
        case "$2" in
            "all")     clean_input ; clean_output ;;
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
