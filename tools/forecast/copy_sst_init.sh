#!/bin/bash

# errors are fatal
set -e

TMPFILE=temp.`basename $0`.$$

function usage {
        echo "`basename $0` copies urban SST field (temperature) from Rijkswaterstaat"
        echo "observations, as prepared by the prepare_sst.sh script, to a wrfinput file"
        echo 
        echo "Usage:"
        echo "`basename $0` sstfile wrfinput"
        echo
        echo "  sstfile    File containing the 'temperature' variable"
        echo "  wrfinput   wrfinput_d?? file"
        echo 
        echo "Run this for the hi-res domains over the Netherlands."
        exit -1
}


if [[ $# != 2 ]]; then
    echo "Got $# arguments"
    usage
fi

ncks -C -O -o "${TMPFILE}" -v temperature "$1"
ncrename -v temperature,SST "${TMPFILE}"
ncks -A -o "$2" -v SST "${TMPFILE}"
rm -f "${TMPFILE}"
