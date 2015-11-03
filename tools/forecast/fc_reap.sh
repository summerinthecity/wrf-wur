#!/bin/sh

TOP=/home/jattema/WRF/WRFV3

cd $TOP
RUNS="2015-??-??"

for R in $RUNS; do
    export RUNDIR=$TOP/$R

    OUT=${RUNDIR}/rsl.out.0000
    if [ ! -e "${OUT}" ]; then
        echo "Skipping run ${RUN} (not run or already archived)"
        continue
    fi

    if grep -q ' wrf: SUCCESS COMPLETE WRF' ${OUT}; then
        echo "Run: ${R} finished"

        JOB=$RUNDIR/post.job
        cat $TOP/tools/forecast/post.template | sed "s=%RUNDIR%=$RUNDIR=" > "$JOB"
        sbatch "$JOB"
    else
        echo "Skipping run ${R} (not finished)"
    fi
done
