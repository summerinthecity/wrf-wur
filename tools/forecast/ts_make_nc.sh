#!/bin/bash

HELP="""
Create a netCDF4 file to hold the WRF timeseries output data.
The 'tslist' file in the current directory is parsed to find the number of stations.

Usage: ts_make_nc.sh <filename>
"""


if [ $# != 1 ]; then
    echo $HELP
fi

# get the number of stations
L=`cat tslist | wc -l`
NSTATIONS=$(( L - 3 ))
echo "Found $NSTATIONS in tslist"

ncgen -k 4 -o $1 << EOF
netcdf timeseries {

dimensions:
        time = UNLIMITED ;  // (10 currently)
        station = $NSTATIONS ;
        level = 15 ;
        strln = 50 ; 

variables:
        float time(time) ; 
            time:calendar = "standard" ;
            time:units = "hours since start of run" ;
            time:least_significant_digit = 3 ;

        int level(level) ;
            level:standard_name = "level" ;
            level:long_name = "level" ;

        int station(station) ;
            station:standard_name = "statoin" ;
            station:long_name = "stadion id" ; 

        char name(station,strln) ;                                         // from tslist
            name:standard_name = "name" ; 
            name:long_name = "Station name"  ;

        float lat(station) ;                                               // from tslist
            lat:standard_name = "latitude" ; 
            lat:long_name = "Station latitude"  ;

        float lon(station) ;                                               // from tslist
            lon:standard_name = "longitude" ; 
            lon:long_name = "Station longitude"  ;

        float glat(station) ;                                              // from tslist
            glat:standard_name = "latitude" ; 
            glat:long_name = "Actual station latitude"  ;

        float glon(station) ;                                              // from tslist
            glon:standard_name = "longitude" ; 
            glon:long_name = "Actual station longitude"  ;

        float gi(station) ;                                                // from tslist
            gi:standard_name = "gi" ; 
            gi:long_name = "Actual station grid i-index"  ;

        float gj(station) ;                                                // from tslist
            gj:standard_name = "gj" ; 
            gj:long_name = "Actual station grid j-index"  ;

        float elevation(station) ;                                         // from tslist
            elevation:standard_name = "elevation" ;
            elevation:long_name = "Station height above sea level" ;

        char prefix(station,strln) ;                                       // from tslist
            prefix:standard_name = "pfx" ;
            prefix:long_name = "prefix" ;
        
        int gridi(station) ;                                               // from tslist
            gridi:standard_name = "i" ;

        int gridj(station) ;                                               // from tslist
            gridj:standard_name = "j" ;


        float height(time, station, level) ;                               // from *.d??.PH
            height:units = "m" ;
            height:standard_name = "height" ;
            height:long_name = "height above ground surface" ;

        float qv(time, station, level) ;                                   // from *.d??.QV
            qv:units = "k/kg" ;
            qv:standard_name = "qv";         
            qv:long_name = "specific humidity";         
            qv:_DeflateLevel = 9;
            qv:_Shuffle = 1;

        float th(time, station, level) ;                                   // from *.d??.TH
            th:units = "K" ;
            th:standard_name = "theta";         
            th:long_name = "potential temperature";         
            th:_DeflateLevel = 9;
            th:_Shuffle = 1;

        float uu(time, station, level) ;                                   // from *.d??.UU
            uu:units = "m/s" ;
            uu:standard_name = "uu";         
            uu:long_name = "U component of wind";         
            uu:_DeflateLevel = 9;
            uu:_Shuffle = 1;

        float vv(time, station, level) ;                                   // from *.d??.VV
            vv:units = "m/s" ;
            vv:standard_name = "vv";         
            vv:long_name = "V component of wind";         
            vv:_DeflateLevel = 9;
            vv:_Shuffle = 1;

        // *.d??.TS 

        float T2m(time,station) ;
            T2m:units = "K" ;
            T2m:standard_name = "t2m" ;
            T2m:long_name = "2 m Temperature" ;
            T2m:least_significant_digit = 2 ;
            T2m:_DeflateLevel = 9;
            T2m:_Shuffle = 1;

        float Q2m(time,station) ;
            Q2m:untis = "kg/kg" ;
            Q2m:standard_name = "q2m" ;
            Q2m:long_name = "2 m Vapor mixing ratio" ;
            Q2m:_DeflateLevel = 9;
            Q2m:_Shuffle = 1;

        float U10m(time, station ) ; 
            U10m:untis = "m/s" ;
            U10m:standard_name = "u10m" ;
            U10m:long_name = "10 m U wind (earth-relative)" ;
            U10m:least_significant_digit = 2 ;
            U10m:_DeflateLevel = 9;
            U10m:_Shuffle = 1;

        float V10m(time, station ) ; 
            V10m:untis = "m/s" ;
            V10m:standard_name = "v10m" ; 
            V10m:long_name = "10 m V wind (earth-relative)" ;
            V10m:least_significant_digit = 2 ;
            V10m:_DeflateLevel = 9;
            V10m:_Shuffle = 1;
            
        float psfc(time,station) ;
            psfc:units = "Pa" ; 
            psfc:add_offset = 101325.0 ;
            psfc:standard_name = "psfc" ;
            psfc:long_name = "surface pressure" ;
            psfc:least_significant_digit = 2 ;
            psfc:_DeflateLevel = 9;
            psfc:_Shuffle = 1;

        float glw(time,station) ;
            glw:units = "W/m2" ;
            glw:standard_name = "glw" ;
            glw:long_name = "downward longwave radiation flux at the ground, downward is positive" ;
            glw:least_significant_digit = 2 ;
            glw:_DeflateLevel = 9;
            glw:_Shuffle = 1;

        float gsw(time,station) ;
            gsw:units = "W/m2" ;
            gsw:standard_name = "gsw" ;
            gsw:long_name = "net shortwave radiation flux at the ground, downward is positive" ;
            gsw:least_significant_digit = 2 ;
            gsw:_DeflateLevel = 9;
            gsw:_Shuffle = 1;

        float hfx(time, station ) ; 
            hfx:units = "W/m2" ;
            hfx:standard_name = "hfx" ;
            hfx:long_name = "surface sensible heat flux, upward is positive" ;
            hfx:least_significant_digit = 2 ;
            hfx:_DeflateLevel = 9;
            hfx:_Shuffle = 1;

        float lh(time,station); 
            lh:units = "W/m2" ;
            lh:standard_name = "lh" ;
            lh:long_name = "surface latent heat flux, upward is positive" ;
            lh:least_significant_digit = 2 ;
            lh:_DeflateLevel = 9;
            lh:_Shuffle = 1;

        float tsk(time,station) ;
            tsk:units = "K" ; 
            tsk:standard_name = "tsk" ;
            tsk:long_name = "skin temperature" ;
            tsk:least_significant_digit = 2 ;
            tsk:_DeflateLevel = 9;
            tsk:_Shuffle = 1;

        float tslb1(time,station) ;
            tslb1:units = "K" ; 
            tslb1:standard_name = "tslb1" ; 
            tslb1:long_name = "top soil layer temperature" ;
            tslb1:least_significant_digit = 2 ;
            tslb1:_DeflateLevel = 9;
            tslb1:_Shuffle = 1;

        float rainc(time,station) ;
            rainc:units = "mm" ; 
            rainc:standard_name = "rainc" ;
            rainc:long_name = "rainfall from a cumulus scheme" ;
            rainc:least_significant_digit = 2 ;
            rainc:_DeflateLevel = 9;
            rainc:_Shuffle = 1;

        float rainnc(time,station) ;
            rainnc:units = "mm" ; 
            rainnc:standard_name = "rainnc" ;
            rainnc:long_name = "rainfall from an explicit scheme" ;
            rainnc:least_significant_digit = 2 ;
            rainnc:_DeflateLevel = 9;
            rainnc:_Shuffle = 1;

        float clw(time,station) ;
            clw:units = "kg/m2" ; 
            clw:standard_name = "clw" ;
            clw:long_name = "total column-integrated water vapor and cloud variables" ;
            clw:least_significant_digit = 2 ;
            clw:_DeflateLevel = 9;
            clw:_Shuffle = 1;

        float tp2m(time,station) ;
            tp2m:units = "K" ; 
            tp2m:standard_name = "tp2m" ; 
            tp2m:long_name = "Park 2m temperature from urban module" ;
            tp2m:least_significant_digit = 2 ;
            tp2m:_DeflateLevel = 9;
            tp2m:_Shuffle = 1;

        float tc2m(time,station) ;
            tc2m:units = "K" ; 
            tc2m:standard_name = "tc2m" ; 
            tc2m:long_name = "Canyon 2m temperature from urban module" ;
            tc2m:least_significant_digit = 2 ;
            tc2m:_DeflateLevel = 9;
            tc2m:_Shuffle = 1;

// global attributes:
                :Conventions = "CF-1.5" ;
                :ModelName = "WRF3.5 WUR-Urban" ;
                :_Format = "netCDF-4 classic model" ;

}
EOF

