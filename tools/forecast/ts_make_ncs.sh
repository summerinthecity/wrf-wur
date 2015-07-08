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

ncgen -k 4 -o "${1}.TS.nc" << EOF
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


        // *.d??.TS 

        float T2m(time,station) ;
            T2m:units = "K" ;
            T2m:standard_name = "t2m" ;
            T2m:long_name = "2 m Temperature" ;
            T2m:least_significant_digit = 2 ;
            T2m:_DeflateLevel = 3;
            T2m:_Shuffle = 1;

        float Q2m(time,station) ;
            Q2m:untis = "kg/kg" ;
            Q2m:standard_name = "q2m" ;
            Q2m:long_name = "2 m Vapor mixing ratio" ;
            Q2m:_DeflateLevel = 3;
            Q2m:_Shuffle = 1;

        float U10m(time, station ) ; 
            U10m:untis = "m/s" ;
            U10m:standard_name = "u10m" ;
            U10m:long_name = "10 m U wind (earth-relative)" ;
            U10m:least_significant_digit = 2 ;
            U10m:_DeflateLevel = 3;
            U10m:_Shuffle = 1;

        float V10m(time, station ) ; 
            V10m:untis = "m/s" ;
            V10m:standard_name = "v10m" ; 
            V10m:long_name = "10 m V wind (earth-relative)" ;
            V10m:least_significant_digit = 2 ;
            V10m:_DeflateLevel = 3;
            V10m:_Shuffle = 1;
            
        float psfc(time,station) ;
            psfc:units = "Pa" ; 
            psfc:add_offset = 101325.0 ;
            psfc:standard_name = "psfc" ;
            psfc:long_name = "surface pressure" ;
            psfc:least_significant_digit = 2 ;
            psfc:_DeflateLevel = 3;
            psfc:_Shuffle = 1;

        float glw(time,station) ;
            glw:units = "W/m2" ;
            glw:standard_name = "glw" ;
            glw:long_name = "downward longwave radiation flux at the ground, downward is positive" ;
            glw:least_significant_digit = 2 ;
            glw:_DeflateLevel = 3;
            glw:_Shuffle = 1;

        float gsw(time,station) ;
            gsw:units = "W/m2" ;
            gsw:standard_name = "gsw" ;
            gsw:long_name = "net shortwave radiation flux at the ground, downward is positive" ;
            gsw:least_significant_digit = 2 ;
            gsw:_DeflateLevel = 3;
            gsw:_Shuffle = 1;

        float hfx(time, station ) ; 
            hfx:units = "W/m2" ;
            hfx:standard_name = "hfx" ;
            hfx:long_name = "surface sensible heat flux, upward is positive" ;
            hfx:least_significant_digit = 2 ;
            hfx:_DeflateLevel = 3;
            hfx:_Shuffle = 1;

        float lh(time,station); 
            lh:units = "W/m2" ;
            lh:standard_name = "lh" ;
            lh:long_name = "surface latent heat flux, upward is positive" ;
            lh:least_significant_digit = 2 ;
            lh:_DeflateLevel = 3;
            lh:_Shuffle = 1;

        float tsk(time,station) ;
            tsk:units = "K" ; 
            tsk:standard_name = "tsk" ;
            tsk:long_name = "skin temperature" ;
            tsk:least_significant_digit = 2 ;
            tsk:_DeflateLevel = 3;
            tsk:_Shuffle = 1;

        float tslb1(time,station) ;
            tslb1:units = "K" ; 
            tslb1:standard_name = "tslb1" ; 
            tslb1:long_name = "top soil layer temperature" ;
            tslb1:least_significant_digit = 2 ;
            tslb1:_DeflateLevel = 3;
            tslb1:_Shuffle = 1;

        float rainc(time,station) ;
            rainc:units = "mm" ; 
            rainc:standard_name = "rainc" ;
            rainc:long_name = "rainfall from a cumulus scheme" ;
            rainc:least_significant_digit = 2 ;
            rainc:_DeflateLevel = 3;
            rainc:_Shuffle = 1;

        float rainnc(time,station) ;
            rainnc:units = "mm" ; 
            rainnc:standard_name = "rainnc" ;
            rainnc:long_name = "rainfall from an explicit scheme" ;
            rainnc:least_significant_digit = 2 ;
            rainnc:_DeflateLevel = 3;
            rainnc:_Shuffle = 1;

        float clw(time,station) ;
            clw:units = "kg/m2" ; 
            clw:standard_name = "clw" ;
            clw:long_name = "total column-integrated water vapor and cloud variables" ;
            clw:least_significant_digit = 2 ;
            clw:_DeflateLevel = 3;
            clw:_Shuffle = 1;

        float tp2m(time,station) ;
            tp2m:units = "K" ; 
            tp2m:standard_name = "tp2m" ; 
            tp2m:long_name = "Park 2m temperature from urban module" ;
            tp2m:least_significant_digit = 2 ;
            tp2m:_DeflateLevel = 3;
            tp2m:_Shuffle = 1;

        float tc2m(time,station) ;
            tc2m:units = "K" ; 
            tc2m:standard_name = "tc2m" ; 
            tc2m:long_name = "Canyon 2m temperature from urban module" ;
            tc2m:least_significant_digit = 2 ;
            tc2m:_DeflateLevel = 3;
            tc2m:_Shuffle = 1;

// global attributes:
                :Conventions = "CF-1.5" ;
                :ModelName = "WRF3.5 WUR-Urban" ;
                :_Format = "netCDF-4 classic model" ;

}
EOF


ncgen -k 4 -o "$1.UU.nc" << EOF
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


        float UU(time, station, level) ;                                   // from *.d??.UU
            UU:units = "m/s" ;
            UU:standard_name = "UU";         
            UU:long_name = "U component of wind";         
            UU:_DeflateLevel = 3;
            UU:_Shuffle = 1;

// global attributes:
                :Conventions = "CF-1.5" ;
                :ModelName = "WRF3.5 WUR-Urban" ;
                :_Format = "netCDF-4 classic model" ;

}
EOF


ncgen -k 4 -o "$1.VV.nc" << EOF
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


        float VV(time, station, level) ;                                   // from *.d??.VV
            VV:units = "m/s" ;
            VV:standard_name = "vv";         
            VV:long_name = "V component of wind";         
            VV:_DeflateLevel = 3;
            VV:_Shuffle = 1;

// global attributes:
                :Conventions = "CF-1.5" ;
                :ModelName = "WRF3.5 WUR-Urban" ;
                :_Format = "netCDF-4 classic model" ;

}
EOF


ncgen -k 4 -o "$1.TH.nc" << EOF
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


        float TH(time, station, level) ;                                   // from *.d??.TH
            TH:units = "K" ;
            TH:standard_name = "theta";         
            TH:long_name = "potential temperature";         
            TH:_DeflateLevel = 3;
            TH:_Shuffle = 1;

// global attributes:
                :Conventions = "CF-1.5" ;
                :ModelName = "WRF3.5 WUR-Urban" ;
                :_Format = "netCDF-4 classic model" ;

}
EOF


ncgen -k 4 -o "$1.QV.nc" << EOF
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


        float QV(time, station, level) ;                                   // from *.d??.QV
            QV:units = "k/kg" ;
            QV:standard_name = "qv";         
            QV:long_name = "specific humidity";         
            QV:_DeflateLevel = 3;
            QV:_Shuffle = 1;

// global attributes:
                :Conventions = "CF-1.5" ;
                :ModelName = "WRF3.5 WUR-Urban" ;
                :_Format = "netCDF-4 classic model" ;

}
EOF



ncgen -k 4 -o "$1.PH.nc" << EOF
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

        float PH(time, station, level) ;                               // from *.d??.PH
            PH:units = "m" ;
            PH:standard_name = "height" ;
            PH:long_name = "height above ground surface" ;

// global attributes:
                :Conventions = "CF-1.5" ;
                :ModelName = "WRF3.5 WUR-Urban" ;
                :_Format = "netCDF-4 classic model" ;

}
EOF
