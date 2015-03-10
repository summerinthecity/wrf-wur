#!/usr/bin/env python

import netCDF4 as cdf
import numpy   as np
import os
import re
import logging
import argparse

logging.basicConfig(level=logging.DEBUG)

ncfile = None
nstations = None
ntimes = None
nlevels = 15

def main():
    global ncfile

    parser = argparse.ArgumentParser(description="A commandline tool to convert WRF timeseries files to netCDF4")
    parser.add_argument('netcdf', metavar="netCDF4 file", type=str, nargs=1, help="NetCDF file name. Create it first with ts_make_nc.sh")
    parser.add_argument('-d', '--domain', metavar="domain", type=int, nargs=1, help="WRF domain number", default=0)
    args = parser.parse_args()

    ncfile = cdf.Dataset( args.netcdf[0], "r+" )
    do_tslist()

    prefix  = ncfile.variables['prefix']
    for stationi in range(nstations):
        filename = "{}.d{:02d}.".format( cdf.chartostring( prefix[stationi] ), args.domain[0] )

        if os.path.isfile( filename + "TS" ):
            do_tsfile ( filename + "TS", stationi )
            do_profile( filename + "UU", stationi, 'uu' )
            do_profile( filename + "VV", stationi, 'vv' )
            do_profile( filename + "TH", stationi, 'th' )
            do_profile( filename + "QV", stationi, 'qv' )
            do_profile( filename + "PH", stationi, 'height' )
        else:
            logging.info( "Skipping station %s", cdf.chartostring( prefix[stationi] ) )

    ncfile.close()


def simplecount(filename):
    lines = 0
    for line in open(filename):
        lines += 1
    return lines

def do_tslist():
    global nstations

    # Parse tslist

    station = ncfile.variables['station']
    name    = ncfile.variables['name']
    prefix  = ncfile.variables['prefix']
    lat     = ncfile.variables['lat']
    lon     = ncfile.variables['lon']

    strln = len( ncfile.dimensions['strln'] )

    filetslist = open( 'tslist', 'r' )

    # Header
    # #-----------------------------------------------#
    # # 24 characters for name | pfx |  LAT  |   LON  |
    # #-----------------------------------------------#

    filetslist.next()
    filetslist.next()
    filetslist.next()

    # Body
    # veenkampen                veenk 51.98101  5.61957
    stationi = -1
    for line in filetslist:
        stationi += 1
        fields = line.split()
        station[stationi] = stationi
        name[stationi]    = cdf.stringtoarr( fields[ 0], strln )
        prefix[stationi]  = cdf.stringtoarr( fields[ 1], strln )
        lat[stationi]     = fields[ 2]
        lon[stationi]     = fields[ 3]
    nstations = stationi + 1

    filetslist.close()
    logging.info( "tslist done, found %i stations", nstations )

def do_profile(filename, stationi, varname):

    logging.debug( "{} starting".format( filename ) )
    profile = np.zeros([ntimes, nlevels] )

    # Parse UU file

    fileUU = open( filename, 'r' )

    # Header
    # veenkampen                 1  1 veenk ( 51.981,   5.620) (  60,  60) ( 51.965,   5.663)   15.4 meters

    fields = fileUU.next().split()

    # Body
    # each line starting with the model time in hours, followed by the variable at model level 1,2,3, ... 
    # up to the highest model level of interest 

    timei = -1 
    for line in fileUU:
        timei += 1
        fields = line.split()
        profile[timei,:] = fields[1:]

    fileUU.close()

    logging.debug( "writing to netcdf file" )
    ncfile.variables[varname][:,stationi,:]     = profile   [:,:]

    logging.info( "{} done".format( filename ) )


def do_tsfile(filename, stationi):
    global ntimes

    ntimes = simplecount(filename) - 1 # remove header

    logging.debug( "{} starting".format( filename ) )

    logging.debug( "parsing ts file" )
    time    = np.zeros([ntimes])            
    T2m     = np.zeros([ntimes] )
    Q2m     = np.zeros([ntimes] )
    U10m    = np.zeros([ntimes] )
    V10m    = np.zeros([ntimes] )
    psfc    = np.zeros([ntimes] )
    glw     = np.zeros([ntimes] )
    gsw     = np.zeros([ntimes] )
    hfx     = np.zeros([ntimes] )
    lh      = np.zeros([ntimes] )
    tsk     = np.zeros([ntimes] )
    tslb1   = np.zeros([ntimes] )
    rainc   = np.zeros([ntimes] )
    rainnc  = np.zeros([ntimes] )
    clw     = np.zeros([ntimes] )
    tc2m    = np.zeros([ntimes] )
    tp2m    = np.zeros([ntimes] )


    # Parse TS file, see the README.tslist in the WRF run directory for details

    fileTS = open( filename, 'r' )

    # Header
    # NZCM McMurdo               2  7 mcm   (-77.850, 166.710) ( 153, 207) (-77.768, 166.500)   81.8 meters
    # 
    # Those are name of the station, grid ID, time-series ID, station lat/lon, grid indices (nearest grid point to
    # the station location), grid lat/lon, elevation.

    line = fileTS.next()
    fields = re.search(r"\( *(\d+) *, *(\d+) *\) *\( *(-?\d+\.\d+) *, *(-?\d+\.\d+) *\) * (-?\d+(\.\d+)?) *meters", line).groups()

    ncfile.variables['gj'  ][stationi]      = fields[0]
    ncfile.variables['gi'  ][stationi]      = fields[1]
    ncfile.variables['glat'][stationi]      = fields[2]
    ncfile.variables['glon'][stationi]      = fields[3]
    ncfile.variables['elevation'][stationi] = fields[4]

    # Body
    # 0   1           2       3   4  5  6  7  8  9     10   11   12   13  14  15        16     17      18
    # id, ts_hour, id_tsloc, ix, iy, t, q, u, v, psfc, glw, gsw, hfx, lh, tsk, tslb(1), rainc, rainnc, clw
    timei = -1 
    for line in fileTS:
        timei += 1
        fields = line.split()
        time  [ timei ] = fields[ 1]
        T2m   [ timei ] = fields[ 5]
        Q2m   [ timei ] = fields[ 6]
        U10m  [ timei ] = fields[ 7]
        V10m  [ timei ] = fields[ 8]
        psfc  [ timei ] = fields[ 9]
        glw   [ timei ] = fields[10]
        gsw   [ timei ] = fields[11]
        hfx   [ timei ] = fields[12]
        lh    [ timei ] = fields[13]
        tsk   [ timei ] = fields[14]
        tslb1 [ timei ] = fields[15]
        rainc [ timei ] = fields[16]
        rainnc[ timei ] = fields[17]
        clw   [ timei ] = fields[18]
        tc2m  [ timei ] = fields[19]
        tp2m  [ timei ] = fields[20]

    fileTS.close()

    logging.debug( "writing to netcdf file" )
    

    ncfile.variables['time'][:]             = time  [:]
    ncfile.variables['T2m'][:,stationi]     = T2m   [:]
    ncfile.variables['Q2m'][:,stationi]     = Q2m   [:]
    ncfile.variables['U10m'][:,stationi]    = U10m  [:]
    ncfile.variables['V10m'][:,stationi]    = V10m  [:]
    ncfile.variables['psfc'][:,stationi]    = psfc  [:]  - ncfile.variables['psfc'].add_offset
    ncfile.variables['glw'][:,stationi]     = glw   [:]
    ncfile.variables['gsw'][:,stationi]     = gsw   [:]
    ncfile.variables['hfx'][:,stationi]     = hfx   [:]
    ncfile.variables['lh'][:,stationi]      = lh    [:]
    ncfile.variables['tsk'][:,stationi]     = tsk   [:]
    ncfile.variables['tslb1'][:,stationi]   = tslb1 [:]
    ncfile.variables['rainc'][:,stationi]   = rainc [:]
    ncfile.variables['rainnc'][:,stationi]  = rainnc[:]
    ncfile.variables['clw'][:,stationi]     = clw   [:]
    ncfile.variables['tc2m'][:,stationi]    = tc2m  [:]
    ncfile.variables['tp2m'][:,stationi]    = tp2m  [:]

    logging.info( "{} done".format( filename ) )
 

if __name__ == "__main__":
    main()
   
