#!/usr/bin/env python

import netCDF4 as cdf
import numpy   as np
import os
import re
import logging
import argparse

logging.basicConfig(level=logging.INFO)

ncfile = None
nstations = None
ntimes = None
nlevels = 15

time    = None
T2m     = None
Q2m     = None
U10m    = None
V10m    = None
psfc    = None
glw     = None
gsw     = None
hfx     = None
lh      = None
tsk     = None
tslb1   = None
rainc   = None
rainnc  = None
clw     = None
tc2m    = None
tp2m    = None
profile = None


def main():
    global ncfile, prefix

    parser = argparse.ArgumentParser(description="A commandline tool to convert WRF timeseries files to netCDF4")
    parser.add_argument('netcdf', metavar="netCDF4 file", type=str, nargs=1, help="NetCDF file base name. Create it first with ts_make_ncs.sh")
    parser.add_argument('-d', '--domain', metavar="domain", type=int, nargs=1, help="WRF domain number", default=0)
    args = parser.parse_args()

    for varname in ['TS']: # ,'UU','VV','TH','QV','PH']:
        logging.info( "{}".format(varname) )
        ncfile = cdf.Dataset( args.netcdf[0] + "." + varname + ".nc", "r+" )
        do_tslist()
        prefix  = ncfile.variables['prefix']

        for stationi in range(nstations):

            filename = "{}.d{:02d}.".format( cdf.chartostring( prefix[stationi] ), args.domain[0] )
            if not os.path.isfile( filename + "TS" ):
                continue 

            if varname in ['TS',]:
                if not ntimes:
                    simplecount(filename + "TS", args.domain[0])

                do_tsfile ( filename + "TS", stationi )
            else:
                do_profile( filename + varname, stationi, varname )

        if varname in ['TS',]:
            flush_tsfile()
        else:
            ncfile.variables[varname][:,:,:]     = profile   [:,:,:]

        ncfile.close()


def simplecount(filename, domain=-1):
    global    ntimes , time   , T2m    , Q2m    , U10m   , V10m   , psfc   , glw    , gsw    , hfx    , lh     , tsk    , tslb1  , rainc  , rainnc , clw    , tc2m   , tp2m   , profile 

    # lines = 0
    # for line in open(filename):
    #     lines += 1
    # ntimes = lines - 1 # remove header

    # override ntimes
    # d4: ntimes = 360000
    # d3: ntimes = 72000
    # d2: ntimes = 14400
    # d1: ntimes = 2880
    override = [-1, 2880, 14400, 72000, 360000]
    ntimes = override[ int(domain) ]

    time    = np.zeros([ntimes])            
    T2m     = np.zeros([ntimes,nstations])
    Q2m     = np.zeros([ntimes,nstations])
    U10m    = np.zeros([ntimes,nstations])
    V10m    = np.zeros([ntimes,nstations])
    psfc    = np.zeros([ntimes,nstations])
    glw     = np.zeros([ntimes,nstations])
    gsw     = np.zeros([ntimes,nstations])
    hfx     = np.zeros([ntimes,nstations])
    lh      = np.zeros([ntimes,nstations])
    tsk     = np.zeros([ntimes,nstations])
    tslb1   = np.zeros([ntimes,nstations])
    rainc   = np.zeros([ntimes,nstations])
    rainnc  = np.zeros([ntimes,nstations])
    clw     = np.zeros([ntimes,nstations])
    tc2m    = np.zeros([ntimes,nstations])
    tp2m    = np.zeros([ntimes,nstations])
    profile = np.zeros([ntimes,nstations, nlevels])

    logging.info( "Number of times: %s" % ntimes )

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

def do_profile(filename, stationi, varname):
    if not os.path.isfile(filename):
        logging.info( "Skipping station: %s : %s", filename, varname )
        return

    logging.debug( "{} {}".format(filename, varname) )

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
        profile[timei,stationi,:] = fields[1:]

    fileUU.close()


def do_tsfile(filename, stationi):

    if not os.path.isfile( filename):
        logging.info( "Skipping station %s : TS", cdf.chartostring( prefix[stationi] ) )
        return

    logging.debug( "{} starting".format( filename ) )

    logging.debug( "parsing ts file" )


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
        try:
            T2m   [ timei,stationi ] = fields[ 5]
            Q2m   [ timei,stationi ] = fields[ 6]
            U10m  [ timei,stationi ] = fields[ 7]
            V10m  [ timei,stationi ] = fields[ 8]
            psfc  [ timei,stationi ] = fields[ 9]
            glw   [ timei,stationi ] = fields[10]
            gsw   [ timei,stationi ] = fields[11]
            hfx   [ timei,stationi ] = fields[12]
            lh    [ timei,stationi ] = fields[13]
            tsk   [ timei,stationi ] = fields[14]
            tslb1 [ timei,stationi ] = fields[15]
            rainc [ timei,stationi ] = fields[16]
            rainnc[ timei,stationi ] = fields[17]
            clw   [ timei,stationi ] = fields[18]
            tc2m  [ timei,stationi ] = fields[19]
            tp2m  [ timei,stationi ] = fields[20]
        except:
            print "Error parsing for time: ", timei

    fileTS.close()


    logging.info( "{} done".format( filename ) )

def flush_tsfile():
    logging.debug( "writing to netcdf file" )
    ncfile.variables['time'][:]      = time  [:]
    ncfile.variables['T2m'][:,:]     = T2m   [:]
    ncfile.variables['Q2m'][:,:]     = Q2m   [:]
    ncfile.variables['U10m'][:,:]    = U10m  [:]
    ncfile.variables['V10m'][:,:]    = V10m  [:]
    ncfile.variables['psfc'][:,:]    = psfc  [:]  - ncfile.variables['psfc'].add_offset
    ncfile.variables['glw'][:,:]     = glw   [:]
    ncfile.variables['gsw'][:,:]     = gsw   [:]
    ncfile.variables['hfx'][:,:]     = hfx   [:]
    ncfile.variables['lh'][:,:]      = lh    [:]
    ncfile.variables['tsk'][:,:]     = tsk   [:]
    ncfile.variables['tslb1'][:,:]   = tslb1 [:]
    ncfile.variables['rainc'][:,:]   = rainc [:]
    ncfile.variables['rainnc'][:,:]  = rainnc[:]
    ncfile.variables['clw'][:,:]     = clw   [:]
    ncfile.variables['tc2m'][:,:]    = tc2m  [:]
    ncfile.variables['tp2m'][:,:]    = tp2m  [:]

if __name__ == "__main__":
    main()
   
