#!/usr/bin/python
# vim: set fileencoding=utf-8 :

import f90nml
import pyproj
import argparse
import sys

#
#      500m
#  X . . . . X . . . . X . . . . X   e_we = 4
#
#      100m
#  Y Y Y Y Y Y Y Y Y Y Y Y Y Y Y Y   e_we = 5 * 3 + 1 (refinement * 3 parent grid points + 1 )
#

def parsenl( namelist ):
    """Parse a WRF namelist, typically namelist.wps, for the grid description
    Returns:
    arrays of west, east, north, south coordinates per domain in degrees, 
    the domain grid spacing dx, dy in meters,
    and the grid projection function (PROJ)"""

    geogrid = namelist['geogrid']
    ndoms = namelist['share']['max_dom']

    west = [0] * ndoms
    east = [0] * ndoms
    north = [0] * ndoms
    south = [0] * ndoms
    dx = [0] * ndoms
    dy = [0] * ndoms

    # Get a projection object from the namelist    
    # lambert conformal conical only supported projection for now
    assert namelist['geogrid']['map_proj'] == 'lambert'

    projstring  = "+proj=lcc +lat_1={truelat1} +lat_2={truelat2} +lat_0={ref_lat} +lon_0={ref_lon} ".format( **geogrid )
    projstring += "+x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
    lambert = pyproj.Proj( projstring )

    # iterate over domains
    for d in range( 0, ndoms ):
        # Get the coordinates of the corners in projected coordinates
        if d==0:
            dx[0]    =   geogrid['dx']
            dy[0]    =   geogrid['dy']
            west[0]  = - dx[0] * (geogrid['e_we'][0] - 1) * 0.5
            east[0]  =   dx[0] * (geogrid['e_we'][0] - 1) * 0.5
            north[0] =   dy[0] * (geogrid['e_sn'][0] - 1) * 0.5
            south[0] = - dy[0] * (geogrid['e_sn'][0] - 1) * 0.5
        else:
            p = geogrid['parent_id'][d] - 1 # fortran to c indexing
            west[d]  = west[p]  + ( geogrid['i_parent_start'][d] - 1 ) * dx[p]
            south[d] = south[p] + ( geogrid['j_parent_start'][d] - 1 ) * dy[p]

            dx[d] = dx[p] / ( 1.0 * geogrid['parent_grid_ratio'][d] )
            dy[d] = dy[p] / ( 1.0 * geogrid['parent_grid_ratio'][d] )

            east[d]  = west[d]  + (geogrid['e_we'][d] - 1.0) * dx[d]
            north[d] = south[d] + (geogrid['e_sn'][d] - 1.0) * dy[d]

    return west, east, north, south, dx, dy, lambert

def printgrids(namelist):
    """Print the domains as defined in the WRF namelist, listing domain extend and center"""
    west, east, north, south, dx, dy, projection = parsenl( namelist )

    geogrid = namelist['geogrid']
    for d in range(0, namelist['share']['max_dom']):
        print
        print "--------------- Domain {:2} ({:2})---------- ".format(d + 1, geogrid['parent_id'][d])
        print "Parent start     {:>12.0f} {:>12.0f}".format( geogrid['i_parent_start'][d], geogrid['j_parent_start'][d] )

        ee = geogrid['i_parent_start'][d] + (geogrid['e_we'][d] - 1) / ( geogrid['parent_grid_ratio'][d])
        en = geogrid['j_parent_start'][d] + (geogrid['e_sn'][d] - 1) / ( geogrid['parent_grid_ratio'][d])

        print "Parent end       {:>12.0f} {:>12.0f}".format( ee, en )
        print "dx dy     m      {:>12.0f} {:>12.0f}".format( dx[d], dy[d] )
        print "Extent   ij      {:>12} {:>12}".format( geogrid['e_we'][d], geogrid['e_sn'][d] )

        print "NW       km      {:>12.6f} {:>12.6f}".format( north[d] * 0.001, west[d] * 0.001 )
        print "SE       km      {:>12.6f} {:>12.6f}".format( south[d] * 0.001, east[d] * 0.001 )
        print "Center   km      {:>12.6f} {:>12.6f}".format( (north[d]+south[d])*0.5, (west[d]+east[d])*0.5 )
        print "Extent   km      {:>12.2f} {:>12.2f}".format( dx[d] * (geogrid['e_we'][d]-1) * 0.001, 
                                                             dy[d] * (geogrid['e_sn'][d]-1) * 0.001 )

        w, n = projection( west[d], north[d], inverse=True)
        e, s = projection( east[d], south[d], inverse=True)
        cx, cy = projection( (west[d]+east[d])*0.5, (north[d]+south[d])*0.5 , inverse=True )
        print "NW        °      {:>12.6f} {:>12.6f}".format( n, w )
        print "SE        °      {:>12.6f} {:>12.6f}".format( s, e )
        print "Center    °      {:>12.6f} {:>12.6f}".format( cy, cx )
        print "Extent    °      {:>12.2f} {:>12.2f}".format( n - s , e - w )

def appendnest(namelist, parent_id, parent_grid_ratio, i_parent_start, j_parent_start, e_we, e_sn):
    """Add a nested grid to the namelist using given parameters:
    parent_id (starts counting a 1, fortran style), parent_grid_ratio, i_parent_start, j_parent_start,
    e_we, e_sn (in gridpoitns)
    NOTE: also copies the geog_data_res,start_date, and end_date from the parent grid.
    """

    share = namelist['share']
    share['max_dom'] += 1
    share['start_date'].append( share['start_date'][ parent_id - 1] )
    share['end_date'].append( share['end_date'][ parent_id - 1] )

    geogrid = namelist['geogrid']
    geogrid['parent_id'].append( parent_id )
    geogrid['parent_grid_ratio'].append( parent_grid_ratio )
    geogrid['i_parent_start'].append( i_parent_start )
    geogrid['j_parent_start'].append( j_parent_start )
    geogrid['e_we'].append( e_we )
    geogrid['e_sn'].append( e_sn )
    geogrid['geog_data_res'].append( geogrid['geog_data_res'][ parent_id - 1] )


def add_centered_nest( namelist, parent_id, parent_grid_ratio, clon, clat, sizex, sizey):
    """Add a nested grid with its centered at (clon, clat) in degrees, with an extend of sizex by sizey kilometers.
    The size is adjusted as necessary to match the parent grid"""
    
    west, east, north, south, dx, dy, projection = parsenl( namelist )

    p = parent_id - 1 # forant -> C indexing

    # Get the projected coordinates of the center point, and the corners
    x, y = projection( clon, clat )
    xs = x - 0.5 * sizex * 1000.0
    xe = x + 0.5 * sizex * 1000.0
    ys = y - 0.5 * sizey * 1000.0
    ye = y + 0.5 * sizey * 1000.0
    print "X Center: ", x, "start: ", xs, "end: ", xe, "size: ", sizex
    print "Y Center: ", y, "start: ", ys, "end: ", ye, "size: ", sizey

    # Translate to WRF grid coordinates, including the fortran index offset
    starti = int( round( ( (xs - west[p])  / dx[p] ) ) ) + 1
    startj = int( round( ( (ys - south[p]) / dy[p] ) ) ) + 1
    endi   = int( round( ( (xe - west[p])  / dx[p] ) ) ) + 1
    endj   = int( round( ( (ye - south[p]) / dy[p] ) ) ) + 1

    # WRF requirement on nested grids
    e_we = int( (endi - starti) * parent_grid_ratio + 1 )
    e_sn = int( (endj - startj) * parent_grid_ratio + 1 )
    print "e_we:     ", e_we
    print "e_we:     ", e_sn

    appendnest( namelist, parent_id, parent_grid_ratio, starti, startj, e_we, e_sn)

def add_rectangular_nest( namelist, parent_id, parent_grid_ratio, newn, neww, news, newe ):
    """Add a nested grid with corners (newn, neww) and (news, newe) in degrees
    The size is adjusted as necessary to match the parent grid"""
    
    west, east, north, south, dx, dy, projection = parsenl( namelist )

    p = parent_id - 1 # forant -> C indexing

    # Get the projected coordinates of the center point, and the corners
    xe,ys = projection( newe, news )
    xs,ye = projection( neww, newn )

    # Translate to WRF grid coordinates, including the fortran index offset
    starti = int( round( ( (xs - west[p])  / dx[p] ) ) ) + 1
    startj = int( round( ( (ys - south[p]) / dy[p] ) ) ) + 1
    endi   = int( round( ( (xe - west[p])  / dx[p] ) ) ) + 1
    endj   = int( round( ( (ye - south[p]) / dy[p] ) ) ) + 1

    print "i_parent_start: ", starti
    print "j_parent_start: ", startj

    # WRF requirement on nested grids
    e_we = int( (endi - starti) * parent_grid_ratio + 1 )
    e_sn = int( (endj - startj) * parent_grid_ratio + 1 )
    print "e_we:           ", e_we
    print "e_sn:           ", e_sn

    appendnest( namelist, parent_id, parent_grid_ratio, starti, startj, e_we, e_sn)

def main():
    parser = argparse.ArgumentParser(description="Add a nested grid to an existing WRF namelist")
    parser.add_argument("-o", "--out", type=str, help="The output namelist, defaults to the input namelist.wps" )
    parser.add_argument("-c", "--center", help="Add a centered nested grid", nargs=2, metavar=('longitude','latitude'), default=False )
    parser.add_argument("-b", "--box", help="Add a nested grid defined by its corners", nargs=4, default=False,
                        metavar=('north','west','south','east' )  )
    parser.add_argument("-p", "--parent_id", type=int, help="The parent_id", default=1 )
    parser.add_argument("-r", "--ratio", type=int, help="The parent_grid_ratio, default is 5", default=5 )
    parser.add_argument("-x", "--sizex", type=float, help="Size of the domain in km", default=10 )
    parser.add_argument("-y", "--sizey", type=float, help="Size of the domain in km", default=10 )
    parser.add_argument('namelist', metavar="namelist",  type=str, nargs=1, help="WRF namelist containing 'share' and 'geogrid' sections")
    args = parser.parse_args()
    print args

    namelist = f90nml.read( args.namelist[0] )
    if args.center:
        if not args.out:
            args.out = args.namelist
            print args.out, "update"

        add_centered_nest( namelist, args.parent_id, args.ratio, args.center[0], args.center[1], args.sizex, args.sizey )
        namelist.write( args.out, force=True)
    elif args.box:
        if not args.out:
            args.out = args.namelist
            print args.out, "update"

        add_rectangular_nest( namelist, args.parent_id, args.ratio, args.box[0], args.box[1], args.box[2], args.box[3] )
        namelist.write( args.out, force=True)
    else:
        print args.namelist[0]
        printgrids( namelist )
        sys.exit()

if __name__ == "__main__":
    main()


