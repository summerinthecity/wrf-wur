#!/usr/bin/env python

import f90nml
import argparse
import sys


def main():
    parser = argparse.ArgumentParser(description="A commandline tool to read and write Fortran 90 namelist.")
    parser.add_argument('namelist', metavar="namelist",  type=str, nargs=1, help="Namelist to parse")
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-s','--set', metavar=("key","value",), required=False, type=str, nargs=2, help="Set namelist variable")
    group.add_argument('-g','--get', metavar=("key",),         required=False, type=str, nargs=1, help="Get namelist variable")
    # parser.add_argument('key',      metavar="key",       type=str, nargs=1, help="Name of the variable, including the section. For instance the variable MAX_DOM in the section DOMAINS would be MAX_DOM:DOMAINS")
    # parser.add_argument('value',    metavar="value",     type=str, nargs=1, help="New value.")
    args = parser.parse_args()

    namelist = f90nml.read( args.namelist[0] )

    if args.get:
        path = args.get[0].split ( ':' )
        crumb = namelist
        while len(path) > 1:
            crumb = crumb[ path[0] ]
            path.pop(0)
        if type( crumb ) == type( [] ):
            print crumb[ int(path[0]) ]
        else:
            print crumb[ path[0] ]
        
    elif args.set:
        path = args.set[0].split ( ':' )
        crumb = namelist
        while len(path) > 1:
            crumb = crumb[ path[0] ]
            path.pop(0)

        # dealing with different types..
        t = type(crumb[path[0]])
       
        if t == type( 1 ):
            crumb[ path[0] ] = int(args.set[1])
        elif t == type( 1.0 ):
            crumb[ path[0] ] = float(args.set[1])
        elif t == type( '' ):
            crumb[ path[0] ] = args.set[1]
        elif t == type( [] ):
            t = type( crumb[path[0]][0] )

            # deal with trailing ',' leading to empty string, crashing int() and float()
            l = args.set[1].split(',')
            while l[-1] == "":
                l.pop()

            if t == type( 1 ):
                crumb[ path[0] ] = [int(i) for i in l]
            if t == type( 1.0 ):
                crumb[ path[0] ] = [float(i) for i in l]
            if t == type( '' ):
                crumb[ path[0] ] = l
        elif t == type( True ):
            if args.set[1] == '.true.':
                crumb[ path[0] ] = True
            elif args.set[1] == '.false.':
                crumb[ path[0] ] = False
            else:
                print "Cannot parse boolean, use .true. or .false."
        else:
            print "Unsupported type: ", t

        f90nml.write( namelist, args.namelist[0], force=True )

if __name__ == "__main__":
    main()

