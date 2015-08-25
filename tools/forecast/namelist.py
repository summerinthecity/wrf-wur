#!/usr/bin/env python2

import f90nml
import argparse

def main(args):
    # check argparse arguments and call the appropriate function
    if args.get:
      # get namelist variable
      # verbose=True to print results to screen
      namelist_get(args.namelist[0], args.get[0], verbose=True)
    elif args.set:
      # set namelist variable
      namelist_set(args.namelist[0], args.set[0], args.set[1], verbose=True)

def namelist_get(filename, getvariable, verbose=False):
    '''
    Get the value of a variable in a namelist
      input arguments:
        filename: filename of namelist
        getvariable: GROUP_NAME:VARIABLE_NAME to get from namelist
        verbose: optional boolean argument if results should be printed to screen
    '''
    namelist = f90nml.read( filename )
    path = getvariable.split ( ':' )
    crumb = namelist
    while len(path) > 1:
        crumb = crumb[ path[0] ]
        path.pop(0)
    if isinstance(crumb, list):
        if verbose:
            print crumb[ int(path[0]) ]
        return crumb[ int(path[0]) ]
    else:
        if verbose:
            print crumb[ path[0] ]
        return crumb[ path[0] ]

def namelist_set(filename, setvariable, setvalue, verbose=False):
    '''
    Set a variable from a namelist to a value
      input arguments:
        filename: filename of namelist
        setvariable: GROUP_NAME:VARIABLE_NAME to set from namelist
        setvalue: value to set setvariable to
        verbose: optional boolean argument if results should be printed to screen
    '''
    namelist = f90nml.read( filename )
    path = setvariable.split ( ':' )
    crumb = namelist
    while len(path) > 1:
        crumb = crumb[ path[0] ]
        path.pop(0)
    # dealing with different types..
    t = type(crumb[path[0]])
    if isinstance(crumb[path[0]], int):  # integer
        crumb[ path[0] ] = int(setvalue)
    elif isinstance(crumb[path[0]], float):  # float
        crumb[ path[0] ] = float(setvalue)
    elif isinstance(crumb[path[0]], str):  # string
        crumb[ path[0] ] = setvalue
    elif isinstance(crumb[path[0]], list):  # list
        t = type( crumb[path[0]][0] )

        # deal with trailing ',' leading to empty string, crashing int() and float()
        l = setvalue.split(',')
        while l[-1] == "":
            l.pop()

        if isinstance(crumb[path[0]][0], int):  # integer
            crumb[ path[0] ] = [int(i) for i in l]
        if isinstance(crumb[path[0]][0], float):  # float
            crumb[ path[0] ] = [float(i) for i in l]
        if isinstance(crumb[path[0]][0], str):  # string
            crumb[ path[0] ] = l
    elif isinstance(crumb[path[0]], bool):  # boolean
        if setvalue == '.true.':
            crumb[ path[0] ] = True
        elif setvalue == '.false.':
            crumb[ path[0] ] = False
        else:
            print "Cannot parse boolean, use .true. or .false."
    else:
        print "Unsupported type: ", t
    # write namelist variable
    f90nml.write( namelist, filename, force=True )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="A commandline tool to read and write Fortran 90 namelist.")
    parser.add_argument('namelist', metavar="namelist",  type=str, nargs=1, help="Namelist to parse")
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-s','--set', metavar=("key","value",), required=False, type=str, nargs=2, help="Set namelist variable")
    group.add_argument('-g','--get', metavar=("key",),         required=False, type=str, nargs=1, help="Get namelist variable")
    args = parser.parse_args()
    # get/set namelist attribute
    main(args)

