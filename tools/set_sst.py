#/share/apps/python-2.7.3/bin/python

import netCDF4 as cdf

newval = 23.5 + 272.15  # in Kelvin

files=[ 'wrfinput_d03' , 'wrfinput_d04' ]

for f in files:
    print "Setting SST of ", f, " to ", newval, "K "
    d = cdf.Dataset( f, 'r+', clobber=False )
    v = d.variables['SST']
    v[0,:,:] = newval
    d.close()
