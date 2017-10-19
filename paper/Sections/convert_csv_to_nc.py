#!/usr/bin/env python

    ## LOTS of assumptions about how you want to store the data in netcdf!
    ## I'm assuming that each colum nis a separate array here
    ## and that it's a time series, but it could be anything
    ## a Nx# single array, for instance..
    ## take a look at: http://cfconventions.org/latest.html
    ## to see how you should store netcdf data


import numpy as np
import netCDF4
## docs for netcdf lib here: 
#http://netcdf4-python.googlecode.com/svn/trunk/docs/netCDF4-module.html
# load the data

# this load the file into a Nx3 array (three columns)
data = np.loadtxt('beta_iter4.csv', delimiter=',')


# create a netcdf Data object

with netCDF4.Dataset('beta_iter4.nc', mode="w", format='NETCDF4') as ds:
    # some file-level meta-data attributes:
    ds.Conventions = "CF-1.6" # if you comply with the convension -- which you should!
    ds.title = 'Beta and basal friction from Ymir for 4th (non-converged) iteration'
    ds.institution = 'UTIG'
    ds.source = 'Output from Ymir 0.1.1478-26c3 with bedmap2, shapiro&ritzwoller, greater Thwaites domain'
    ds.history = 'Converted to netcdf from csv output from paraview from ant_beta.face0.pvtu Ymir output'

    # defining the dimensions of your arrays:
    time = ds.createDimension('time', data.shape[0])

    # variables for the columns -- you should use real names
    names = ['beta','basal friction (Pa/(km/a))','Easting','Northing','z']
    for i in range(data.shape[4]):
        print(names[i])
        var = ds.createVariable('var%i'%i,data.dtype, (names[i],))
        var[:] = data[:,i] 
        ## adds some attibutes
        var.units = 'the_proper_unit_string'
        var.long_name = 'a nice long name that describes the data'
        var.standard_name = 'a_CF_standard_name'
    print ds