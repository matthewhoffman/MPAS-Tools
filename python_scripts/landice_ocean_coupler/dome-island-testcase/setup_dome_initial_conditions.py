#!/usr/bin/python
# Generate initial conditions for dome land ice test case

import sys, numpy
try:
  from Scientific.IO.NetCDF import NetCDFFile
  netCDF_module = 'Scientific.IO.NetCDF'
except ImportError:
  try:
    from netCDF4 import Dataset as NetCDFFile
    netCDF_module = 'netCDF4'
  except ImportError:
      print 'Unable to import any of the following python modules:'
      print '  Scientific.IO.NetCDF \n  netcdf4 '
      print 'One of them must be installed.'
      raise ImportError('No netCDF module found')
from math import sqrt

# Check to see if a grid file was specified on the command line.
# If not, land_ice_grid.nc is used.
if len(sys.argv) > 1:
  if sys.argv[1][0] == '-': # The filename can't begin with a hyphen
    print '\nUsage:  python setup_dome_initial_conditions.py [GRID.NC]\nIf no filename is supplied, land_ice_grid.nc will be used.'
    sys.exit(0)
  else:
    gridfilename = sys.argv[1]
else:
  gridfilename = 'land_ice_grid.nc'



# Open the file, get needed dimensions
try:
    gridfile = NetCDFFile(gridfilename,'r+')
    if (netCDF_module == 'Scientific.IO.NetCDF'):
         nVertLevels = gridfile.dimensions['nVertLevels']
    else:
         nVertLevels = len(gridfile.dimensions['nVertLevels'])
    if nVertLevels != 9:
         print 'nVerLevels in the supplied file was ', nVertLevels, '.  Were you expecting 9?'
    # Get variables
    xCell = gridfile.variables['xCell'][:]
    yCell = gridfile.variables['yCell'][:]
    thickness = gridfile.variables['thickness'][:]
    bedTopography = gridfile.variables['bedTopography'][:]
    normalVelocity = gridfile.variables['normalVelocity'][:]
    layerThicknessFractions = gridfile.variables['layerThicknessFractions'][:]
    temperature = gridfile.variables['temperature'][:]
    # Get b.c. variables
    beta = gridfile.variables['betaTimeSeries'][:]
    SMB = gridfile.variables['sfcMassBalTimeSeries'][:]
    Tsfc = gridfile.variables['sfcAirTempTimeSeries'][:]
    G = gridfile.variables['basalHeatFluxTimeSeries'][:]
    BMB = gridfile.variables['marineBasalMassBalTimeSeries'][:]
except:
    sys.exit('Error: The grid file specified is either missing or lacking needed dimensions/variables.')


# Assign variable values for dome
# Define dome dimensions - all in meters
r0 = 100000.0 * sqrt(0.125)
h0 = 5000.0 * sqrt(0.125)
#x0 = 30000.0
#y0 = 30000.0
x0 = xCell.min() + 0.5 * (xCell.max() - xCell.min() )
y0 = yCell.min() + 0.5 * (yCell.max() - yCell.min() )
# Calculate distance of each cell center from dome center
r = ((xCell - x0)**2 + (yCell - y0)**2)**0.5
# Set default value for non-dome cells
thickness[:] = 0.0
# Calculate the dome thickness for cells within the desired radius (thickness will be NaN otherwise)
thickness[0, r<r0] = h0 * (1.0 - (r[r<r0] / r0)**2)**0.5
# zero velocity everywhere
normalVelocity[:] = 0.0
# flat bed at sea level
bedTopography[:] = -1000.0
# constant, arbitrary temperature, degrees C
temperature[:] = -20.0 
# Setup layerThicknessFractions
layerThicknessFractions[:] = 1.0 / nVertLevels

# boundary conditions
# Sample values to use, or comment these out for them to be 0.
beta[:] = 10000.
#SMB[:] = 2.0/1000.0 * (thickness[:] + bedTopography[:]) - 1.0  # units: m/yr, lapse rate of 1 m/yr with 0 at 500 m
Tsfc[:,0] = -5.0/1000.0 * (thickness[0,:] + bedTopography[0,:]) # lapse rate of 5 deg / km
G = 0.01
#BMB[:] = -20.0  # units: m/yr

# Reassign the modified numpy array values back into netCDF variable objects 
gridfile.variables['thickness'][:] = thickness
gridfile.variables['normalVelocity'][:] = normalVelocity
gridfile.variables['bedTopography'][:] = bedTopography
gridfile.variables['temperature'][:] = temperature
gridfile.variables['layerThicknessFractions'][:] = layerThicknessFractions
gridfile.variables['betaTimeSeries'][:] = beta
gridfile.variables['sfcMassBalTimeSeries'][:] = SMB
gridfile.variables['sfcAirTempTimeSeries'][:] = Tsfc
gridfile.variables['basalHeatFluxTimeSeries'][:] = G
gridfile.variables['marineBasalMassBalTimeSeries'][:] = BMB

gridfile.close()
print 'Successfully added dome initial conditions to ', gridfilename

