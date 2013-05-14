#!/usr/bin/python
import sys, os, glob, shutil, numpy, math, subprocess

from netCDF4 import *
from netCDF4 import Dataset as NetCDFFile
from pylab import *

from optparse import OptionParser

parser = OptionParser()
parser.add_option("--lice", dest="landice_grid", help="Path to land ice grid file", metavar="GRID")
parser.add_option("--ocean", dest="ocean_grid", help="Path to ocean grid file", metavar="GRID")
parser.add_option("--output", dest="output_grid", help="Path to output grid file", metavar="GRID")

options, args = parser.parse_args()

if not options.landice_grid:
	parser.error("A land ice grid file is required.")
else:
	landice_grid_path = options.landice_grid

if not options.ocean_grid:
	parser.error("A land ice grid file is required.")
else:
	ocean_grid_path = options.ocean_grid

if not options.output_grid:
	output_grid_path = "new_ocean.nc"
	print "Output grid not specified. Defaulting to %s"%output_grid_path
else:
	output_grid = options.output_grid

# ---------- PARAMETERS --------------- #
gravity = 9.8101
oceanThicknessMinimum = 50.0 # meters
densityIce = 900.0  # land ice model's ice density (kg m^-3)

subprocess.check_call('cp %s %s'%(ocean_grid_path, output_grid_path), shell=True, executable='/bin/bash');

lice_grid = NetCDFFile(landice_grid_path, 'r')
ocean_grid = NetCDFFile(output_grid_path, 'a')

lice_nCells = len(lice_grid.dimensions['nCells'])
ocean_nCells = len(ocean_grid.dimensions['nCells'])

if lice_nCells != ocean_nCells:
	print "Land ice and ocean grids have different dimensions. Exiting."
	lice_grid.close()
	ocean_grid.close()
	quit()


nCells = len(ocean_grid.dimensions['nCells'])
nVertLevels = len(ocean_grid.dimensions['nVertLevels'])

iceThickness = lice_grid.variables['thickness'][0,:]

oceanThickness_full = ocean_grid.variables['layerThickness']

oceanDensity = ocean_grid.variables['density'][0,:,:]
oceanThickness = oceanThickness_full[0,:,:]
bottomDepth = ocean_grid.variables['bottomDepth'][:]

try:
	oceanSurfacePressure = ocean_grid.createVariable('seaSurfacePressure', 'f8', ( 'nCells' ,) )
except:
	oceanSurfacePressure = ocean_grid.variables['seaSurfacePressure']

for i in arange(0, nCells):
	bottomPressure = oceanDensity[i][0] * gravity * 0.5 * oceanThickness[i][0]
	avgDensity = oceanDensity[i][0]
	for j in arange(1, nVertLevels):
		avgDensity += oceanDensity[i][j]
		bottomPressure += 0.5 * gravity * (oceanDensity[i][j-1] * oceanThickness[i][j-1] + oceanDensity[i][j] * oceanThickness[i][j])

	avgDensity /= nVertLevels
	pressureLimit = bottomPressure - gravity * oceanDensity[i][nVertLevels-1] * oceanThicknessMinimum

	surfPressure = np.minimum(densityIce * gravity * iceThickness[i], pressureLimit)

	if surfPressure > 0.0:
		columnThickness = surfPressure / ( gravity * avgDensity )
	else:
		columnThickness = bottomDepth[i]

	print 'Column Thickness = ', columnThickness, ' pressureLimit = ', pressureLimit, ' surfPressure = ', surfPressure
	ocean_grid.variables['layerThickness'][0,i,:] = columnThickness / nVertLevels
	oceanSurfacePressure[i] = surfPressure

lice_grid.close()
ocean_grid.close()
