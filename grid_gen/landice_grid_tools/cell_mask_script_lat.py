#!/usr/bin/python
import sys, os, glob, shutil, numpy, math
import random;

from netCDF4 import *
from netCDF4 import Dataset as NetCDFFile
from pylab import *

import matplotlib
import matplotlib.pyplot as plt

from optparse import OptionParser

parser = OptionParser()
parser.add_option("-f", "--file", dest="filename", help="Path to grid file", metavar="FILE")

options, args = parser.parse_args()

if not options.filename:
	parser.error("A grid file is required.")

grid = NetCDFFile(options.filename, 'a')
#grid = NetCDFFile(options.filename, 'a', format='NETCDF3_64BIT')
#grid = NetCDFFile(options.filename, 'a', format='NETCDF3_CLASSIC')

nCells = len(grid.dimensions['nCells'])
nEdges = len(grid.dimensions['nEdges'])
nVertices = len(grid.dimensions['nVertices'])
vertexDegree = len(grid.dimensions['vertexDegree'])

latCell_full = grid.variables['latCell']
cellsOnVertex_full = grid.variables['cellsOnVertex']
#verticesOnEdge_full = grid.variables['verticesOnEdge']
#cellsOnEdge_full = grid.variables['cellsOnEdge']

print "cellsOnVertex"
print cellsOnVertex_full[:]


try:
	cullCell_full = grid.createVariable('cullCell', 'i4', ( 'nCells' ,) ) 
except:
	cullCell_full = grid.variables['cullCell']

#try:
#	vertexMask_full = grid.createVariable('vertexMask', 'f8', ( 'nVertices' ,) ) 
#except:
#	vertexMask_full = grid.variables['vertexMask']
#
#try:
#	edgeMask_full = grid.createVariable('edgeMask', 'f8', ( 'nEdges' ,) ) 
#except:
#	edgeMask_full = grid.variables['edgeMask']


removed = 0
for i in arange(0, nCells):
	#if(abs(latCell_full[i]) < 5.0 * 3.14159 / 180.0 or abs(latCell_full[i]) > 10.0 * 3.14159 /180.0 ):
	if( latCell_full[i] > -60.0 * 3.14159 / 180.0):
		cullCell_full[i] = 1
		removed = removed + 1
	else:
		cullCell_full[i] = 0

print "Removed %d cells"%removed

#removed = 0
#for i in arange(0, nVertices):
#	keep = False
#	for j in arange(0, vertexDegree):
#		if cellMask_full[cellsOnVertex_full[i,j]-1] == 1:
#			keep =True
#	
#	if keep:
#		vertexMask_full[i] = 1
#	else:
#		vertexMask_full[i] = 0
#		removed = removed+1
#
#print "Removed %d vertices"%removed
#
#removed = 0
#for i in arange(0, nEdges):
#	keep = vertexMask_full[verticesOnEdge_full[i,0]-1] == 1 and vertexMask_full[verticesOnEdge_full[i,1]-1] == 1
#	cell1 = cellsOnEdge_full[i, 0]-1
#	cell2 = cellsOnEdge_full[i, 1]-1
#
#	keep = cellMask_full[cell1] == 1
#	if cell2 != -1:
#		keep = keep or cellMask_full[cell2] == 1
#
#	if keep:
#		edgeMask_full[i] = 1
#	else:
#		edgeMask_full[i] = 0
#		removed = removed+1

#print "Removed %d edges"%removed

print "cellsOnVertex"
print cellsOnVertex_full[:]

grid.close()
