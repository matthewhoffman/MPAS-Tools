# Script to take CISM output and convert to ascii files of the format:
# latitude, longitude, thickness
# This involves reprojecting from polar stereographic to lat long
# and then building 1-d arrays to save as text files for each time slice.

import netCDF4
import numpy as np
import pyproj
#import sys


# Standard Lat/Long 
projLL   = pyproj.Proj(proj='latlong', datum='WGS84')

# Possible output projections
projGIS_Bamber = pyproj.Proj('+proj=stere +lat_ts=71.0 +lat_0=90 +lon_0=321.0 +k_0=1.0 +x_0=800000.0 +y_0=3400000.0 +geoidgrids=./egm08_25.gtx')

projAIS_BEDMAP2 = pyproj.Proj('+proj=stere +lat_ts=-71.0 +lat_0=-90 +lon_0=0.0 +k_0=1.0 +x_0=0.0 +y_0=0.0 +geoidgrids=./egm08_25.gtx')

projAIS_BEDMAP2_no_geoid = pyproj.Proj('+proj=stere +lat_ts=-71.0 +lat_0=-90 +lon_0=0.0 +k_0=1.0 +x_0=0.0 +y_0=0.0')



####################### SET THIS STUFF ###################

#infile = 'grid.1000.mpas.nc'

infile = '../culled_mesh.nc'


#outfilePrefix = 'grid.1000.mpas.plane.'
projOut = projAIS_BEDMAP2_no_geoid
saveLatMax = -00.0
saveLatMin = -91.0

##########################################


#======== Notes about pyproj =========
# The documentation for the transform function says:"In addition to converting between cartographic and geographic projection coordinates, this function can take care of datum shifts (which cannot be done using the __call__ method of the Proj instances). "
# Note this approach becomes an x,y,z transform instead of just x,y.
# ===================================


d2r = np.pi / 180.0
r2d = 180.0 / np.pi


# Get stuff from main file
f = netCDF4.Dataset(infile,'r')
lat = f.variables['latCell'][:] * r2d
lon = f.variables['lonCell'][:] * r2d
cellsOnVertex = f.variables['cellsOnVertex'][:]
nVertices = len(f.dimensions['nVertices'])
nCells = len(f.dimensions['nCells'])
meshDensity = f.variables['meshDensity'][:]



d2r = np.pi / 180.0
keepCells = np.nonzero( np.logical_and(lat<saveLatMax, lat>saveLatMin) )[0]


[xOut, yOut] = projOut(lon[keepCells], lat[keepCells])


# Build lookup table from full indices to new pruned indices
cellIndicesFullToPruned = np.zeros( ( nCells, ) )
counter = 0 # make the new indices 1-based
for c in keepCells:
    counter = counter + 1
    cellIndicesFullToPruned[c] = counter

print cellIndicesFullToPruned

##Optional - plot out the lat long positions to make sure they look ok
import matplotlib.pyplot as plt
fig = plt.figure(1)
plt.plot(xOut, yOut, '.')
plt.show()


# find the complete triangles that exist in cellsOnVertex after culling - we will only export those.

keepVerticesFullArray = np.zeros( (nVertices,) )
for v in range(nVertices):
  if (cellsOnVertex[v,0] in keepCells+1 and
      cellsOnVertex[v,1] in keepCells+1 and
      cellsOnVertex[v,2] in keepCells+1) :
    keepVerticesFullArray[v] = 1

keepVertices = np.nonzero( keepVerticesFullArray == 1 )[0]

# Now need to replace the indices in cellsOnVertex with the new condensed indices.
cellsOnVertexPrune = cellsOnVertex[keepVertices,:]

for v in range(len(keepVertices)):
    for c in range(3):
        oldindex = cellsOnVertexPrune[v,c]-1  # source info is 1-based, python is 0
        cellsOnVertexPrune[v,c] = cellIndicesFullToPruned[oldindex]

print cellsOnVertexPrune.max()
print cellIndicesFullToPruned.max()
 
print cellsOnVertex[keepVertices,:]
#print keepVertices
#print len(keepVertices), len(keepCells)

# ============================================
# ============================================
# FOR TRIANGLE 
# ============================================
# ============================================
# Build Triangle compatible node and ele files

# .node files
#    First line: <# of vertices> <dimension (must be 2)> <# of attributes> <# of boundary markers (0 or 1)>
#    Remaining lines: <vertex #> <x> <y> [attributes] [boundary marker] 
# Blank lines and comments prefixed by `#' may be placed anywhere. Vertices must be numbered consecutively, starting from one or zero.
# The attributes, which are typically floating-point values of physical quantities (such as mass or conductivity) associated with the nodes of a finite element mesh, are copied unchanged to the output mesh. If -q, -a, -u, or -s is selected, each new Steiner point added to the mesh will have quantities assigned to it by linear interpolation.

#.ele files
#    First line: <# of triangles> <nodes per triangle> <# of attributes>
#    Remaining lines: <triangle #> <node> <node> <node> ... [attributes] 
#Blank lines and comments prefixed by `#' may be placed anywhere. Triangles must be numbered consecutively, starting from one or zero. Nodes are indices into the corresponding .node file. The first three nodes are the corner vertices, and are listed in counterclockwise order around each triangle. (The remaining nodes, if any, depend on the type of finite element used.)
#As in .node files, the attributes are typically floating-point values of physical quantities (such as mass or conductivity) associated with the elements (triangles) of a finite element mesh. Because there is no simple mapping from input to output triangles, an attempt is made to interpolate attributes, which may result in a good deal of diffusion of attributes among nearby triangles as the triangulation is refined. Attributes do not diffuse across segments, so attributes used to identify segment-bounded regions remain intact.
#In output .ele files, all triangles have three nodes each

#np.savetxt('mesh.node', np.column_stack( (np.linspace(1, len(xOut), num=len(xOut), endpoint=True), xOut, yOut) ), delimiter=' ', comments='', header="%i 2 0 0"%len(xOut) , fmt='%i %.8f %.8f')

#np.savetxt('mesh.ele', np.column_stack( ( np.linspace(1, len(keepVertices), num=len(keepVertices), endpoint=True), cellsOnVertexPrune ) ),   delimiter=' ', fmt='%i %i %i %i', comments='', header="%i 3 0"%len(keepVertices) )




# ============================================
# ============================================
# FOR ascii_to_netcdf_packager.cpp
# ============================================
# ============================================

# TODO: ascii_to_netcdf_packager only works on spherical meshes right now?

#	This program reads three ascii files from the current directory:
#		* end_points.dat - This file should contain the x, y, and z coordinates
#						   for every cell center in the mesh. Each row is a
#						   point, and the columns are order x y z.

#		* triangles.dat - This file contains the indices for cells that make up
#						  each triangle. Each row is a triangle listing the
#						  indices for each cell that is a vertex of the
#						  triangle. Each column is an index of a triangle
#						  vertex.

#		* point_density.dat - This file contains the value of the density
#							  function evaluated at each cell center.


np.savetxt('end_points.dat', np.column_stack( (xOut, yOut, np.zeros((len(xOut),))) ), delimiter=' ')  #, header='x, y, z'  , fmt='%.8f'

np.savetxt('triangles.dat', cellsOnVertex[keepVertices,:], delimiter=' ', fmt='%i')  #, header='x, y, z')

np.savetxt('point_density.dat', meshDensity[keepCells], delimiter=' ')  #, header='x, y, z')


