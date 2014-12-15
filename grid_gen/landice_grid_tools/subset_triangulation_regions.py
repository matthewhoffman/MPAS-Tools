# Script to take CISM output and convert to ascii files of the format:
# latitude, longitude, thickness
# This involves reprojecting from polar stereographic to lat long
# and then building 1-d arrays to save as text files for each time slice.

import netCDF4
import numpy as np



####################### SET THIS STUFF ###################

infile = 'grid.100.nc.mpas'

saveLatMax = -40.0

#saveLatMax = -45.0
#saveLatMin = -91.0

keepCriterion = 1  


useAnyVertex = True   # True if any vertex of the cell should be below the threshold.  False if just the cell center should be.

keepOnlyFullTriangles = True # True if triangles are kept only if all cells in the triangles are kept.  False if triangles are kept if ANY cells in the triangles are kept.

##########################################




d2r = np.pi / 180.0
r2d = 180.0 / np.pi


# Get stuff from main file
f = netCDF4.Dataset(infile,'r')
lat = f.variables['latCell'][:]
lon = f.variables['lonCell'][:]
latv = f.variables['latVertex'][:]
lonv = f.variables['lonVertex'][:]
cellsOnVertex = f.variables['cellsOnVertex'][:]
verticesOnCell = f.variables['verticesOnCell'][:]
nVertices = len(f.dimensions['nVertices'])
nCells = len(f.dimensions['nCells'])
maxEdges = len(f.dimensions['maxEdges'])
meshDensity = f.variables['meshDensity'][:]



# 1. find all cells that have any vertex below the lat threshold
keepCellFullList = np.zeros( (nCells,) )
for c in range(nCells):
  if useAnyVertex==True:
    ind = np.nonzero(verticesOnCell[c,:]>0)
    minlat = latv[verticesOnCell[c,ind]-1].min() * r2d
    if minlat < saveLatMax:
      keepCellFullList[c] = 1
  else:
    if lat[c] * r2d < saveLatMax:
      keepCellFullList[c] = 1
keepCells = np.nonzero( keepCellFullList )[0]
print 'Found %i of %i original cells below latitude threshold'%(len(keepCells), nCells)

# 2. find any triangles that have any of the cells we are keeping
# For some unknown reason, RegionTriangulation has to include complete triangles.
# If this restriction is understood and relaxed, steps 2 and 3 would not be necessary.
keepVerticesFullArray = np.zeros( (nVertices,) )
if keepOnlyFullTriangles == False:
  for v in range(nVertices):
    if (cellsOnVertex[v,0] in keepCells+1 or
        cellsOnVertex[v,1] in keepCells+1 or
        cellsOnVertex[v,2] in keepCells+1) :
      keepVerticesFullArray[v] = 1
else:
  for v in range(nVertices):
    if (cellsOnVertex[v,0] in keepCells+1 and
        cellsOnVertex[v,1] in keepCells+1 and
        cellsOnVertex[v,2] in keepCells+1) :
      keepVerticesFullArray[v] = 1
keepVertices = np.nonzero( keepVerticesFullArray == 1 )[0]
print 'Keeping %i of %i original triangles'%(len(keepVertices), nVertices)

# 3. Now rebuild the keepCells array to include all cells in all triangles that we have now identified
# (if keepOnlyFullTriangles==True, this shouldn't do anything.)
keepCellFullList = np.zeros( (nCells,) )
for v in keepVertices:
  for c in range(3):
    # Note that each cell will be set to 1 multiple times - that's ok here
    keepCellFullList[cellsOnVertex[v,c]-1] = 1  # cellsOnVertex is 1-based, but we are indexing a 0-based array
keepCells = np.nonzero( keepCellFullList )[0]
print 'Needed %i of %i original cells to fill out triangulation'%(len(keepCells), nCells)





lonOut = lon[keepCells]
latOut = lat[keepCells]

# --------
# Now rebuild the triangulation with indices from 1 to len(keepCells)
# --------

# Build lookup table from full indices to new pruned indices
cellIndicesFullToPruned = np.zeros( ( nCells, ) )
counter = 0 # make the new indices 1-based
for c in keepCells:
    counter = counter + 1
    cellIndicesFullToPruned[c] = counter

#print cellIndicesFullToPruned

##Optional - plot out the lat long positions to make sure they look ok
import matplotlib.pyplot as plt
fig = plt.figure(1)
plt.plot(lon*r2d, lat*r2d, '.')
plt.plot(lonOut*r2d, latOut*r2d, 'o')
for v in keepVertices:
    plt.plot( np.hstack( (lon[cellsOnVertex[v,:]-1], lon[cellsOnVertex[v,0]-1]) )*r2d, np.hstack( (lat[cellsOnVertex[v,:]-1], lat[cellsOnVertex[v,0]-1]) )*r2d, '-k')
#plt.plot(xv, yv, 'xr')
plt.plot((0, 360.0), (saveLatMax, saveLatMax), ':r')
plt.show()


# Now need to replace the indices in cellsOnVertex with the new condensed indices.
cellsOnVertexPrune = cellsOnVertex[keepVertices,:]
for v in range(len(keepVertices)):
    for c in range(3):
        oldindex = cellsOnVertexPrune[v,c]-1  # source info is 1-based, python is 0
        if oldindex in keepCells:
          cellsOnVertexPrune[v,c] = cellIndicesFullToPruned[oldindex]
        else:
          cellsOnVertexPrune[v,c] = -1  # mark this as a degenerate cell location

#print cellsOnVertexPrune.max()
#print cellIndicesFullToPruned.max()
 
#print cellsOnVertex[keepVertices,:]
#print keepVertices
#print len(keepVertices), len(keepCells)

# Now deal with degenerate triangles
for v in range(len(keepVertices)):
    for c in range(3):
        if cellsOnVertexPrune[v,c] == -1:
            print 'ERROR: degenerate triangle found!  There should not be any anymore!'
#            cellsOnVertexPrune[v,c] = cellsOnVertexPrune[v,:].max()  # assign the maximum cell value on this vertex here to get somthing valid


# convert cell lat/lon to unit circle x,y,z
xOut = np.cos(lonOut) * np.cos(latOut);
yOut = np.sin(lonOut) * np.cos(latOut);
zOut = np.sin(latOut);

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


np.savetxt('end_points.dat', np.column_stack( (xOut, yOut, zOut) ), delimiter=' ')  #, header='x, y, z'  , fmt='%.8f'

#np.savetxt('triangles.dat', cellsOnVertex[keepVertices,:]-1, delimiter=' ', fmt='%i')  #, header='x, y, z')

np.savetxt('triangles.dat', cellsOnVertexPrune-1, delimiter=' ', fmt='%i')  #, header='x, y, z')

np.savetxt('point_density.dat', meshDensity[keepCells], delimiter=' ')  #, header='x, y, z')


