#!/usr/bin/env python3

# TODO: 
# adjust normalization for region area


import argparse
import numpy as np
import os
import subprocess
import sys
import xarray as xr
from pyremap.descriptor.lat_lon_grid_descriptor import LatLonGridDescriptor
from pyremap.descriptor.mpas_mesh_descriptor import MpasMeshDescriptor
from pyremap.remapper import Remapper

def create_scrip_from_latlon(input_file, lat_var='latitude', lon_var='longitude'):
    fname = os.path.basename(input_file)
    fprefix = os.path.splitext(fname)[0]
    grid = LatLonGridDescriptor.read(input_file, latVarName=lat_var, lonVarName=lon_var)
    # create the scrip file in the run dir
    scrip_filename = f"{fprefix}_scrip.nc"
    grid.to_scrip(scrip_filename)
    return grid, scrip_filename

def create_scrip_from_mpas(mesh_file):
    fname = os.path.basename(mesh_file)
    fprefix = os.path.splitext(fname)[0]
    mesh = MpasMeshDescriptor(mesh_file, meshName=fprefix)
    # create the scrip file in the run dir
    scrip_filename = f"{fprefix}_scrip.nc"
    mesh.to_scrip(scrip_filename)
    print(f"Saved scrip file {scrip_filename}")
    return mesh, scrip_filename

def create_mapfile_base(src_scrip, dest_scrip):
    mapfile_base = "map_src_to_dest.nc"
    # delete file if it exists because Remapper won't clobber
    if os.path.isfile(mapfile_base):
        os.remove(mapfile_base)
    remapper = Remapper(
        sourceDescriptor=src_scrip,
        destinationDescriptor=dest_scrip,
        mappingFileName=mapfile_base,
    )
    remapper.esmf_build_map(method='bilinear')  # method doesn't matter
    ds = xr.open_dataset(mapfile_base)
    # drop the actual mapping, keeping all the other fields
    ds = ds.drop_vars(['col', 'row', 'S'])
    print(ds)
    return ds
 
def remap_files(climatology_file, ocn_mesh_file, glc_mesh_file):
    ds = xr.open_dataset(climatology_file)
    regions = ds.sizes['region']

    # Create scrip files
    source_scrip, src_scrip_file  = create_scrip_from_latlon(climatology_file)
    glc_scrip, glc_scrip_file = create_scrip_from_mpas(glc_mesh_file)
    ocn_scrip, ocn_scrip_file = create_scrip_from_mpas(ocn_mesh_file)

    # add ice shelf mask to ocn scrip
    ds_ocn = xr.open_dataset(ocn_mesh_file)
    imask = (ds_ocn['landIceMask'][0,:].values == 0).astype('i')
    ds_ocn_scrip = xr.open_dataset(ocn_scrip_file)
    ds_ocn_scrip['grid_imask'].data = imask
    # Define encoding to suppress _FillValue for all variables
    encoding = {var: {"_FillValue": None} for var in ds_ocn_scrip.data_vars}
    ds_ocn_scrip.to_netcdf(ocn_scrip_file+'2', encoding=encoding, format="NETCDF4")
    ds_ocn_scrip.close()

    print(f"Creating melt->ocn mapping")
    melt_to_ocn_mapping_file = "map_meltgrid_to_ocn_bilinear.nc"
    # delete file if it exists because Remapper won't clobber
    if os.path.isfile(melt_to_ocn_mapping_file):
        os.remove(melt_to_ocn_mapping_file)

    # Build mapping file.  Can't use pyremap because it does 
    # not support masks in a scrip file as of now
    cargs = ['ESMF_RegridWeightGen',
             '--source', src_scrip_file,
             '--destination', ocn_scrip_file,
             '--dst_regional',
             '--weight', melt_to_ocn_mapping_file,
             '--method', 'bilinear',
             '--netcdf4']
    subprocess.check_call(cargs)

    # Create melt->ocn map files for each region
    if not args.skip_melt_to_ocn:
        melt_frac = np.zeros(regions)
        for region_idx in range(regions):
            # melt: The spatial integral on the spherical Earth summed over the 12 months and all regions is equal to 1.0
            region_melt = ds['melt'].isel(region=region_idx).sum(dim='time')
            region_ds = region_melt.to_dataset(name='melt')
            region_melt_file = f"region_{region_idx}_annual_mean_melt.nc"
            region_ds.to_netcdf(region_melt_file)
            region_ds.close()

            # perform remapping
            region_melt_file_ocn = f"region_{region_idx}_annual_mean_melt_ocn.nc"
            cargs = ['ncremap',
                     '-m', melt_to_ocn_mapping_file,
                     '-P', 'mpas',
                     '-v', 'melt',
                     region_melt_file, region_melt_file_ocn]
            print(cargs)
            subprocess.check_call(cargs)

            # renormalize
            ds_melt_ocn = xr.open_dataset(region_melt_file_ocn)
            melt = ds_melt_ocn['melt'].values
            melt_sum = np.nansum(melt * ds_ocn['areaCell'].values)
            melt_frac[region_idx] = melt_sum
            print(f'sum of melt on ocn mesh={melt_sum}  Renormalizing.')
            melt /= melt_sum
            print(f'sum of melt on ocn mesh={np.nansum(melt * ds_ocn["areaCell"].values)}')
            ds_melt_ocn['melt'].data = melt
            ds_melt_ocn.to_netcdf(region_melt_file_ocn, format="NETCDF3_CLASSIC")
    ds_ocn.close()
    print(f'Sum of melt acroos all regions on MPASO mesh before normalization (should be close but not exactly equal to 1.0): {melt_frac.sum()}')

    # Create rgn->glc map file
    print("\nCreate mapping file and remap region->glc")
    rgn_to_glc_mapping_file = "map_regions_to_glc_nstod.nc"
    # delete file if it exists because Remapper won't clobber
    if os.path.isfile(rgn_to_glc_mapping_file):
        os.remove(rgn_to_glc_mapping_file)
    remapper = Remapper(
        sourceDescriptor=source_scrip,
        destinationDescriptor=glc_scrip,
        mappingFileName=rgn_to_glc_mapping_file,
    )
    remapper.esmf_build_map(method='neareststod')
    region_file_glc = f"regions_glc.nc"
    remapper.remap_file(climatology_file, region_file_glc,
                        variableList=['region_map_expanded'], overwrite=True)

    # Build custom mapping file
    # Create generic mapping file using ESMF to get the metadata fields
    print("\nBuilding custom mapping file")
    ds_custom_map = create_mapfile_base(glc_scrip, ocn_scrip)
    ds_custom_map.to_netcdf('tmp.nc')
    ds_rgn_glc = xr.open_dataset(region_file_glc)
    first_time = True
    for region_idx in range(regions):
        # get glc cells in this region
        glc_rgn_idx = np.nonzero(ds_rgn_glc.region_map_expanded.values == region_idx + 1)[0]
        n_glc_idx = len(glc_rgn_idx)
        # get ocn cells for this region
        ds_ocn = xr.open_dataset(f"region_{region_idx}_annual_mean_melt_ocn.nc")
        melt = ds_ocn.melt.values
        ocn_rgn_idx = np.nonzero(melt > 0)[0]
        n_ocn_idx = len(ocn_rgn_idx)
        print(f'Region={region_idx + 1}, n_glc_idx={n_glc_idx}, n_ocn_idx={n_ocn_idx}')

        # now append the relevant entries for this region
        new_col = np.repeat(glc_rgn_idx + 1, n_ocn_idx)
        new_row = np.tile(ocn_rgn_idx + 1, n_glc_idx)
        new_S   = np.tile(melt[ocn_rgn_idx], n_glc_idx)
        if first_time:
            col = new_col
            row = new_row
            S = new_S
            first_time = False
        else:
            col = np.concatenate([col, new_col])
            row = np.concatenate([row, new_row])
            S   = np.concatenate([S,   new_S])
    # save to mapping file
    print(f"n_s={len(S)}")
    ds_custom_map['col'] = (['n_s',], col)
    ds_custom_map['row'] = (['n_s',], row)
    ds_custom_map['S']   = (['n_s',], S)
    # Save to NetCDF (ESMF prefers classic format)
    encoding = {
        "row": {"dtype": "int32"},
        "col": {"dtype": "int32"},
        "S": {"dtype": "float64"},
    }
    filename='final_map.nc'
    ds_custom_map.to_netcdf(filename, encoding=encoding)
    print(f"Created dummy ESMF mapping file: {filename}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create remapping files for iceberg melt data.")
    parser.add_argument("--ocn_mesh", required=True, help="Path to ocean mesh file.")
    parser.add_argument("--glc_mesh", required=True, help="Path to glacier (ice-sheet) mesh file.")
    parser.add_argument("--melt_clim", required=True, help="Path to iceberg melt climatology NetCDF file. i.e. AQ_iceberg_melt.nc")
    parser.add_argument("--skip_melt_to_ocn", help="if the melt to ocn step should be skipped (assumes has already been run)", action='store_true')
    args = parser.parse_args()

    remap_files(args.melt_clim, args.ocn_mesh, args.glc_mesh)
