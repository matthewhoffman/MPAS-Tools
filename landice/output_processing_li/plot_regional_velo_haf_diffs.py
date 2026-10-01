#!/usr/bin/env python3

import argparse

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np
import xarray as xr

import mosaic

SEC_PER_YEAR = 365.0 * 24.0 * 60.0 * 60.0
RHO_I = 910.0
RHO_W = 1028.0


def at_time(da, index=0):
    if "Time" in da.dims:
        return da.isel(Time=index)
    return da


def boundary_segments(ds_mesh, mask, region_mask=None):
    """
    Build line segments along cell-cell boundaries where mask changes
    from True to False.

    Uses verticesOnCell / cellsOnCell connectivity from the MPAS mesh.
    """

    mask = np.asarray(mask, dtype=bool)

    if region_mask is not None:
        region_mask = np.asarray(region_mask, dtype=bool)

    n_cells = ds_mesh.sizes["nCells"]

    n_edges_on_cell = ds_mesh["nEdgesOnCell"].values
    cells_on_cell = ds_mesh["cellsOnCell"].values.astype(int) - 1
    edges_on_cell = ds_mesh["edgesOnCell"].values.astype(int) - 1
    vertices_on_edge = ds_mesh["verticesOnEdge"].values.astype(int) - 1

    x_vertex = ds_mesh["xVertex"].values
    y_vertex = ds_mesh["yVertex"].values

    segments = []

    for i_cell in range(n_cells):
        if region_mask is not None and not region_mask[i_cell]:
            continue
        n_edges = n_edges_on_cell[i_cell]
        for i_edge in range(n_edges):
            j_cell = cells_on_cell[i_cell, i_edge]
            # Skip mesh exterior / invalid neighbor.
            if j_cell < 0:
                continue
            # Only process each shared edge once.
            if j_cell < i_cell:
                continue
            if region_mask is not None:
                if not region_mask[j_cell]:
                    continue
            if mask[i_cell] == mask[j_cell]:
                continue
            # Edge i connects vertex i to vertex i+1 around the cell.
            found_edge = edges_on_cell[i_cell, i_edge]
            v0 = vertices_on_edge[found_edge, 0]
            v1 = vertices_on_edge[found_edge, 1]
            if v0 < 0 or v1 < 0:
                continue
            segments.append([(x_vertex[v0], y_vertex[v0]), (x_vertex[v1], y_vertex[v1])])
    return segments


def plot_geometry(ax, ds_mesh, thickness, bed, region_mask, edge_color, gl_color,  edge_ls="-", gl_ls="--"):
    thickness = np.asarray(thickness)
    ice_mask = thickness > 0.0
    flotation = np.asarray(bed) + (RHO_I / RHO_W) * thickness
    grounded_mask = (thickness > 0.0) & (flotation > 0.0)
    ice_segments = boundary_segments(ds_mesh, ice_mask, region_mask=region_mask)
    gl_segments = boundary_segments(ds_mesh, grounded_mask, region_mask=region_mask)

    if len(ice_segments) > 0:
        ax.add_collection(LineCollection(ice_segments, colors=edge_color, linewidths=0.5, linestyles=edge_ls))

    if len(gl_segments) > 0:
        ax.add_collection(LineCollection(gl_segments, colors=gl_color, linewidths=0.5, linestyles=gl_ls))


parser = argparse.ArgumentParser(description="Plot MALI 2-D thickness change and surface-speed misfit")

parser.add_argument("initial_file")
parser.add_argument("output_file_1")
parser.add_argument("time_index_1", type=int)
parser.add_argument("output_file_2")
parser.add_argument("time_index_2", type=int)
parser.add_argument("region_mask_file", nargs="?", default=None)
parser.add_argument("region", nargs="?", default=None)

args = parser.parse_args()

if (args.region_mask_file is None) != (args.region is None):
    parser.error("region_mask_file and region must either both be given " "or both be omitted")

ds_init = xr.open_dataset(args.initial_file)
ds1 = xr.open_dataset(args.output_file_1)
ds2 = xr.open_dataset(args.output_file_2)

h1 = at_time(ds1["thickness"], args.time_index_1)
h2 = at_time(ds2["thickness"], args.time_index_2)

# Bed always comes from the initial-condition file.
bed = at_time(ds_init["bedTopography"], 0)

# Build the Mosaic descriptor from the initial-condition file.
descriptor = mosaic.Descriptor(ds_init, use_latlon=False)

# -----------------------------------------------------------------------------
# Optional region subset
# -----------------------------------------------------------------------------

region_mask = None

xmin = float(ds_init.xCell.min())
xmax = float(ds_init.xCell.max())
ymin = float(ds_init.yCell.min())
ymax = float(ds_init.yCell.max())

region_label = "full mesh"


if args.region_mask_file is not None:

    ds_regions = xr.open_dataset(args.region_mask_file)
    # Allow either zero-based region index or region name.
    try:
        region_index = int(args.region)
    except ValueError:
        names = None
        for name_var in ("regionNames", "regionMaskNames"):
            if name_var not in ds_regions:
                continue
            raw = ds_regions[name_var].values
            if raw.ndim == 1:
                names = [str(value.decode() if isinstance(value, bytes) else value).strip() for value in raw]
            else:
                names = []
                for row in raw:
                    chars = []
                    for char in row:
                        if isinstance(char, bytes):
                            chars.append(char.decode("utf-8"))
                        else:
                            chars.append(str(char))
                    names.append("".join(chars).replace("\x00", "").strip())
            break
        if names is None:
            raise ValueError(
                "Region was specified by name, but " "the region file has neither " "regionNames nor regionMaskNames"
            )
        if args.region not in names:
            raise ValueError(f"Region '{args.region}' not found.\n" "Available regions are:\n" + "\n".join(names))
        region_index = names.index(args.region)

    region_mask_da = ds_regions["regionCellMasks"].isel(nRegions=region_index).astype(bool)
    region_mask = region_mask_da.values
    region_label = args.region
    x_region = ds_init.xCell.where(region_mask_da, drop=True)
    y_region = ds_init.yCell.where(region_mask_da, drop=True)

    xmin = float(x_region.min())
    xmax = float(x_region.max())
    ymin = float(y_region.min())
    ymax = float(y_region.max())

dx = xmax - xmin
dy = ymax - ymin
pad_x = 0.03 * dx if dx > 0.0 else 1000.0
pad_y = 0.03 * dy if dy > 0.0 else 1000.0


# =============================================================================
# Thickness difference
# =============================================================================

dh = h2 - h1
if region_mask is not None:
    dh = dh.where(xr.DataArray(region_mask, dims=("nCells",)))
max_abs_dh = float(np.nanmax(np.abs(dh.values)))
if max_abs_dh == 0.0:
    max_abs_dh = 1.0
fig, ax = plt.subplots(figsize=(10, 8), constrained_layout=True)
max_abs_dh = 50
pc = mosaic.polypcolor(ax, descriptor, dh, cmap="RdBu_r", vmin=-max_abs_dh, vmax=max_abs_dh, edgecolors="none")

bdy1_color = 'b'
gl1_color = 'g'
bdy2_color = 'c'
gl2_color = 'lime'
plot_geometry(ax, ds_init, h1.values, bed.values, region_mask, edge_color=bdy1_color, gl_color=gl1_color, edge_ls="-", gl_ls="-")
plot_geometry(ax, ds_init, h2.values, bed.values, region_mask, edge_color=bdy2_color, gl_color=gl2_color, edge_ls="-", gl_ls="-")
ax.set_xlim(xmin - pad_x, xmax + pad_x)
ax.set_ylim(ymin - pad_y, ymax + pad_y)
ax.set_aspect("equal")
ax.set_xlabel("x [m]")
ax.set_ylabel("y [m]")
ax.set_title(
    f"Thickness change: "
    f"{args.output_file_2}"
    f"[{args.time_index_2}] - "
    f"{args.output_file_1}"
    f"[{args.time_index_1}]\n"
    f"Region: {region_label}"
)

fig.colorbar(pc, ax=ax, label="Thickness difference [m]")
ax.plot([], [], color=bdy1_color, ls="-", label="Ice edge, time 1")
ax.plot([], [], color=gl1_color, ls="-", label="Grounding line, time 1")
ax.plot([], [], color=bdy2_color, ls="-", label="Ice edge, time 2")
ax.plot([], [], color=gl2_color, ls="-", label="Grounding line, time 2")
ax.legend(loc="best")

fig.savefig("thickness_difference.png", dpi=300)

#plt.close(fig)
print("Wrote thickness_difference.png")

# =============================================================================
# Surface-speed difference from observations
#
# Uses the SECOND requested model time slice.
# =============================================================================
obs_u = at_time(ds_init["observedSurfaceVelocityX"], 0)
obs_v = at_time(ds_init["observedSurfaceVelocityY"], 0)
obs_speed = np.sqrt(obs_u**2 + obs_v**2) * SEC_PER_YEAR
model_speed = at_time(ds2["surfaceSpeed"], args.time_index_2) * SEC_PER_YEAR
speed_diff = model_speed - obs_speed
# Don't plot velocity error
# where the model has no ice.
speed_diff = speed_diff.where(h2 > 0.0)
if region_mask is not None:
    speed_diff = speed_diff.where(xr.DataArray(region_mask, dims=("nCells",)))
max_abs_du = float(np.nanmax(np.abs(speed_diff.values)))
if max_abs_du == 0.0:
    max_abs_du = 1.0
fig2, ax2 = plt.subplots(figsize=(10, 8), constrained_layout=True)

bdy1_color = 'b'
gl1_color = 'g'
max_abs_du = 400
pc = mosaic.polypcolor(ax2, descriptor, speed_diff, cmap="RdBu_r", vmin=-max_abs_du, vmax=max_abs_du, edgecolors="none")
plot_geometry(ax2, ds_init, h2.values, bed.values, region_mask, edge_color=bdy1_color, gl_color=gl1_color, edge_ls="-", gl_ls="-")
ax2.set_xlim(xmin - pad_x, xmax + pad_x)
ax2.set_ylim(ymin - pad_y, ymax + pad_y)
ax2.set_aspect("equal")
ax2.set_xlabel("x [m]")
ax2.set_ylabel("y [m]")
ax2.set_title(f"Modeled - observed surface speed: " f"{args.output_file_2}" f"[{args.time_index_2}]\n" f"Region: {region_label}")
fig2.colorbar(pc, ax=ax2, label=("Surface speed difference " "[m yr$^{-1}$]"))
ax2.plot([], [], color=bdy1_color, ls="-", label="Ice edge")
ax2.plot([], [], color=gl1_color, ls="-", label="Grounding line")
ax2.legend(loc="best")
fig2.savefig("surface_speed_difference.png", dpi=300)
#plt.close(fig2)
print("Wrote surface_speed_difference.png")

#plt.show()
