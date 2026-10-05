#!/usr/bin/env python3

import argparse

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm
import numpy as np
import xarray as xr

import mosaic

SEC_PER_YEAR = 365.0 * 24.0 * 60.0 * 60.0
RHO_I = 910.0
RHO_W = 1028.0


def height_above_flotation(thickness, bed):
    """
    Compute height above flotation (HAF) given ice thickness and bed
    topography (both in meters, with bed relative to sea level).
    """

    flotation_thickness = np.maximum(0.0, -bed * (RHO_W / RHO_I))
    return np.maximum(thickness - flotation_thickness, 0.0)


def at_time(da, index=0):
    if "Time" in da.dims:
        return da.isel(Time=index)
    return da


def xtime_ymd(ds, index):
    """
    Return the 'YYYY-MM-DD' portion of xtime at the given time index.
    """

    xtime = ds["xtime"].isel(Time=index).values
    if hasattr(xtime, "tobytes"):
        xtime_str = xtime.tobytes().decode("utf-8")
    else:
        xtime_str = str(xtime)
    return xtime_str.strip().split("_")[0]


def boundary_segments(ds_mesh, mask):
    """
    Build line segments along cell-cell boundaries where mask changes
    from True to False.

    ``ds_mesh`` is expected to be a mosaic ``Descriptor.ds``-style dataset
    (zero-indexed connectivity arrays, with ``-1`` denoting a missing / mesh
    exterior neighbor). If ``ds_mesh`` has been culled down to a region via
    ``mosaic.utils.cull_mesh``, edges that were cut at the region boundary
    already have a ``-1`` neighbor, so they are naturally excluded below
    without needing a separate region mask.

    Uses verticesOnCell / cellsOnCell connectivity from the MPAS mesh.
    """

    mask = np.asarray(mask, dtype=bool)

    # cellsOnEdge gives the (at most) two cells adjacent to each edge, so
    # each shared cell-cell boundary is represented exactly once per edge
    # and this can be evaluated with array operations instead of a
    # per-cell/per-edge Python loop.
    cells_on_edge = ds_mesh["cellsOnEdge"].values.astype(int)
    vertices_on_edge = ds_mesh["verticesOnEdge"].values.astype(int)

    x_vertex = ds_mesh["xVertex"].values
    y_vertex = ds_mesh["yVertex"].values

    i_cell = cells_on_edge[:, 0]
    j_cell = cells_on_edge[:, 1]

    # Skip mesh exterior / invalid neighbors on either side of the edge.
    valid = (i_cell >= 0) & (j_cell >= 0)

    is_boundary = valid & (mask[i_cell.clip(min=0)] != mask[j_cell.clip(min=0)])

    v0 = vertices_on_edge[:, 0]
    v1 = vertices_on_edge[:, 1]
    is_boundary &= (v0 >= 0) & (v1 >= 0)

    v0 = v0[is_boundary]
    v1 = v1[is_boundary]

    segments = np.stack(
        [
            np.stack([x_vertex[v0], y_vertex[v0]], axis=-1),
            np.stack([x_vertex[v1], y_vertex[v1]], axis=-1),
        ],
        axis=1,
    ).tolist()
    return segments


def plot_geometry(ax, ds_mesh, thickness, bed, edge_color, gl_color, edge_ls="-", gl_ls="--"):
    thickness = np.asarray(thickness)
    ice_mask = thickness > 0.0
    flotation = np.asarray(bed) + (RHO_I / RHO_W) * thickness
    grounded_mask = (thickness > 0.0) & (flotation > 0.0)
    ice_segments = boundary_segments(ds_mesh, ice_mask)
    gl_segments = boundary_segments(ds_mesh, grounded_mask)

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
parser.add_argument(
    "--haf-range",
    type=float,
    default=None,
    help="One-sided colorbar range (max abs value, in m) for the height-above-flotation "
    "difference plot. Default: 90th percentile of abs(dhaf).",
)
parser.add_argument(
    "--speed-range",
    type=float,
    default=None,
    help="One-sided colorbar range (max abs value, in m/yr) for the surface-speed "
    "difference plot. Default: 90th percentile of abs(speed_diff).",
)
parser.add_argument(
    "--over-color",
    default="magenta",
    help="Color used for values above the colorbar range (default: %(default)s).",
)
parser.add_argument(
    "--under-color",
    default="purple",
    help="Color used for values below the colorbar range (default: %(default)s).",
)

args = parser.parse_args()

if (args.region_mask_file is None) != (args.region is None):
    parser.error("region_mask_file and region must either both be given " "or both be omitted")

ds_init = xr.open_dataset(args.initial_file)
ds1 = xr.open_dataset(args.output_file_1)
ds2 = xr.open_dataset(args.output_file_2)

h1 = at_time(ds1["thickness"], args.time_index_1)
h2 = at_time(ds2["thickness"], args.time_index_2)

date1 = xtime_ymd(ds1, args.time_index_1)
date2 = xtime_ymd(ds2, args.time_index_2)

# Bed always comes from the initial-condition file.
bed = at_time(ds_init["bedTopography"], 0)

# Build the Mosaic descriptor from the initial-condition file.
descriptor = mosaic.Descriptor(ds_init, use_latlon=False)

# -----------------------------------------------------------------------------
# Optional region subset
# -----------------------------------------------------------------------------

region_mask = None
index_to_cell_id = None

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

    if not np.any(region_mask):
        raise ValueError(f"Region '{region_label}' selects zero cells; nothing to plot.")

    x_region = ds_init.xCell.where(region_mask_da, drop=True)
    y_region = ds_init.yCell.where(region_mask_da, drop=True)

    xmin = float(x_region.min())
    xmax = float(x_region.max())
    ymin = float(y_region.min())
    ymax = float(y_region.max())

    # Physically cull the mesh down to just the selected region's cells
    # *before* any patches/boundary segments are built, so mosaic.polypcolor
    # and the ice-edge/grounding-line LineCollections only ever have to
    # process the (typically much smaller) region instead of the full mesh.
    # ``descriptor.ds`` is already the zero-indexed minimal mesh dataset that
    # mosaic.utils.cull_mesh expects (mirrors how Descriptor culls internally
    # for reprojection). Edges cut at the region boundary end up with a
    # cellsOnEdge value of -1, which boundary_segments already treats as a
    # mesh exterior, so no separate region-mask filtering is needed there.
    #
    # NOTE: ``descriptor.sizes`` (the *original*, un-culled mesh dimension
    # sizes) must be left untouched: mosaic.polypcolor/_get_array_location
    # uses it to detect that a data array is still full-mesh-sized and to
    # auto-subset it via the culled ds's ``indexToCellID`` lookup table. So
    # full-size data arrays (dhaf, speed_diff, etc.) should be passed to
    # mosaic.polypcolor unmodified; only arrays indexed directly against the
    # culled connectivity (i.e. in plot_geometry/boundary_segments) need to
    # be explicitly subset to the culled mesh via indexToCellID below.
    cells_to_cull = ~region_mask
    descriptor.ds = mosaic.utils.cull_mesh(descriptor.ds, cells_to_cull)
    index_to_cell_id = descriptor.ds["indexToCellID"].values

region_suffix = f"_region{args.region}" if args.region is not None else ""

dx = xmax - xmin
dy = ymax - ymin
pad_x = 0.03 * dx if dx > 0.0 else 1000.0
pad_y = 0.03 * dy if dy > 0.0 else 1000.0


# =============================================================================
# Height-above-flotation difference
# =============================================================================

haf1 = height_above_flotation(h1, bed)
haf2 = height_above_flotation(h2, bed)
dhaf = haf2 - haf1
if args.haf_range is not None:
    max_abs_dhaf = args.haf_range
else:
    max_abs_dhaf = float(np.nanpercentile(np.abs(dhaf.values), 90))
if max_abs_dhaf == 0.0:
    max_abs_dhaf = 1.0
fig, ax = plt.subplots(figsize=(10, 8), constrained_layout=True)
cmap_haf = plt.get_cmap("RdBu_r").copy()
cmap_haf.set_over(args.over_color)
cmap_haf.set_under(args.under_color)
pc = mosaic.polypcolor(ax, descriptor, dhaf, cmap=cmap_haf, vmin=-max_abs_dhaf, vmax=max_abs_dhaf, edgecolors="none")

bdy1_color = 'b'
gl1_color = 'g'
bdy2_color = 'c'
gl2_color = 'lime'
# plot_geometry indexes directly into the (possibly culled) mesh
# connectivity, so thickness/bed must match descriptor.ds's current size.
h1_mesh = h1.isel(nCells=index_to_cell_id).values if index_to_cell_id is not None else h1.values
h2_mesh = h2.isel(nCells=index_to_cell_id).values if index_to_cell_id is not None else h2.values
bed_mesh = bed.isel(nCells=index_to_cell_id).values if index_to_cell_id is not None else bed.values
plot_geometry(ax, descriptor.ds, h1_mesh, bed_mesh, edge_color=bdy1_color, gl_color=gl1_color, edge_ls="-", gl_ls="-")
plot_geometry(ax, descriptor.ds, h2_mesh, bed_mesh, edge_color=bdy2_color, gl_color=gl2_color, edge_ls="-", gl_ls="-")
ax.set_xlim(xmin - pad_x, xmax + pad_x)
ax.set_ylim(ymin - pad_y, ymax + pad_y)
ax.set_aspect("equal")
ax.set_xlabel("x [m]")
ax.set_ylabel("y [m]")
ax.set_title(
    f"Height above flotation change: "
    f"{date2} - "
    f"{date1}\n"
    f"Region: {region_label}"
)

fig.colorbar(pc, ax=ax, label="Height above flotation difference [m]", extend="both")
ax.plot([], [], color=bdy1_color, ls="-", label="Ice edge, time 1")
ax.plot([], [], color=gl1_color, ls="-", label="Grounding line, time 1")
ax.plot([], [], color=bdy2_color, ls="-", label="Ice edge, time 2")
ax.plot([], [], color=gl2_color, ls="-", label="Grounding line, time 2")
ax.legend(loc="best")

haf_filename = f"haf_difference{region_suffix}.png"
fig.savefig(haf_filename, dpi=300)

#plt.close(fig)
print(f"Wrote {haf_filename}")

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
if args.speed_range is not None:
    max_abs_du = args.speed_range
else:
    max_abs_du = float(np.nanpercentile(np.abs(speed_diff.values), 90))
if max_abs_du == 0.0:
    max_abs_du = 1.0
fig2, ax2 = plt.subplots(figsize=(10, 8), constrained_layout=True)

bdy1_color = 'b'
gl1_color = 'g'
cmap_speed = plt.get_cmap("RdBu_r").copy()
cmap_speed.set_over(args.over_color)
cmap_speed.set_under(args.under_color)
pc = mosaic.polypcolor(ax2, descriptor, speed_diff, cmap=cmap_speed, vmin=-max_abs_du, vmax=max_abs_du, edgecolors="none")
plot_geometry(ax2, descriptor.ds, h2_mesh, bed_mesh, edge_color=bdy1_color, gl_color=gl1_color, edge_ls="-", gl_ls="-")
ax2.set_xlim(xmin - pad_x, xmax + pad_x)
ax2.set_ylim(ymin - pad_y, ymax + pad_y)
ax2.set_aspect("equal")
ax2.set_xlabel("x [m]")
ax2.set_ylabel("y [m]")
ax2.set_title(f"Modeled - observed surface speed: " f"{date2}\n" f"Region: {region_label}")
fig2.colorbar(pc, ax=ax2, label=("Surface speed difference " "[m yr$^{-1}$]"), extend="both")
ax2.plot([], [], color=bdy1_color, ls="-", label="Ice edge")
ax2.plot([], [], color=gl1_color, ls="-", label="Grounding line")
ax2.legend(loc="best")
speed_filename = f"surface_speed_difference{region_suffix}.png"
fig2.savefig(speed_filename, dpi=300)
#plt.close(fig2)
print(f"Wrote {speed_filename}")

# =============================================================================
# Modeled vs. observed surface speed, 1:1 heatmap, split by grounded/floating
# =============================================================================
axis_lo, axis_hi = -2.0, 5.0

flotation_full = bed.values + (RHO_I / RHO_W) * h2.values
grounded_mask_full = (h2.values > 0.0) & (flotation_full > 0.0)
floating_mask_full = (h2.values > 0.0) & (flotation_full <= 0.0)

obs_speed_vals = obs_speed.values
model_speed_vals = model_speed.values
valid = (
    np.isfinite(obs_speed_vals)
    & np.isfinite(model_speed_vals)
    & (obs_speed_vals > 0.0)
    & (model_speed_vals > 0.0)
)
# Restrict the heatmap to the selected region, if one was given, so it
# matches the subset shown in the map plots above.
if region_mask is not None:
    valid = valid & region_mask

heatmap_bins = np.linspace(axis_lo, axis_hi, 141)
fig3, (ax3_grounded, ax3_floating) = plt.subplots(1, 2, figsize=(16, 8), constrained_layout=True)

for ax3, ice_mask, panel_label in (
    (ax3_grounded, grounded_mask_full, "Grounded ice"),
    (ax3_floating, floating_mask_full, "Floating ice"),
):
    panel_valid = valid & ice_mask
    log_obs_speed = np.log10(obs_speed_vals[panel_valid])
    log_model_speed = np.log10(model_speed_vals[panel_valid])

    _, _, _, heatmap_img = ax3.hist2d(
        log_obs_speed, log_model_speed, bins=heatmap_bins, cmap="viridis", norm=LogNorm()
    )
    fig3.colorbar(heatmap_img, ax=ax3, label="Count")
    ax3.plot([axis_lo, axis_hi], [axis_lo, axis_hi], color="k", linewidth=0.8, linestyle="--", label="1:1")
    ax3.set_xlim(axis_lo, axis_hi)
    ax3.set_ylim(axis_lo, axis_hi)
    ax3.set_aspect("equal")
    ax3.set_xlabel("log10(observed speed) [log10(m yr$^{-1}$)]")
    ax3.set_ylabel("log10(modeled speed) [log10(m yr$^{-1}$)]")
    ax3.set_title(panel_label)
    ax3.legend(loc="best")

fig3.suptitle(f"Modeled vs. observed surface speed: " f"{date2}\n" f"Region: {region_label}")
hist_filename = f"surface_speed_heatmap{region_suffix}.png"
fig3.savefig(hist_filename, dpi=300)
#plt.close(fig3)
print(f"Wrote {hist_filename}")

#plt.show()
