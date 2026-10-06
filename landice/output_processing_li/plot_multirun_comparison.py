#!/usr/bin/env python3
"""
plot_multirun_comparison.py

Compare many MALI Antarctic Ice Sheet runs (each a subdirectory of the
directory this script is run from) across regional grounded-mass-balance /
volume-above-flotation (VAF) time series and regional surface-velocity
misfit (grounded/floating speed heatmaps + a modeled-minus-observed speed
map), producing one comparison figure per ISMIP6 region plus one
whole-ice-sheet summary figure.

This is a prototype script assembled directly from logic/data in two
existing scripts in this directory (intentionally *not* generalized or
refactored into shared modules):
  * plot_regionalStats2.py           - grounded-MB / VAF time series plotting
                                        and the Rignot et al. (2019) regional
                                        net-mass-balance obs dataset.
  * plot_regional_velo_haf_diffs.py  - mosaic-based mesh culling/plotting,
                                        grounded/floating speed heatmaps, and
                                        the modeled-vs-observed speed-
                                        difference map.

Assumed directory structure (run this script from the parent directory
containing the run subdirectories):

    <rundir>/streams.landice             (XML; stream "input" gives the MALI
                                           mesh/IC file; stream "regionsInput"
                                           gives the region-mask file)
    <rundir>/output/regionalStats.nc
    <rundir>/output/output_2d_2005.nc

All runs are assumed to share the identical mesh and region-mask file, so
the mesh/IC file and region-mask file are resolved only once -- from the
first (alphabetically sorted) run subdirectory -- and reused for every run.

At most 8 run subdirectories are supported (one row per run in each
per-region figure, plus a summary row); if more are found the script exits
with an error asking the user to reduce the number of run directories.

Usage:
    python plot_multirun_comparison.py [--dpi DPI]
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import argparse
import copy
import os
import re
import sys
import xml.etree.ElementTree as ET

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import xarray as xr

import mosaic

plt.rcParams["text.hinting"] = "no_hinting"

RHO_I = 910.0
RHO_W = 1028.0
SEC_PER_YEAR = 365.0 * 24.0 * 60.0 * 60.0
GT_PER_M3_ICE = RHO_I / 1.0e12  # volume (m3) -> mass (Gt), using ice density

MAX_RUNS = 8
OUTPUT_DIR = "multirun_comparison_figs"
OUTPUT_2D_RELPATH = os.path.join("output", "output_2d_2005.nc")
REGIONAL_STATS_RELPATH = os.path.join("output", "regionalStats.nc")

HEATMAP_AXIS_LO, HEATMAP_AXIS_HI = -2.0, 5.0
HEATMAP_BINS = np.linspace(HEATMAP_AXIS_LO, HEATMAP_AXIS_HI, 141)

SPEED_DIFF_MAX_ABS = 50.0  # hard-coded one-sided colorbar range [m/yr] for speed-error maps

# ============================================================================
# Rignot et al. (2019) regional net mass-balance obs (PNAS 116(4), 1095-1103;
# ice discharge D09-17 for 2009-2017, net mass balance uses Rignot 2008 SMB),
# copied directly from the OBSERVATIONAL_DATASETS registry in
# plot_regionalStats2.py. [mean, 1-sigma uncertainty] in Gt/yr.
# ============================================================================
RIGNOT_2019_NET_MB = {
    'ISMIP6BasinAAp': [-1.6, 1.6],
    'ISMIP6BasinApB': [3.8, 1.4],
    'ISMIP6BasinBC': [-1.5, 5.7],
    'ISMIP6BasinCCp': [-8.3, 2.3],
    'ISMIP6BasinCpD': [-20.8, 2.3],
    'ISMIP6BasinDDp': [-7.2, 1.9],
    'ISMIP6BasinDpE': [-2.3, 0.4],
    'ISMIP6BasinEF': [33.2, 7.1],
    'ISMIP6BasinFG': [-18.4, 2.7],
    'ISMIP6BasinGH': [-100.6, 3.4],
    'ISMIP6BasinHHp': [-15.5, 1.3],
    'ISMIP6BasinHpI': [-6.4, 2.2],
    'ISMIP6BasinIIpp': [-27.8, 4.4],
    'ISMIP6BasinIppJ': [0.0, 1.5],
    'ISMIP6BasinJK': [2.4, 14.5],
    'ISMIP6BasinKA': [3.0, 1.3],
}

# Copied directly from plot_regionalStats2.py's OBSERVATIONAL_DATASETS['basin_names'].
BASIN_NAMES = {
    'ISMIP6BasinAAp': 'Dronning Maud Land',
    'ISMIP6BasinApB': 'Enderby Land',
    'ISMIP6BasinBC': 'Amery-Lambert',
    'ISMIP6BasinCCp': 'Phillipi, Denman',
    'ISMIP6BasinCpD': 'Totten',
    'ISMIP6BasinDDp': 'Mertz',
    'ISMIP6BasinDpE': 'Victoria Land',
    'ISMIP6BasinEF': 'Ross',
    'ISMIP6BasinFG': 'Getz',
    'ISMIP6BasinGH': 'Thwaites/PIG',
    'ISMIP6BasinHHp': 'Bellingshausen',
    'ISMIP6BasinHpI': 'George VI',
    'ISMIP6BasinIIpp': 'Larsen A-C',
    'ISMIP6BasinIppJ': 'Larsen E',
    'ISMIP6BasinJK': 'FRIS',
    'ISMIP6BasinKA': 'Brunt-Stancomb',
}


# ============================================================================
# Helpers adapted from plot_regional_velo_haf_diffs.py
# ============================================================================

def sanitize_filename_component(label):
    """
    Make a region label safe to embed in a filename by replacing anything
    that isn't alphanumeric, a dash, or an underscore with an underscore.
    """

    return re.sub(r"[^A-Za-z0-9_-]+", "_", label.strip())


def boundary_segments(ds_mesh, mask):
    """
    Build line segments along cell-cell boundaries where mask changes
    from True to False. Uses verticesOnCell / cellsOnCell connectivity from
    the MPAS mesh. Copied directly from plot_regional_velo_haf_diffs.py.
    """

    mask = np.asarray(mask, dtype=bool)

    cells_on_edge = ds_mesh["cellsOnEdge"].values.astype(int)
    vertices_on_edge = ds_mesh["verticesOnEdge"].values.astype(int)

    x_vertex = ds_mesh["xVertex"].values
    y_vertex = ds_mesh["yVertex"].values

    i_cell = cells_on_edge[:, 0]
    j_cell = cells_on_edge[:, 1]

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
    from matplotlib.collections import LineCollection

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


def alnum_region_key(raw_name):
    """
    Convert a raw region name string (e.g. "ISMIP6 Basin A-Ap") to an
    alphanumeric-only key (e.g. "ISMIP6BasinAAp"), matching the convention
    used in plot_regionalStats2.py to look names up in BASIN_NAMES /
    RIGNOT_2019_NET_MB.
    """

    return ''.join(filter(str.isalnum, raw_name))


# ============================================================================
# streams.landice parsing
# ============================================================================

def parse_stream_filename(run_dir, stream_name):
    """
    Read streams.landice (XML) in run_dir and return the absolute path
    (resolved relative to run_dir) of the filename_template for the stream
    (either <stream> or <immutable_stream>) with the given name attribute.
    """

    streams_path = os.path.join(run_dir, "streams.landice")
    tree = ET.parse(streams_path)
    root = tree.getroot()
    for tag in ("immutable_stream", "stream"):
        for elem in root.findall(tag):
            if elem.get("name") == stream_name:
                filename_template = elem.get("filename_template")
                if filename_template is None:
                    raise ValueError(f"Stream '{stream_name}' in {streams_path} has no filename_template")
                return os.path.normpath(os.path.join(run_dir, filename_template))
    raise ValueError(f"Stream named '{stream_name}' not found in {streams_path}")


def find_run_dirs(parent_dir):
    """
    Return a sorted list of subdirectories of parent_dir that look like MALI
    run directories (containing streams.landice and output/regionalStats.nc).
    """

    run_dirs = []
    for entry in sorted(os.listdir(parent_dir)):
        full_path = os.path.join(parent_dir, entry)
        if not os.path.isdir(full_path):
            continue
        streams_file = os.path.join(full_path, "streams.landice")
        stats_file = os.path.join(full_path, REGIONAL_STATS_RELPATH)
        if os.path.isfile(streams_file) and os.path.isfile(stats_file):
            run_dirs.append(full_path)
        else:
            print(f"Skipping '{entry}': not a run directory "
                  "(missing streams.landice and/or output/regionalStats.nc)")
    return run_dirs


# ============================================================================
# Shared mesh / region setup (resolved once, from the first run)
# ============================================================================

def load_shared_mesh_and_regions(first_run_dir):
    mesh_file = parse_stream_filename(first_run_dir, "input")
    regions_file = parse_stream_filename(first_run_dir, "regionsInput")
    print(f"Using shared mesh/IC file: {mesh_file}")
    print(f"Using shared region-mask file: {regions_file}")

    ds_init = xr.open_dataset(mesh_file)
    ds_regions = xr.open_dataset(regions_file)

    descriptor_pristine = mosaic.Descriptor(ds_init, use_latlon=False)

    bed = ds_init["bedTopography"].isel(Time=0)
    obs_u = ds_init["observedSurfaceVelocityX"].isel(Time=0)
    obs_v = ds_init["observedSurfaceVelocityY"].isel(Time=0)
    obs_speed_vals = (np.sqrt(obs_u ** 2 + obs_v ** 2) * SEC_PER_YEAR).values

    n_regions = ds_regions.sizes["nRegions"]
    raw_names_arr = ds_regions["regionNames"].values
    region_names_raw = [
        str(value.decode() if isinstance(value, bytes) else value).strip() for value in raw_names_arr
    ]
    region_keys = [alnum_region_key(name) for name in region_names_raw]

    # dims are (nCells, nRegions); transpose to (nRegions, nCells) for easy
    # per-region indexing below.
    region_cell_masks = ds_regions["regionCellMasks"].values.astype(bool).T

    return {
        "ds_init": ds_init,
        "descriptor_pristine": descriptor_pristine,
        "bed_vals": bed.values,
        "obs_speed_vals": obs_speed_vals,
        "n_regions": n_regions,
        "region_names_raw": region_names_raw,
        "region_keys": region_keys,
        "region_cell_masks": region_cell_masks,
        "x_cell": ds_init["xCell"].values,
        "y_cell": ds_init["yCell"].values,
    }


def build_region_mesh_info(mesh_shared, region_idx):
    """
    Cull the shared pristine mesh descriptor down to a single region. This is
    done once per region and reused across all runs, since all runs share
    the same mesh.
    """

    region_mask = mesh_shared["region_cell_masks"][region_idx]
    descriptor = copy.deepcopy(mesh_shared["descriptor_pristine"])
    index_to_cell_id = None
    if np.any(~region_mask):
        cells_to_cull = ~region_mask
        descriptor.ds = mosaic.utils.cull_mesh(descriptor.ds, cells_to_cull)
        index_to_cell_id = descriptor.ds["indexToCellID"].values

    bed_mesh = (
        mesh_shared["bed_vals"][index_to_cell_id]
        if index_to_cell_id is not None
        else mesh_shared["bed_vals"]
    )

    x_region = mesh_shared["x_cell"][region_mask]
    y_region = mesh_shared["y_cell"][region_mask]
    xmin, xmax = float(x_region.min()), float(x_region.max())
    ymin, ymax = float(y_region.min()), float(y_region.max())
    dx = xmax - xmin
    dy = ymax - ymin
    pad_x = 0.03 * dx if dx > 0.0 else 1000.0
    pad_y = 0.03 * dy if dy > 0.0 else 1000.0

    return {
        "region_mask": region_mask,
        "descriptor": descriptor,
        "index_to_cell_id": index_to_cell_id,
        "bed_mesh": bed_mesh,
        "bounds": (xmin - pad_x, xmax + pad_x, ymin - pad_y, ymax + pad_y),
    }


# ============================================================================
# Per-run data loading
# ============================================================================

def load_run_data(run_dir, mesh_shared):
    run_label = os.path.basename(os.path.normpath(run_dir))
    print(f"  Loading run '{run_label}'...")

    # --- Regional time series (from plot_regionalStats2.py's unit handling) ---
    ds_stats = xr.open_dataset(os.path.join(run_dir, REGIONAL_STATS_RELPATH))
    n_regions_local = ds_stats.sizes["nRegions"]
    if n_regions_local != mesh_shared["n_regions"]:
        sys.exit(
            f"ERROR: run '{run_label}' has {n_regions_local} regions in regionalStats.nc, "
            f"but the shared region-mask file has {mesh_shared['n_regions']}. "
            "All runs must use the same mesh/regions."
        )

    yr = ds_stats["daysSinceStart"].values / 365.0
    vol_ground = ds_stats["regionalGroundedIceVolume"].values * GT_PER_M3_ICE
    vol_ground_chg = vol_ground - vol_ground[0, :]
    vaf = ds_stats["regionalVolumeAboveFloatation"].values * GT_PER_M3_ICE
    vaf_chg = vaf - vaf[0, :]

    # --- 2-D velocity/thickness snapshot (first Time index) ---
    ds_2d = xr.open_dataset(os.path.join(run_dir, OUTPUT_2D_RELPATH))
    h2 = ds_2d["thickness"].isel(Time=0)
    model_speed = ds_2d["surfaceSpeed"].isel(Time=0) * SEC_PER_YEAR
    date2 = xtime_ymd(ds_2d, 0)

    bed_vals = mesh_shared["bed_vals"]
    h2_vals = h2.values
    model_speed_vals = model_speed.values
    obs_speed_vals = mesh_shared["obs_speed_vals"]

    flotation_full = bed_vals + (RHO_I / RHO_W) * h2_vals
    grounded_mask_full = (h2_vals > 0.0) & (flotation_full > 0.0)
    floating_mask_full = (h2_vals > 0.0) & (flotation_full <= 0.0)

    speed_diff_vals = model_speed_vals - obs_speed_vals
    speed_diff_vals = np.where(h2_vals > 0.0, speed_diff_vals, np.nan)

    valid_base = (
        np.isfinite(obs_speed_vals)
        & np.isfinite(model_speed_vals)
        & (obs_speed_vals > 0.0)
        & (model_speed_vals > 0.0)
    )

    return {
        "label": run_label,
        "yr": yr,
        "vol_ground_chg": vol_ground_chg,
        "vaf_chg": vaf_chg,
        "h2_vals": h2_vals,
        "model_speed_vals": model_speed_vals,
        "speed_diff_vals": speed_diff_vals,
        "grounded_mask_full": grounded_mask_full,
        "floating_mask_full": floating_mask_full,
        "valid_base": valid_base,
        "date2": date2,
    }


# ============================================================================
# Figure generation
# ============================================================================

def _plot_heatmap_panel(ax, fig, obs_speed_vals, model_speed_vals, panel_valid, title):
    if not np.any(panel_valid):
        ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes)
        ax.set_title(title)
        return

    log_obs = np.log10(obs_speed_vals[panel_valid])
    log_model = np.log10(model_speed_vals[panel_valid])
    _, _, _, heatmap_img = ax.hist2d(
        log_obs, log_model, bins=HEATMAP_BINS, cmap="viridis", norm=LogNorm()
    )
    fig.colorbar(heatmap_img, ax=ax, label="Count")
    ax.plot(
        [HEATMAP_AXIS_LO, HEATMAP_AXIS_HI], [HEATMAP_AXIS_LO, HEATMAP_AXIS_HI],
        color="k", linewidth=0.8, linestyle="--",
    )
    ax.set_xlim(HEATMAP_AXIS_LO, HEATMAP_AXIS_HI)
    ax.set_ylim(HEATMAP_AXIS_LO, HEATMAP_AXIS_HI)
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=9)


def make_region_figure(region_idx, region_label, region_key, run_labels, run_colors,
                        run_data_list, mesh_info, obs_speed_vals, args):
    n_runs = len(run_data_list)
    n_rows = 1 + n_runs
    fig, axs = plt.subplots(n_rows, 3, figsize=(15, 3.3 * n_rows))
    if n_rows == 1:
        axs = axs.reshape(1, 3)
    fig.suptitle(f"Region: {region_label}", fontsize=12)

    # --- Row 0: grounded MB / VAF summary across runs, + legend panel ---
    ax_grd, ax_vaf, ax_leg = axs[0, 0], axs[0, 1], axs[0, 2]

    lines = []
    for label, color, run_data in zip(run_labels, run_colors, run_data_list):
        line, = ax_grd.plot(run_data["yr"], run_data["vol_ground_chg"][:, region_idx],
                             color=color, label=label)
        lines.append(line)
        ax_vaf.plot(run_data["yr"], run_data["vaf_chg"][:, region_idx], color=color, label=label)

    obs_handle = None
    if region_key in RIGNOT_2019_NET_MB:
        mn, sig = RIGNOT_2019_NET_MB[region_key]
        yr_ref = run_data_list[0]["yr"]
        obs_handle = ax_grd.fill_between(
            yr_ref, yr_ref * (mn - sig), yr_ref * (mn + sig),
            color="gray", alpha=0.3, label="Rignot 2019 net obs",
        )
        ax_vaf.fill_between(
            yr_ref, yr_ref * (mn - sig), yr_ref * (mn + sig),
            color="gray", alpha=0.3, label="Rignot 2019 net obs",
        )

    ax_grd.set_xlabel("Year")
    ax_grd.set_ylabel("Grounded volume change (Gt)")
    ax_grd.set_title("Grounded mass balance")
    ax_grd.grid(True)

    ax_vaf.set_xlabel("Year")
    ax_vaf.set_ylabel("VAF change (Gt)")
    ax_vaf.set_title("Volume above flotation")
    ax_vaf.grid(True)

    ax_leg.axis("off")
    handles = list(lines)
    labels = list(run_labels)
    if obs_handle is not None:
        handles.append(obs_handle)
        labels.append("Rignot 2019 net obs")
    ax_leg.legend(handles, labels, loc="center", fontsize=9, frameon=False)

    # --- Rows 1..N: per-run grounded/floating heatmaps + speed-diff map ---
    region_mask = mesh_info["region_mask"]
    xmin, xmax, ymin, ymax = mesh_info["bounds"]
    index_to_cell_id = mesh_info["index_to_cell_id"]

    for i, (label, color, run_data) in enumerate(zip(run_labels, run_colors, run_data_list)):
        row = i + 1
        ax_grd_hm, ax_flt_hm, ax_map = axs[row, 0], axs[row, 1], axs[row, 2]

        valid = run_data["valid_base"] & region_mask
        _plot_heatmap_panel(
            ax_grd_hm, fig, obs_speed_vals, run_data["model_speed_vals"],
            valid & run_data["grounded_mask_full"], f"{label}: grounded",
        )
        _plot_heatmap_panel(
            ax_flt_hm, fig, obs_speed_vals, run_data["model_speed_vals"],
            valid & run_data["floating_mask_full"], f"{label}: floating",
        )
        ax_grd_hm.set_ylabel("log10(model speed) [m yr$^{-1}$]", fontsize=8)
        ax_flt_hm.set_ylabel("", fontsize=8)
        if row == n_rows - 1:
            ax_grd_hm.set_xlabel("log10(obs speed) [m yr$^{-1}$]", fontsize=8)
            ax_flt_hm.set_xlabel("log10(obs speed) [m yr$^{-1}$]", fontsize=8)

        h2_mesh = (
            run_data["h2_vals"][index_to_cell_id] if index_to_cell_id is not None else run_data["h2_vals"]
        )
        bed_mesh = mesh_info["bed_mesh"]
        cmap_speed = plt.get_cmap("RdBu_r").copy()
        cmap_speed.set_over("magenta")
        cmap_speed.set_under("purple")
        pc = mosaic.polypcolor(
            ax_map, mesh_info["descriptor"], run_data["speed_diff_vals"],
            cmap=cmap_speed, vmin=-SPEED_DIFF_MAX_ABS, vmax=SPEED_DIFF_MAX_ABS, edgecolors="none",
        )
        plot_geometry(ax_map, mesh_info["descriptor"].ds, h2_mesh, bed_mesh,
                       edge_color="b", gl_color="g", edge_ls="-", gl_ls="-")
        ax_map.set_xlim(xmin, xmax)
        ax_map.set_ylim(ymin, ymax)
        ax_map.set_aspect("equal")
        ax_map.set_title(f"{label}: model - obs speed ({run_data['date2']})", fontsize=9)
        fig.colorbar(pc, ax=ax_map, label="Speed diff [m yr$^{-1}$]", extend="both")

    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out_name = os.path.join(
        OUTPUT_DIR, f"region{region_idx:02d}_{sanitize_filename_component(region_label)}.png"
    )
    fig.savefig(out_name, dpi=args.dpi)
    plt.close(fig)
    print(f"  Wrote {out_name}")


def make_summary_figure(run_labels, run_colors, run_data_list, args):
    fig, (ax_grd, ax_vaf) = plt.subplots(1, 2, figsize=(12, 5.5))
    fig.suptitle("Whole ice sheet summary (all regions combined)", fontsize=12)

    for label, color, run_data in zip(run_labels, run_colors, run_data_list):
        ax_grd.plot(run_data["yr"], run_data["vol_ground_chg"].sum(axis=1), color=color, label=label)
        ax_vaf.plot(run_data["yr"], run_data["vaf_chg"].sum(axis=1), color=color, label=label)

    mn_tot = sum(v[0] for v in RIGNOT_2019_NET_MB.values())
    sig_tot = float(np.sqrt(sum(v[1] ** 2 for v in RIGNOT_2019_NET_MB.values())))
    yr_ref = run_data_list[0]["yr"]
    ax_grd.fill_between(
        yr_ref, yr_ref * (mn_tot - sig_tot), yr_ref * (mn_tot + sig_tot),
        color="gray", alpha=0.3, label="Rignot 2019 net obs (AIS total)",
    )
    ax_vaf.fill_between(
        yr_ref, yr_ref * (mn_tot - sig_tot), yr_ref * (mn_tot + sig_tot),
        color="gray", alpha=0.3, label="Rignot 2019 net obs (AIS total)",
    )

    ax_grd.set_xlabel("Year")
    ax_grd.set_ylabel("Grounded volume change (Gt)")
    ax_grd.set_title("Grounded mass balance")
    ax_grd.grid(True)
    ax_grd.legend(fontsize=8)

    ax_vaf.set_xlabel("Year")
    ax_vaf.set_ylabel("VAF change (Gt)")
    ax_vaf.set_title("Volume above flotation")
    ax_vaf.grid(True)
    ax_vaf.legend(fontsize=8)

    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out_name = os.path.join(OUTPUT_DIR, "whole_ice_sheet_summary.png")
    fig.savefig(out_name, dpi=args.dpi)
    plt.close(fig)
    print(f"Wrote {out_name}")


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--dpi", type=int, default=75,
        help="DPI for saved figures (default: %(default)s, favoring fast generation)",
    )
    args = parser.parse_args()

    run_dirs = find_run_dirs(".")
    if len(run_dirs) == 0:
        sys.exit(
            "ERROR: no run subdirectories found (each must contain streams.landice "
            "and output/regionalStats.nc)"
        )
    if len(run_dirs) > MAX_RUNS:
        sys.exit(
            f"ERROR: found {len(run_dirs)} run subdirectories, but at most {MAX_RUNS} "
            "are supported. Please reduce the number of run directories (e.g. run this "
            "script from a directory containing only a subset of runs)."
        )

    run_labels = [os.path.basename(os.path.normpath(d)) for d in run_dirs]
    print(f"Found {len(run_dirs)} run(s): {', '.join(run_labels)}")

    cmap = plt.get_cmap("tab10")
    run_colors = [cmap(i) for i in range(len(run_dirs))]

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print(f"Resolving shared mesh/region info from first run: {run_labels[0]}")
    mesh_shared = load_shared_mesh_and_regions(run_dirs[0])
    obs_speed_vals = mesh_shared["obs_speed_vals"]

    print("Loading per-run data...")
    run_data_list = [load_run_data(d, mesh_shared) for d in run_dirs]

    n_regions = mesh_shared["n_regions"]
    for r in range(n_regions):
        region_key = mesh_shared["region_keys"][r]
        region_label = BASIN_NAMES.get(region_key, mesh_shared["region_names_raw"][r])
        print(f"Making figure for region {r + 1}/{n_regions}: {region_label}")
        mesh_info = build_region_mesh_info(mesh_shared, r)
        make_region_figure(r, region_label, region_key, run_labels, run_colors,
                            run_data_list, mesh_info, obs_speed_vals, args)

    print("Making whole-ice-sheet summary figure")
    make_summary_figure(run_labels, run_colors, run_data_list, args)

    print(f"Done. Figures written to {OUTPUT_DIR}/")


if __name__ == "__main__":
    main()
