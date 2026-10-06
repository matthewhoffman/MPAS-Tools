#!/usr/bin/env python3
"""
Region-by-region calibration of muFriction based on the modeled/observed
surface-speed ratio.

For each region in ``region_mask_file``:
  1. Grounded-ice cells are binned by log10(modeled speed) (4 bins per log
     cycle by default, overridable via --bins-per-decade), and the
     area-weighted mean modeled/observed speed ratio is computed per bin,
     restricted to each bin's interquartile range of ratio values to
     reduce sensitivity to outliers. By default, "modeled speed"/
     "observed speed" are ``surfaceSpeed``/observed surface speed; with
     ``--use-basal-speed-ratio``, the deformation velocity (modeled
     ``surfaceSpeed - basalSpeed``) is removed from both, so the ratio
     is instead (modeled ``basalSpeed``) / (observed surface speed minus
     modeled deformation velocity).
  2. A curve is fit to log10(area-weighted mean ratio) vs. log10(modeled
     speed) bin centers, using either a low-degree polynomial (``polyfit``,
     the default) or a piecewise-linear interpolant directly through the
     bin means (``piecewise-linear``); see ``--fit-type``.
  3. Each grounded cell's correction factor is the fitted ratio, evaluated at
     that cell's own modeled speed (clipped to the fit's bin range to avoid
     extrapolation blowups), raised to the ``--exponent`` power (default
     1/3).
  4. A height-above-flotation (HAF) cubic-smoothstep taper (the same for
     every region) blends the correction between "none" (exponent 0, no
     change) at/below ``haf_begin`` and "full" (the complete ``^exponent``
     correction) at/above ``haf_end``.

The corrected ``muFriction`` field is written to a new copy of
``mesh_file``. A per-region two-panel diagnostic PNG is also produced:
(left) the heatmap-with-fit, showing the area-weighted bin means as well
as the fitted curve, and (right) a spatial map of the mu ratio (i.e.
after/before), using a fixed log-scale colorbar range of [0.1, 100] with
above/below-range indicators.

Reference for the relevant calculations/variable names:
``plot_regional_velo_haf_diffs.py``.
"""

import argparse
import copy
import os
import re
import subprocess

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import xarray as xr

import mosaic

SEC_PER_YEAR = 365.0 * 24.0 * 60.0 * 60.0
RHO_I = 910.0
RHO_W = 1028.0

# Default bins per log cycle (decade) in log10(modeled-speed) space;
# overridable via the --bins-per-decade CLI argument.
DEFAULT_BINS_PER_DECADE = 4
# Minimum number of valid cells a bin must contain to be used in the fit.
MIN_CELLS_PER_BIN = 10


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


def sanitize_filename_component(label):
    """
    Make a region label safe to embed in a filename by replacing anything
    that isn't alphanumeric, a dash, or an underscore with an underscore.
    """

    return re.sub(r"[^A-Za-z0-9_-]+", "_", label.strip())


def get_region_names(ds_regions):
    """
    Return a list of region names decoded from ``regionNames`` or
    ``regionMaskNames`` in ``ds_regions``, or ``None`` if neither variable
    is present.
    """

    for name_var in ("regionNames", "regionMaskNames"):
        if name_var not in ds_regions:
            continue
        raw = ds_regions[name_var].values
        if raw.ndim == 1:
            return [str(value.decode() if isinstance(value, bytes) else value).strip() for value in raw]
        names = []
        for row in raw:
            chars = []
            for char in row:
                if isinstance(char, bytes):
                    chars.append(char.decode("utf-8"))
                else:
                    chars.append(str(char))
            names.append("".join(chars).replace("\x00", "").strip())
        return names
    return None


def boundary_segments(ds_mesh, mask):
    """
    Build line segments along cell-cell boundaries where mask changes
    from True to False. Reproduced verbatim from
    plot_regional_velo_haf_diffs.py.
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


def plot_geometry(ax, ds_mesh, thickness, bed, edge_color, gl_color):
    thickness = np.asarray(thickness)
    ice_mask = thickness > 0.0
    flotation = np.asarray(bed) + (RHO_I / RHO_W) * thickness
    grounded_mask = (thickness > 0.0) & (flotation > 0.0)
    ice_segments = boundary_segments(ds_mesh, ice_mask)
    gl_segments = boundary_segments(ds_mesh, grounded_mask)

    from matplotlib.collections import LineCollection

    if len(ice_segments) > 0:
        ax.add_collection(LineCollection(ice_segments, colors=edge_color, linewidths=0.5))
    if len(gl_segments) > 0:
        ax.add_collection(LineCollection(gl_segments, colors=gl_color, linewidths=0.5, linestyles="--"))


def cubic_smoothstep(t):
    """
    Cubic smoothstep, mapping t (clipped to [0, 1]) to 3t^2 - 2t^3.
    """

    t = np.clip(t, 0.0, 1.0)
    return 3.0 * t**2 - 2.0 * t**3


def compute_bin_edges(log_model_valid, bins_per_decade):
    """
    Compute global log10(modeled-speed) bin edges, at bins_per_decade bins
    per decade, spanning the full range of valid data (snapped outward to
    the bin-width grid so every data point falls within some bin).
    """

    bin_width = 1.0 / bins_per_decade
    lo = np.floor(log_model_valid.min() / bin_width) * bin_width
    hi = np.ceil(log_model_valid.max() / bin_width) * bin_width
    if hi <= lo:
        hi = lo + bin_width
    n_bins = int(round((hi - lo) / bin_width))
    return np.linspace(lo, hi, n_bins + 1)


def fit_region_ratio(log_model, ratio, area, bin_edges, poly_degree, fit_type):
    """
    Bin (log_model, ratio) pairs into bin_edges, compute the area-weighted
    mean ratio per bin (restricted to each bin's interquartile range of
    ratio values, to reduce sensitivity to outliers), and fit a curve to
    log10(area-weighted mean ratio) vs. bin center, using only bins with
    at least MIN_CELLS_PER_BIN cells (before IQR filtering).

    ``fit_type`` is either ``"polyfit"`` (fit a ``poly_degree`` polynomial,
    weighted by each usable bin's total IQR-filtered area) or
    ``"piecewise-linear"`` (linearly interpolate directly between the
    usable bins' area-weighted means).

    Returns (fit, bin_centers, bin_mean_log_ratio, bin_counts, fit_lo,
    fit_hi) where ``fit`` is a dict describing the fitted curve (see
    ``evaluate_log_ratio_fit``), or ``None`` if there were not enough
    usable bins to fit. bin_centers/bin_mean_log_ratio/bin_counts cover
    every bin with at least 1 cell (for heatmap/diagnostic plotting),
    regardless of whether that bin was used in the fit.
    """

    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bin_index = np.digitize(log_model, bin_edges) - 1
    n_bins = len(bin_centers)

    bin_mean_log_ratio = np.full(n_bins, np.nan)
    bin_counts = np.zeros(n_bins, dtype=int)
    bin_total_area = np.zeros(n_bins)

    for i in range(n_bins):
        in_bin = bin_index == i
        count = int(np.count_nonzero(in_bin))
        bin_counts[i] = count
        if count > 0:
            bin_ratio = ratio[in_bin]
            bin_area = area[in_bin]
            # Restrict to the bin's interquartile range of ratio values
            # before averaging, to reduce sensitivity to outliers.
            q1, q3 = np.percentile(bin_ratio, [25.0, 75.0])
            iqr_mask = (bin_ratio >= q1) & (bin_ratio <= q3)
            if not np.any(iqr_mask):
                iqr_mask = np.ones_like(bin_ratio, dtype=bool)
            total_area = np.sum(bin_area[iqr_mask])
            bin_total_area[i] = total_area
            bin_mean_log_ratio[i] = np.log10(np.sum(bin_ratio[iqr_mask] * bin_area[iqr_mask]) / total_area)

    usable = bin_counts >= MIN_CELLS_PER_BIN
    n_usable = int(np.count_nonzero(usable))
    min_bins_needed = (poly_degree + 1) if fit_type == "polyfit" else 2
    if n_usable < min_bins_needed:
        return None, bin_centers, bin_mean_log_ratio, bin_counts, None, None

    fit_lo = bin_centers[usable].min()
    fit_hi = bin_centers[usable].max()

    if fit_type == "polyfit":
        weights = np.sqrt(bin_total_area[usable])
        coeffs = np.polyfit(bin_centers[usable], bin_mean_log_ratio[usable], deg=poly_degree, w=weights)
        fit = {"type": "polyfit", "coeffs": coeffs}
    elif fit_type == "piecewise-linear":
        fit = {"type": "piecewise-linear", "xp": bin_centers[usable], "fp": bin_mean_log_ratio[usable]}
    else:
        raise ValueError(f"Unknown fit_type '{fit_type}'")

    return fit, bin_centers, bin_mean_log_ratio, bin_counts, fit_lo, fit_hi


def evaluate_log_ratio_fit(fit, log_x):
    """
    Evaluate a fit produced by ``fit_region_ratio`` (either a polynomial or
    a piecewise-linear interpolant) at the given log10(modeled speed)
    value(s). Values outside the fit's bin range should be clipped by the
    caller first to avoid extrapolation; ``np.interp`` clamps to the
    endpoint values automatically for the piecewise-linear case.
    """

    if fit["type"] == "polyfit":
        return np.polyval(fit["coeffs"], log_x)
    elif fit["type"] == "piecewise-linear":
        return np.interp(log_x, fit["xp"], fit["fp"])
    raise ValueError(f"Unknown fit type '{fit['type']}'")


def plot_region_summary(
    region_label, region_suffix, log_obs, log_model, fit, bin_centers, bin_mean_log_ratio, bin_counts, fit_lo, fit_hi, date_label,
    descriptor, mu_before_mesh, mu_after_mesh, h_mesh, bed_mesh, xmin, xmax, ymin, ymax, dpi, speed_label,
):
    """
    Produce a single two-panel PNG per region: (left) the modeled-vs-
    observed speed heatmap with the area-weighted bin means and fitted
    ratio curve, and (right) a spatial map of the muFriction ratio
    (after / before correction), using a fixed log-scale colorbar range
    of [0.1, 100] with above/below-range indicators.
    """

    fig, (ax_heatmap, ax_map) = plt.subplots(1, 2, figsize=(16, 8), constrained_layout=True)

    # --- Left panel: heatmap with fitted curve ---
    axis_lo, axis_hi = -2.0, 5.0
    heatmap_bins = np.linspace(axis_lo, axis_hi, 141)

    _, _, _, heatmap_img = ax_heatmap.hist2d(log_obs, log_model, bins=heatmap_bins, cmap="viridis", norm=LogNorm())
    fig.colorbar(heatmap_img, ax=ax_heatmap, label="Count")
    ax_heatmap.plot([axis_lo, axis_hi], [axis_lo, axis_hi], color="k", linewidth=0.8, linestyle="--", label="1:1")

    # Area-weighted mean ratio per populated bin, converted to the same
    # (log_obs, log_model) heatmap coordinates as the fitted curve below:
    # log_obs = log_model - log10(ratio).
    populated = bin_counts > 0
    if np.any(populated):
        y_bins = bin_centers[populated]
        x_bins = y_bins - bin_mean_log_ratio[populated]
        ax_heatmap.plot(
            x_bins, y_bins, linestyle="none", marker="o", markersize=5,
            markerfacecolor="orange", markeredgecolor="k", label="Area-weighted bin mean",
        )

    if fit is not None:
        # The fit is log10(ratio) vs. log10(modeled speed), i.e.
        # y_fit (log_model) is the independent variable; invert
        # ratio = model/obs to get the corresponding log10(observed speed)
        # for each point on the curve: log_obs = log_model - log10(ratio).
        y_fit = np.linspace(fit_lo, fit_hi, 100)
        x_fit = y_fit - evaluate_log_ratio_fit(fit, y_fit)
        ax_heatmap.plot(x_fit, y_fit, color="r", linewidth=1.5, label="Fitted ratio")

    ax_heatmap.set_xlim(axis_lo, axis_hi)
    ax_heatmap.set_ylim(axis_lo, axis_hi)
    ax_heatmap.set_aspect("equal")
    ax_heatmap.set_xlabel(f"log10(observed {speed_label}) [log10(m yr$^{{-1}}$)]")
    ax_heatmap.set_ylabel(f"log10(modeled {speed_label}) [log10(m yr$^{{-1}}$)]")
    ax_heatmap.set_title(f"Modeled vs. observed {speed_label} (grounded ice): {date_label}")
    ax_heatmap.legend(loc="best")

    # --- Right panel: spatial map of the muFriction ratio ---
    dx = xmax - xmin
    dy = ymax - ymin
    pad_x = 0.03 * dx if dx > 0.0 else 1000.0
    pad_y = 0.03 * dy if dy > 0.0 else 1000.0

    mu_ratio_mesh = np.divide(
        mu_after_mesh, mu_before_mesh, out=np.ones_like(mu_after_mesh), where=mu_before_mesh != 0.0
    )

    # Fixed log-scale diverging colormap range centered on 1.0 (no
    # change), with above/below-range indicators for values outside
    # [0.1, 100].
    vmin, vmax = 0.1, 100.0

    pc = mosaic.polypcolor(
        ax_map, descriptor, mu_ratio_mesh, cmap="RdBu_r", norm=LogNorm(vmin=vmin, vmax=vmax), edgecolors="none"
    )
    plot_geometry(ax_map, descriptor.ds, h_mesh, bed_mesh, edge_color="k", gl_color="r")
    ax_map.set_xlim(xmin - pad_x, xmax + pad_x)
    ax_map.set_ylim(ymin - pad_y, ymax + pad_y)
    ax_map.set_aspect("equal")
    ax_map.set_xlabel("x [m]")
    ax_map.set_ylabel("y [m]")
    ax_map.set_title("muFriction ratio (after / before speed-ratio correction)")
    fig.colorbar(pc, ax=ax_map, label="muFriction ratio (after / before)", extend="both")

    fig.suptitle(f"Region: {region_label}")

    filename = f"mu_ratio_summary{region_suffix}.png"
    fig.savefig(filename, dpi=dpi)
    plt.close(fig)
    print(f"Wrote {filename}")


def main():
    parser = argparse.ArgumentParser(
        description="Calibrate muFriction region-by-region based on the modeled/observed surface-speed ratio, "
        "with a height-above-flotation taper."
    )

    parser.add_argument("mesh_file", help="MALI mesh/IC file with mesh info, observed velocity, bedTopography, thickness, and muFriction.")
    parser.add_argument(
        "solution_file",
        help="MALI output file with surfaceSpeed (time level 0 is used); also requires basalSpeed if "
        "--use-basal-speed-ratio is set.",
    )
    parser.add_argument("haf_begin", type=float, help="Height above flotation [m] at/below which the correction is fully tapered off (0).")
    parser.add_argument("haf_end", type=float, help="Height above flotation [m] at/above which the correction is fully applied (1).")
    parser.add_argument("region_mask_file", help="Region mask file with regionCellMasks; every region is processed.")
    parser.add_argument(
        "output_file",
        help="Output NetCDF file: a copy of mesh_file with a corrected muFriction field. "
        "Suggested naming convention: <mesh_file stem>_mucorrected.nc",
    )
    parser.add_argument(
        "--poly-degree",
        type=int,
        default=3,
        help="Degree of the polynomial fit to log10(ratio) vs. log10(modeled speed) (default: %(default)s). "
        "Ignored if --fit-type is piecewise-linear.",
    )
    parser.add_argument(
        "--fit-type",
        choices=("polyfit", "piecewise-linear"),
        default="polyfit",
        help="Type of curve fit to log10(area-weighted mean ratio) vs. log10(modeled speed) bin centers: "
        "a 'polyfit' polynomial (default) or a 'piecewise-linear' interpolant directly through the bin means.",
    )
    parser.add_argument(
        "--bins-per-decade",
        type=int,
        default=DEFAULT_BINS_PER_DECADE,
        help="Number of speed bins per log10 cycle (decade) used for the ratio fit and heatmap bin means "
        "(default: %(default)s).",
    )
    parser.add_argument(
        "--exponent",
        type=float,
        default=1.0 / 3.0,
        help="Power to which the fitted speed ratio is raised to produce the full (fully-tapered-in) muFriction "
        "correction factor (default: %(default)s, i.e. 1/3).",
    )
    parser.add_argument(
        "--dpi",
        type=float,
        default=300,
        help="DPI used when saving all figures (default: %(default)s).",
    )
    parser.add_argument(
        "--use-basal-speed-ratio",
        action="store_true",
        help="Use (modeled basalSpeed) / (observed surface speed minus modeled deformation velocity) as the "
        "ratio instead of (modeled surfaceSpeed) / (observed surface speed), where the modeled deformation "
        "velocity is surfaceSpeed - basalSpeed. Requires basalSpeed in solution_file. Cells where the "
        "deformation-removed observed speed is <= 0 are excluded from the analysis.",
    )

    args = parser.parse_args()

    if not (args.haf_begin < args.haf_end):
        parser.error("haf_begin must be less than haf_end")

    ds_mesh = xr.open_dataset(args.mesh_file)
    ds_soln = xr.open_dataset(args.solution_file)
    ds_regions = xr.open_dataset(args.region_mask_file)

    thickness = at_time(ds_mesh["thickness"], 0)
    bed = at_time(ds_mesh["bedTopography"], 0)
    mu_old_da = at_time(ds_mesh["muFriction"], 0)
    mu_old = mu_old_da.values
    area = ds_mesh["areaCell"].values

    obs_u = at_time(ds_mesh["observedSurfaceVelocityX"], 0)
    obs_v = at_time(ds_mesh["observedSurfaceVelocityY"], 0)
    obs_surface_speed = (np.sqrt(obs_u**2 + obs_v**2) * SEC_PER_YEAR).values
    model_surface_speed = (at_time(ds_soln["surfaceSpeed"], 0) * SEC_PER_YEAR).values

    if args.use_basal_speed_ratio:
        # Remove the modeled deformation velocity (surfaceSpeed - basalSpeed)
        # from both the modeled and observed speed, so the ratio reflects
        # the sliding (basal) component rather than the full surface speed.
        # The observed counterpart has no independently-observed basal/
        # deformation component, so the modeled deformation velocity is
        # used as the best available estimate; cells where this leaves a
        # negative "observed basal speed" are excluded below (valid_base).
        model_basal_speed = (at_time(ds_soln["basalSpeed"], 0) * SEC_PER_YEAR).values
        deformation_speed = model_surface_speed - model_basal_speed
        model_speed = model_basal_speed
        obs_speed = obs_surface_speed - deformation_speed
        speed_label = "basal speed"
    else:
        model_speed = model_surface_speed
        obs_speed = obs_surface_speed
        speed_label = "surface speed"

    haf = height_above_flotation(thickness.values, bed.values)

    flotation = bed.values + (RHO_I / RHO_W) * thickness.values
    grounded_mask = (thickness.values > 0.0) & (flotation > 0.0)

    valid_base = (
        np.isfinite(obs_speed)
        & np.isfinite(model_speed)
        & (obs_speed > 0.0)
        & (model_speed > 0.0)
        & grounded_mask
    )

    if args.use_basal_speed_ratio:
        n_negative_obs_basal = int(np.count_nonzero(grounded_mask & np.isfinite(obs_speed) & (obs_speed <= 0.0)))
        if n_negative_obs_basal > 0:
            print(
                f"Note: {n_negative_obs_basal} grounded cell(s) have a non-positive deformation-removed "
                "observed speed estimate (observed surface speed <= modeled deformation velocity) and are "
                "excluded from the ratio fit/correction."
            )

    ratio_all = np.divide(model_speed, obs_speed, out=np.ones_like(model_speed), where=obs_speed > 0.0)

    # Global bin edges, shared across all regions, from the full valid-data
    # range of log10(modeled speed).
    log_model_all = np.log10(model_speed[valid_base])
    bin_edges = compute_bin_edges(log_model_all, args.bins_per_decade)

    # Descriptor used for the mu ratio map (built once; a fresh deep copy
    # is culled per-region, since mosaic.utils.cull_mesh mutates
    # descriptor.ds in place).
    descriptor_pristine = mosaic.Descriptor(ds_mesh, use_latlon=False)

    n_regions = ds_regions.sizes["nRegions"]
    names = get_region_names(ds_regions)

    # Correction-factor exponent base: ratio_fit per cell, defaults to 1.0
    # (no change) outside any region.
    ratio_fit_field = np.ones_like(obs_speed)

    date_label = "time level 0"

    # Plot inputs that depend on mu_new (computed after the HAF taper,
    # below) are collected here and rendered once the taper has been
    # applied.
    pending_mu_plots = []

    n_succeeded = 0
    n_skipped = 0
    for region_index in range(n_regions):
        region_label = names[region_index] if names is not None else str(region_index)
        region_suffix = f"_region{region_index:02d}_{sanitize_filename_component(region_label)}"
        try:
            region_mask_da = ds_regions["regionCellMasks"].isel(nRegions=region_index).astype(bool)
            region_mask = region_mask_da.values

            if not np.any(region_mask):
                raise ValueError(f"Region '{region_label}' selects zero cells; skipping.")

            region_valid = valid_base & region_mask
            if not np.any(region_valid):
                raise ValueError(f"Region '{region_label}' has no valid grounded-ice cells; skipping.")

            log_obs_region = np.log10(obs_speed[region_valid])
            log_model_region = np.log10(model_speed[region_valid])
            ratio_region = ratio_all[region_valid]
            area_region = area[region_valid]

            fit, bin_centers, bin_mean_log_ratio, bin_counts, fit_lo, fit_hi = fit_region_ratio(
                log_model_region, ratio_region, area_region, bin_edges, args.poly_degree, args.fit_type
            )

            if fit is None:
                min_bins_needed = (args.poly_degree + 1) if args.fit_type == "polyfit" else 2
                raise ValueError(
                    f"Region '{region_label}' does not have enough populated speed bins "
                    f"(>= {MIN_CELLS_PER_BIN} cells/bin, needs >= {min_bins_needed} such bins "
                    f"for fit type '{args.fit_type}'); skipping."
                )

            # Per-cell correction factor for this region: evaluate the fit
            # at each cell's own modeled speed, clipped to the fit's bin
            # range to avoid extrapolation blowups. Cells with zero/invalid
            # modeled speed (e.g. ice-free cells) are left at the default
            # ratio_fit_field value of 1.0 (no correction) to avoid
            # log10(0)/log10(negative) warnings.
            eval_mask = region_mask & np.isfinite(model_speed) & (model_speed > 0.0)
            log_model_region_full = np.log10(model_speed[eval_mask])
            log_model_clipped = np.clip(log_model_region_full, fit_lo, fit_hi)
            log_ratio_fit = evaluate_log_ratio_fit(fit, log_model_clipped)
            ratio_fit_field[eval_mask] = 10.0**log_ratio_fit

            x_region = ds_mesh.xCell.where(region_mask_da, drop=True)
            y_region = ds_mesh.yCell.where(region_mask_da, drop=True)
            xmin = float(x_region.min())
            xmax = float(x_region.max())
            ymin = float(y_region.min())
            ymax = float(y_region.max())

            # Build the culled mesh/descriptor for the mu ratio map.
            descriptor = copy.deepcopy(descriptor_pristine)
            cells_to_cull = ~region_mask
            descriptor.ds = mosaic.utils.cull_mesh(descriptor.ds, cells_to_cull)
            index_to_cell_id = descriptor.ds["indexToCellID"].values

            h_mesh = thickness.isel(nCells=index_to_cell_id).values
            bed_mesh = bed.isel(nCells=index_to_cell_id).values

            # mu_after is computed below (after the HAF taper is applied
            # globally), so defer the combined heatmap+mu-ratio-map plot
            # until then.
            pending_mu_plots.append(
                dict(
                    region_label=region_label,
                    region_suffix=region_suffix,
                    log_obs_region=log_obs_region,
                    log_model_region=log_model_region,
                    fit=fit,
                    bin_centers=bin_centers,
                    bin_mean_log_ratio=bin_mean_log_ratio,
                    bin_counts=bin_counts,
                    fit_lo=fit_lo,
                    fit_hi=fit_hi,
                    descriptor=descriptor,
                    index_to_cell_id=index_to_cell_id,
                    h_mesh=h_mesh,
                    bed_mesh=bed_mesh,
                    xmin=xmin,
                    xmax=xmax,
                    ymin=ymin,
                    ymax=ymax,
                )
            )
            n_succeeded += 1
        except Exception as exc:  # noqa: BLE001 - intentionally broad: warn and continue
            print(f"Warning: skipping region '{region_label}': {exc}")
            n_skipped += 1
            continue

    # -------------------------------------------------------------------
    # HAF taper (same for every region): blend each cell's correction
    # between "none" (exponent 0) at/below haf_begin and "full" (exponent
    # 1) at/above haf_end, using a cubic smoothstep.
    # -------------------------------------------------------------------
    taper_t = (haf - args.haf_begin) / (args.haf_end - args.haf_begin)
    taper = cubic_smoothstep(taper_t)

    mu_new = mu_old * ratio_fit_field ** (taper * args.exponent)

    # Now that mu_new exists, produce the deferred two-panel summary plots.
    for pending in pending_mu_plots:
        mu_before_mesh = mu_old[pending["index_to_cell_id"]]
        mu_after_mesh = mu_new[pending["index_to_cell_id"]]
        plot_region_summary(
            pending["region_label"],
            pending["region_suffix"],
            pending["log_obs_region"],
            pending["log_model_region"],
            pending["fit"],
            pending["bin_centers"],
            pending["bin_mean_log_ratio"],
            pending["bin_counts"],
            pending["fit_lo"],
            pending["fit_hi"],
            date_label,
            pending["descriptor"],
            mu_before_mesh,
            mu_after_mesh,
            pending["h_mesh"],
            pending["bed_mesh"],
            pending["xmin"],
            pending["xmax"],
            pending["ymin"],
            pending["ymax"],
            args.dpi,
            speed_label,
        )

    if n_regions > 1:
        print(f"Processed {n_succeeded} region(s), skipped {n_skipped} region(s).")

    # -------------------------------------------------------------------
    # Write the corrected muFriction field to a new copy of mesh_file.
    # -------------------------------------------------------------------
    mu_da = ds_mesh["muFriction"]
    if "Time" in mu_da.dims:
        mu_new_out = np.broadcast_to(mu_new, mu_da.shape)
    else:
        mu_new_out = mu_new
    ds_mesh["muFriction"] = xr.DataArray(
        mu_new_out,
        dims=mu_da.dims,
        attrs={
            **mu_da.attrs,
            "long_name": (
                f"{mu_da.attrs.get('long_name', 'basal friction coefficient')} "
                "(corrected by region-wise modeled/observed speed ratio via "
                "calibrate_mu_from_velo_ratio.py)"
            ),
        },
    )
    ds_mesh.load()

    # Writing NETCDF3 directly via xarray/netCDF4 is very slow for large
    # files. Instead, write fast as NETCDF4 to a temp file, then use nco's
    # `ncks` to convert to CDF5 (64-bit data, "pnetcdf") format, which MPAS
    # requires. (Same pattern as convert_budd_N_ocean_to_transition2.py.)
    tmp_nc4 = args.output_file + ".nc4.tmp"
    ds_mesh.to_netcdf(tmp_nc4, format="NETCDF4")
    ds_mesh.close()

    try:
        subprocess.run(["ncks", "-O", "--fl_fmt=64bit_data", tmp_nc4, args.output_file], check=True)
    finally:
        os.remove(tmp_nc4)

    print(f"Wrote corrected mesh file: {args.output_file}")


if __name__ == "__main__":
    main()
