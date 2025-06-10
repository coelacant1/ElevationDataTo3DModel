#!/usr/bin/env python3
"""GeoTIFFToSTL.py - Convert a GeoTIFF DEM/bathymetry raster to an STL mesh.

Features
--------
* Down-sampling & Gaussian blur for large rasters
* Datum shift so chosen `--water-level` becomes Z=0
* **--align-shoreline**: shift *land only* so the lowest shoreline cell is 0m
  → removes the vertical cliff while preserving underwater relief
* Optional land/water vertical re-balancing via `--land-ratio`
* Uniform or per-axis scaling, optional centring at origin
* `--debug` prints extensive Z-stats at each pipeline stage

Quick example
-------------
    python GeoTIFFToSTL.py input.tif output.stl \
        --water-level 0 \
        --align-shoreline \
        --decimate 4 --center --scale 100 --debug
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import Affine
from scipy.ndimage import (distance_transform_edt, gaussian_filter, zoom)
import trimesh

# =============================================================================
# Raster helpers
# =============================================================================

def load_elevation(path: str):
    src = rasterio.open(path)
    Z = src.read(1).astype(np.float32)
    return Z, src.transform, src.nodata


def decimate_array(Z: np.ndarray, factor: int):
    if factor <= 1:
        return Z
    return zoom(Z, 1 / factor, order=1)


def decimated_transform(t: Affine, factor: int):
    if factor <= 1:
        return t
    return Affine(t.a * factor, t.b, t.c, t.d, t.e * factor, t.f)

# =============================================================================
# Vertical processing (shoreline, land/water balance)
# =============================================================================

def process_vertical(
    Z: np.ndarray,
    *,
    water_level: float,
    nodata: float | None,
    align_shoreline: bool = False,
    raise_water: bool = False,
    land_ratio: float | None = None,
    blend_px: int = 0,
    curve: str = "linear",
) -> np.ndarray:
    # """Return new Z array with shoreline datum aligned and optional tweaks."""
    Z = Z.astype(np.float32).copy()

    # 1 water surface becomes 0
    Z -= water_level

    # 2 replace nodata under water with 0
    if nodata is not None:
        water_mask0 = (Z <= 0) | (Z == nodata)
        Z[water_mask0 & (Z == nodata)] = 0.0

    water_mask = Z <= 0
    land_mask = ~water_mask

        # 3 datum mismatch fix
    if align_shoreline and raise_water:
        raise ValueError("--align-shoreline and --raise-water are mutually exclusive")

    if align_shoreline and land_mask.any():
        # pull LAND down so its lowest pixel sits at 0
        shoreline_level = np.nanmin(Z[land_mask])  # could be >0
        Z[land_mask] -= shoreline_level
        water_mask = Z <= 0
        land_mask = ~water_mask

    if raise_water and water_mask.any():
        # lift WATER up so its highest pixel touches first positive land value
        # find first significant positive land pixel (>= 1 m) to ignore numeric noise
        positive_land = Z[land_mask & (Z > 1)]
        if positive_land.size:
            offset = np.nanmin(positive_land)
        else:
            offset = np.nanmin(Z[land_mask]) or 0.0
        Z[water_mask] += offset
        # after lifting, recalc masks
        water_mask = Z <= offset * 0.001  # allow tiny epsilon around 0
        land_mask = ~water_mask

    # 4 - optional land/water re-balance
    if land_ratio is not None and 0 < land_ratio < 1 and land_mask.any() and water_mask.any():
        old_land = np.nanmax(Z[land_mask]) or 1e-6
        old_water = abs(np.nanmin(Z[water_mask])) or 1e-6
        total = old_land + old_water
        new_land = land_ratio * total
        new_water = (1 - land_ratio) * total
        Z[land_mask]  *= new_land / old_land
        Z[water_mask] *= new_water / old_water

    # 5 - optional feather of first land pixels
    if blend_px > 0 and land_mask.any():
        dist = distance_transform_edt(water_mask)
        ramp = np.clip(dist / blend_px, 0, 1)
        if curve == "cosine":
            ramp = 0.5 - 0.5 * np.cos(np.pi * ramp)
        Z[land_mask] *= ramp[land_mask]
        Z = gaussian_filter(Z, sigma=0.6)

    return Z

# =============================================================================
# Mesh helpers
# =============================================================================

def to_mesh(Z: np.ndarray, transform: Affine) -> trimesh.Trimesh:
    rows, cols = Z.shape
    xs, ys = np.meshgrid(np.arange(cols), np.arange(rows))
    x_geo, y_geo = rasterio.transform.xy(transform, ys, xs, offset="center")
    verts = np.column_stack([np.asarray(x_geo).ravel(),
                             np.asarray(y_geo).ravel(),
                             Z.ravel()])
    faces = []
    for r in range(rows - 1):
        r0 = r * cols
        r1 = (r + 1) * cols
        for c in range(cols - 1):
            v0 = r0 + c
            v1 = v0 + 1
            v2 = r1 + c
            v3 = v2 + 1
            faces.extend([[v0, v2, v1], [v1, v2, v3]])
    faces = np.asarray(faces, dtype=np.int64)
    return trimesh.Trimesh(vertices=verts, faces=faces, process=False)


def scale_mesh(mesh: trimesh.Trimesh, *, scale=None,
               scale_x=None, scale_y=None, scale_z=None):
    if scale is not None:
        mesh.apply_scale(scale)
    else:
        mesh.apply_scale([
            1.0 if scale_x is None else scale_x,
            1.0 if scale_y is None else scale_y,
            1.0 if scale_z is None else scale_z,
        ])


def center_mesh(mesh: trimesh.Trimesh, at=(0.0, 0.0, 0.0)):
    mesh.apply_translation(np.asarray(at) - mesh.bounds.mean(axis=0))

# =============================================================================
# Debug helper
# =============================================================================

def debug_stats(stage: str, Z: np.ndarray, water_level: float):
    zmin, zmax = float(np.nanmin(Z)), float(np.nanmax(Z))
    water_mask = Z <= water_level
    water_cnt = int(np.count_nonzero(water_mask))
    land_cnt = Z.size - water_cnt
    w_min = float(np.nanmin(Z[water_mask])) if water_cnt else float('nan')
    w_max = float(np.nanmax(Z[water_mask])) if water_cnt else float('nan')
    l_min = float(np.nanmin(Z[~water_mask])) if land_cnt else float('nan')
    l_max = float(np.nanmax(Z[~water_mask])) if land_cnt else float('nan')

    print(f"[DEBUG] {stage}: global({zmin:.3f},{zmax:.3f}) wl={water_level:.3f} | "
          f"water px={water_cnt} ({w_min:.3f},{w_max:.3f}) | "
          f"land px={land_cnt} ({l_min:.3f},{l_max:.3f})")
    if water_level < zmin or water_level > zmax:
        print("[WARN ] water_level outside DEM range - no water/land split.")

# =============================================================================
# CLI parsing
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(description="Convert GeoTIFF DEM to STL mesh")
    p.add_argument("input", help="GeoTIFF path")
    p.add_argument("output", help="Output STL/OBJ/PLY path")

    p.add_argument("--decimate", type=int, default=1, help="Down-sample factor")
    p.add_argument("--sigma", type=float, default=0.0, help="Gaussian blur (px)")

    p.add_argument("--water-level", type=float, required=True, help="Water surface elevation in DEM units")
    p.add_argument("--align-shoreline", action="store_true", help="Lower LAND so its lowest shoreline cell is 0m")
    p.add_argument("--raise-water", action="store_true", help="Raise WATER so its highest pixel meets the first positive land cell")
    p.add_argument("--blend", type=int, default=0, help="Feather width in px (land only)")
    p.add_argument("--blend-curve", choices=["linear", "cosine"], default="linear")
    p.add_argument("--land-ratio", type=float, default=None, help="0-1 fraction of height to allocate to land after squash")

    p.add_argument("--center", action="store_true", help="Center mesh at origin")
    p.add_argument("--scale", type=float, default=None, help="Uniform scale")
    p.add_argument("--scale-x", type=float, default=None)
    p.add_argument("--scale-y", type=float, default=None)
    p.add_argument("--scale-z", type=float, default=None)

    p.add_argument("--debug", action="store_true", help="Print debug stats")
    return p.parse_args()

# =============================================================================
# Main
# =============================================================================

def main():
    a = parse_args()

    # 1 load DEM
    Z, transform, nodata = load_elevation(a.input)
    if a.debug:
        debug_stats("raw", Z, a.water_level)

    # 2 decimate & blur
    Z = decimate_array(Z, a.decimate)
    if a.sigma > 0:
        Z = gaussian_filter(Z, sigma=a.sigma)
    transform = decimated_transform(transform, a.decimate)
    if a.debug:
        debug_stats("after decimate/blur", Z, a.water_level)

    # 3 vertical adjustments (shoreline, align, etc.)
        Z = process_vertical(
        Z,
        water_level=a.water_level,
        nodata=nodata,
        align_shoreline=a.align_shoreline,
        raise_water=a.raise_water,
        land_ratio=a.land_ratio,
        blend_px=a.blend,
        curve=a.blend_curve,
    )
    if a.debug:
        debug_stats("after vertical process", Z, 0.0)  # water now at 0

    # 4 mesh generation
    mesh = to_mesh(Z, transform)

    # 5 scaling / centring
    scale_mesh(mesh, scale=a.scale, scale_x=a.scale_x, scale_y=a.scale_y, scale_z=a.scale_z)
    if a.center:
        center_mesh(mesh)

        # 6 export
    mesh.export(a.output)
    print(f"Saved {a.output} (verts={len(mesh.vertices)})")


if __name__ == "__main__":
    sys.exit(main())
