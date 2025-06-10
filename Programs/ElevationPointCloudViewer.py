import pyvista as pv
import numpy as np
import argparse
import vtk
import os

#elev_m = log(1 + (log_midpoint/shifted) × log_factor) × elevation_scale

def load_data(filename):
    # Use delimiter=None to split on any whitespace
    return np.loadtxt(filename, delimiter=None)

def _derive_scales(mid_frac: float, top_frac: float) -> tuple[float, float]:
    """
    Compute (above_scale, below_scale) so that:

        new_top_height   = top_frac  * H
        new_bottom_height = (1-top_frac) * H
        ───────────────────────────────────
        ⇒ above_scale = top_frac / (1-mid_frac)
           below_scale = (1-top_frac) / mid_frac
    """
    if not (0.0 < mid_frac < 1.0):
        raise ValueError("--mid-frac must be between 0 and 1")
    if not (0.0 < top_frac < 1.0):
        raise ValueError("--top-frac must be between 0 and 1")

    above_scale = top_frac / (1.0 - mid_frac)
    below_scale = (1.0 - top_frac) / mid_frac
    return above_scale, below_scale

def apply_midpoint_scaling_auto(z: np.ndarray, mid_frac: float | None, top_frac: float | None) -> np.ndarray:
    """
    Compress / stretch elevations so that:
      • everything above the midpoint (fraction mid_frac of the height)
        ends up occupying `top_frac` of the *final* height.
      • the final min–max span is identical to the original.
    Pass None for either argument to skip this step.
    """
    if mid_frac is None or top_frac is None:
        return z # nothing to do

    z_min, z_max = float(z.min()), float(z.max())
    H           = z_max - z_min
    z_mid       = z_min + mid_frac * H

    above_scale, below_scale = _derive_scales(mid_frac, top_frac)

    # vectorised transform
    above = z > z_mid
    below = ~above
    z_new = z.copy()
    z_new[above] = z_mid + (z[above] - z_mid) * above_scale
    z_new[below] = z_mid + (z[below] - z_mid) * below_scale
    return z_new

def apply_log_scale(
    elevations: np.ndarray,
    elevation_scale: float,
    log_factor: float = 1.0,
    log_midpoint: float = 1.0
) -> np.ndarray:
    """
    Return a log-scaled version of `elevations` in meters.

        z' = log1p( (e - e_min) / log_midpoint * log_factor ) * elevation_scale

    • Shifts the minimum to zero so the log() is well-behaved.
    • Uses np.log1p for numerical stability with small values.
    """
    log_midpoint = max(log_midpoint, 1e-12) # prevent /0
    shifted = np.clip(elevations - elevations.min(), 0, None)

    if np.all(shifted == 0): # flat terrain edge-case
        return np.zeros_like(shifted)

    return np.log1p((shifted / log_midpoint) * log_factor) * elevation_scale

def geographic_scale(
    data: np.ndarray,
    elevation_scale: float,
    use_log: bool = False,
    log_factor: float = 1.0,
    log_midpoint: float = 1.0,
    use_midpoint: bool = False,
    mid_frac: float = 0.25,
    top_frac: float = 0.10,
):
    """Convert lon/lat (deg) -> meters and return XYZ points + Z array."""
    longitudes, latitudes = data[:, 0], data[:, 1]
    elevations = np.nan_to_num(data[:, 2], nan=0.0)

    # degrees -> meters 
    avg_lat_rad = np.radians(np.mean(latitudes))
    lat_m = latitudes * 111_132
    lon_m = longitudes * (111_320 * np.cos(avg_lat_rad))

    # pick vertical transform
    if use_log:
        elev_m = apply_log_scale(
            elevations,
            elevation_scale=elevation_scale,
            log_factor=log_factor,
            log_midpoint=log_midpoint
        )
    else:
        elev_m = elevations * elevation_scale

    if use_midpoint:
        elev_m = apply_midpoint_scaling_auto(elev_m, mid_frac, top_frac)

    return np.column_stack((lon_m, lat_m, elev_m)), elev_m

def export_transformed_data(points: np.ndarray, src_file: str, args):
    base = os.path.splitext(os.path.basename(src_file))[0]
    export_name = (f"{base}_es{args.elevation_scale}"
                   f"_ls{int(args.use_log)}"
                   f"_lm{args.log_midpoint}"
                   f"_mf{args.mid_frac}_tf{args.top_frac}.xyz")
    np.savetxt(export_name, points, fmt="%.6f")
    print(f"Exported: {export_name}")

def visualize_point_cloud(points, elevations, log_scale=False):
    cloud = pv.PolyData(points)
    cloud['depth'] = elevations

    plotter = pv.Plotter()
    plotter.add_points(
        cloud,
        scalars='depth',
        cmap='viridis',
        point_size=2,
        render_points_as_spheres=True
    )
    plotter.background_color = 'white'
    plotter.add_scalar_bar(title='Depth (transformed m)', vertical=True)

    # Enable 3Dconnexion SpaceMouse support via VTK.
    interactor = plotter.iren.interactor
    style = vtk.vtkInteractorStyle3D()
    interactor.SetInteractorStyle(style)
    interactor.SetUseTDx(True)

    title = "Bathymetric Viewer (Log Scale)" if log_scale else "Bathymetric Viewer"
    plotter.show(title=title)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Geo-accurate bathymetry viewer with optional log-scale transformation to emphasize underwater detail and export the transformed data as an STL mesh.'
    )
    parser.add_argument("filenames", nargs="+", type=str,
                    help="One or more XYZ data files (space-separated)")
    parser.add_argument('--elevation-scale', type=float, default=25,
                        help='Vertical exaggeration factor (default: 25)')
    parser.add_argument('--log-scale', action='store_true',
                        help='Apply logarithmic transformation to elevation values.')
    parser.add_argument('--log-factor', type=float, default=1.0,
                        help='Factor to adjust the log scaling transformation (default: 1.0)')
    parser.add_argument('--log-midpoint', type=float, default=1.0,
                        help='Midpoint adjustment factor for log scaling (default: 1.0)')
    parser.add_argument('--export', action='store_true',
                        help='Export the transformed data as an STL file only (no visualization)')
    parser.add_argument("--use-log",     action="store_true")
    parser.add_argument("--use-midpoint", action="store_true")
    parser.add_argument("--mid-frac",    type=float, default=0.25)
    parser.add_argument("--top-frac",    type=float, default=0.10)

    args = parser.parse_args()

    all_points, all_elev = [], []

    for fname in args.filenames:
        data = np.loadtxt(fname)

        points, elev_m = geographic_scale(
            data,
            args.elevation_scale,
            use_log=args.use_log,
            log_factor=args.log_factor,
            log_midpoint=args.log_midpoint,
            use_midpoint=args.use_midpoint,
            mid_frac=args.mid_frac,
            top_frac=args.top_frac,
        )

        if args.export:
            export_transformed_data(points, fname, args)
        else:
            all_points.append(points)
            all_elev.append(elev_m)

    if not args.export:
        points_combined = np.vstack(all_points)
        elev_combined   = np.concatenate(all_elev)
        visualize_point_cloud(points_combined, elev_combined, log_scale=args.use_log)
