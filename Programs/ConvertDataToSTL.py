#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse
from scipy.spatial import Delaunay
from stl import mesh

def generate_stl(input_file, output_file):
    # Read the point cloud file; it should be space-delimited with headers:
    # x_scaled, y_scaled, z_scaled
    df = pd.read_csv(input_file, sep='\s+')
    
    # Ensure that the required columns exist.
    try:
        points = df[['x_scaled', 'y_scaled', 'z_scaled']].values
    except KeyError:
        raise ValueError("Input file must contain columns 'x_scaled', 'y_scaled', and 'z_scaled'.")

    # Perform 2D Delaunay triangulation using the x and y coordinates.
    tri = Delaunay(points[:, :2])
    
    # Create an array of triangles (n_triangles, 3, 3): each triangle has 3 vertices in 3D.
    triangles = points[tri.simplices]

    # Create the STL mesh using numpy-stl.
    stl_mesh = mesh.Mesh(np.zeros(triangles.shape[0], dtype=mesh.Mesh.dtype))
    for i, triangle in enumerate(triangles):
        stl_mesh.vectors[i] = triangle
        # The normal vectors will be computed automatically by numpy-stl.

    # Save the mesh to an STL file.
    stl_mesh.save(output_file)
    print(f"STL file saved to '{output_file}'.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Convert a point cloud (with x_scaled, y_scaled, z_scaled columns) to an STL file."
    )
    parser.add_argument("input_file", help="Input file name containing the point cloud data")
    parser.add_argument("-o", "--output", default="output.stl",
                        help="Output STL file name (default: output.stl)")
    
    args = parser.parse_args()
    generate_stl(args.input_file, args.output)
