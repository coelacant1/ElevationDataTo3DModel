#!/usr/bin/env python3
import pandas as pd
import argparse

def process_file(input_file, output_file, scaled_width, scaled_height, z_range, resolution):
    # Define column names (no header in the file)
    col_names = ['longitude', 'latitude', 'elevation']
    # Read the space-delimited file
    df = pd.read_csv(input_file, sep='\s+', header=None, names=col_names)
    
    # Calculate the original extents for X and Y
    x_min = df['longitude'].min()
    x_max = df['longitude'].max()
    y_min = df['latitude'].min()
    y_max = df['latitude'].max()
    
    # Scale the X and Y coordinates
    df['x_scaled'] = (df['longitude'] - x_min) / (x_max - x_min) * scaled_width
    df['y_scaled'] = (df['latitude'] - y_min) / (y_max - y_min) * scaled_height
    
    # Calculate the original extents for the Z coordinate
    z_min = df['elevation'].min()
    z_max = df['elevation'].max()
    
    # Scale the elevation to the desired z_range.
    df['z_scaled'] = (df['elevation'] - z_min) / (z_max - z_min) * z_range
    
    # Create grid cell indices for the scaled X and Y coordinates
    df['grid_x'] = (df['x_scaled'] / resolution).astype(int)
    df['grid_y'] = (df['y_scaled'] / resolution).astype(int)
    
    # Group by the grid cell and take the first point as the representative
    simplified_df = df.groupby(['grid_y', 'grid_x']).first().reset_index()
    
    # Save the simplified data to the output file (space-delimited)
    simplified_df.to_csv(output_file, sep=" ", index=False)
    print(f"Simplified data saved to '{output_file}'.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Scale and simplify geospatial data from a space-delimited file, including scaling the Z coordinate."
    )
    parser.add_argument("input_file", help="Input file name containing space-separated data")
    parser.add_argument("-o", "--output", default="simplified_data.xyz",
                        help="Output file name (default: simplified_data.xyz)")
    parser.add_argument("--width", type=float, default=3.5,
                        help="Desired output width in meters (default: 1000)")
    parser.add_argument("--height", type=float, default=1.5,
                        help="Desired output height in meters (default: 800)")
    parser.add_argument("--z_range", type=float, default=0.2,
                        help="Desired output range for elevation (z coordinate) in meters (default: 100)")
    parser.add_argument("--resolution", type=float, default=0.01,
                        help="Grid resolution in meters (default: 1.0)")
    
    args = parser.parse_args()
    
    process_file(args.input_file, args.output, args.width, args.height, args.z_range, args.resolution)
