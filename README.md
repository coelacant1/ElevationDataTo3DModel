# ElevationDataTo3DModel

A suite of Python tools for converting raw elevation data (point clouds, GeoTIFFs, etc.) into STL meshes, with optional visualization and mesh simplification steps.

![Point cloud view of Lake Huron](https://github.com/coelacant1/ElevationDataTo3DModel/blob/main/Pictures/erie_lld_logscale.jpg)



---

WIP Changes Still in Progress
- GeoTIFF Issues: Correcting Shoreline, abrupt change from land to water
- Point Cloud Viewer: Require testing for SpaceMouse input

---

## Table of Contents

* [Features](#features)
* [Getting Started](#getting-started)
  * [Prerequisites](#prerequisites)
  * [Installation](#installation)
* [Usage](#usage)
  * [1. Visualize or Export Point Cloud](#1-visualize-or-export-point-cloud)
  * [2. Simplify Point Cloud](#2-simplify-point-cloud)
  * [3. Convert Point Cloud to STL](#3-convert-point-cloud-to-stl)
  * [4. Convert GeoTIFF to STL](#4-convert-geotiff-to-stl)
* [Data Sources](#data-sources)
* [Example Outputs](#example-outputs)
* [Directory Structure](#directory-structure)
* [License](#license)

---

## Features

* **Interactive visualization** of large elevation point clouds
* **Log-scale elevation** rendering with adjustable midpoint, factor, and vertical scaling
* **Mesh simplification** to control resolution and final physical size
* **Point-cloud -> STL conversion**
* **GeoTIFF -> STL conversion**, with options for shoreline alignment, water-level slicing, decimation, scaling, and centering

---

## Getting Started

### Prerequisites

* Python 3.8+
* Pip Requirements
* Git

### Installation

1. Clone this repository and enter the project folder:

   ```bash
   git clone https://github.com/coelacant1/ElevationDataTo3DModel.git
   cd ElevationDataTo3DModel
   ```

2. Create and activate a virtual environment:

   ```bash
   python3 -m venv venv
   source venv/bin/activate
   ```

3. Install the required Python packages:

   ```bash
   python3 -m pip install -r requirements.txt
   ```

---

## Usage

### 1. Visualize or Export Point Cloud

```bash
# Basic visualization
python3 Programs/ElevationPointCloudViewer.py Data/huron_lld.xyz

# Log-scale elevation, custom factors, and vertical scaling
python3 Programs/ElevationPointCloudViewer.py \
  Data/huron_lld.xyz \
  --log-scale \
  --log-factor 1.0 \
  --log-midpoint 20.0 \
  --elevation-scale 100.0

# Visualization + export to file
python3 Programs/ElevationPointCloudViewer.py \
  Data/huron_lld.xyz \
  --log-scale \
  --log-factor 1.0 \
  --log-midpoint 20.0 \
  --elevation-scale 100.0 \
  --export
```

![Point cloud view of Lake Erie with Log Scale](https://github.com/coelacant1/ElevationDataTo3DModel/blob/main/Pictures/erie_lld_logscale.jpg)

![Point cloud view of Lake Erie with No Scaling](https://github.com/coelacant1/ElevationDataTo3DModel/blob/main/Pictures/erie_lld_noscale.jpg)

![Exported point cloud of Lake Erie in Blender](https://github.com/coelacant1/ElevationDataTo3DModel/blob/main/Pictures/erie_lld_export.jpg)

### 2. Simplify Point Cloud

```bash
# Simplify to target resolution factor (no resizing)
python3 Programs/GenerateSimplifiedData.py \
  huron_lld_es100.0_ls0_lm20.0_mf0.25_tf0.1.xyz

# Simplify + scale to final physical dimensions (width x height x z_range)
python3 Programs/GenerateSimplifiedData.py \
  huron_lld_es100.0_ls0_lm20.0_mf0.25_tf0.1.xyz \
  --width 2.0 \
  --height 1.0 \
  --z_range 0.2 \
  --resolution 0.1
```

### 3. Convert Point Cloud to STL

```bash
python3 Programs/ConvertDataToSTL.py simplified_data.xyz \
  -o output.stl
```

### 4. Convert GeoTIFF to STL

```bash
# Basic GeoTIFF → STL
python3 Programs/GeoTIFFToSTL.py \
  Data/DEMGlobalMosaic.tiff \
  Exports/output5.stl

# With shoreline alignment & custom scaling
python3 Programs/GeoTIFFToSTL.py \
  Data/DEMGlobalMosaic.tiff \
  Exports/output5.stl \
  --align-shoreline \
  --scale-x 100 \
  --scale-y 100 \
  --scale-z 0.02 \
  --center

# Add water‐level slicing, decimation, and debug output
python3 Programs/GeoTIFFToSTL.py \
  Data/DEMGlobalMosaic.tiff \
  Exports/output5.stl \
  --water-level 176 \
  --align-shoreline \
  --decimate 16 \
  --scale-x 100 \
  --scale-y 100 \
  --scale-z 0.02 \
  --center \
  --debug
```

![Contrast Auto-Scaled GeoTIFF Input](https://github.com/coelacant1/ElevationDataTo3DModel/blob/main/Pictures/GeoTIFFInput.png)

![Stress Display GeoTIFF Input](https://github.com/coelacant1/ElevationDataTo3DModel/blob/main/Pictures/GeoTIFFInput2.png)

![GeoTIFF Output Angled in Blender](https://github.com/coelacant1/ElevationDataTo3DModel/blob/main/Pictures/GeoTIFFOutput.png)

![GeoTIFF Output Top-Down in Blender](https://github.com/coelacant1/ElevationDataTo3DModel/blob/main/Pictures/GeoTIFFOutput2.png)

---

## Data Sources

* **Great Lakes Bathymetry XYZ Data** (NOAA)
  [https://www.ncei.noaa.gov/products/great-lakes-bathymetry](https://www.ncei.noaa.gov/products/great-lakes-bathymetry)

* **DEM Global Mosaic** (IHO/NOAA)
  Test Coordinates: `-85.17981, 40.41677, -77.51685, 44.22467`
  [https://www.ncei.noaa.gov/maps/iho\_dcdb/](https://www.ncei.noaa.gov/maps/iho_dcdb/)

---

## License

This project is licensed under the GNU General Public V3 License. See [LICENSE](LICENSE) for details.
