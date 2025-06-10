import rasterio, numpy as np, sys

tif = sys.argv[1]
with rasterio.open(tif) as src:
    Z = src.read(1)
    nodata = src.nodata
    water = Z <= 0

    # Dilate once so we get shoreline neighbours
    from scipy.ndimage import binary_dilation
    edge = binary_dilation(water) & ~water
    print("Pixels classified:")
    print(f"-water     {water.sum():,}")
    print(f"-shoreline {edge.sum():,}")

    # print stats
    for name, mask in [("water", water), ("shoreline", edge)]:
        vals = Z[mask]
        print(f"\n{name}: min {vals.min():.3f}, max {vals.max():.3f}, "
              f"mean {vals.mean():.3f}")
