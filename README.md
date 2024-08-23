# arcpyGeoUtilsX
# arcpyGeoUtilsX

`arcpyGeoUtilsX` is a Python package designed for advanced GIS processing using ArcGIS's `arcpy` module. This package provides utilities for tasks like DEM processing, watershed delineation, and runoff analysis.

## Installation

**Note:** This package requires ArcGIS Pro with the `arcpy` module.

## Usage

Here's a simple example of how to use the package:

```python
from arcpyGeoUtilsX import SmallProjectsGISprocessing

processor = SmallProjectsGISprocessing('path_to_your_file')
clipped_raster = processor.extrdem()
