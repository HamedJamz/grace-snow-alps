"""
Download the ERA5-Land monthly data used in this project from the Copernicus Climate Data Store.

Needs a free CDS account, the ERA5-Land licence accepted on the dataset page, and an API key
in ~/.cdsapirc (see https://cds.climate.copernicus.eu/how-to-api). Then:

    pip install cdsapi
    python data/download_era5_land.py
"""

import cdsapi

request = {
    "product_type": ["monthly_averaged_reanalysis"],
    "variable": [
        "snow_depth_water_equivalent",
        "volumetric_soil_water_layer_1",
        "volumetric_soil_water_layer_2",
        "volumetric_soil_water_layer_3",
        "volumetric_soil_water_layer_4",
    ],
    "year": [str(y) for y in range(2002, 2027)],
    "month": [f"{m:02d}" for m in range(1, 13)],
    "time": ["00:00"],
    "data_format": "netcdf",
    "download_format": "unarchived",
    "area": [49, 4, 43, 18],          # North, West, South, East
}

cdsapi.Client().retrieve("reanalysis-era5-land-monthly-means", request).download(
    "data/raw/era5_land_monthly_alps_2002_2026.nc")
