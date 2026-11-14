# Data

All files in `raw/` are small subsets of public datasets. Values are unchanged.

## GravIS Level-3 TWS (GFZ)

- `GRAVIS-3_GFZOP_0600_TWS_europe_subset.nc` – GFZ GravIS RL06, version 0006,
  doi:10.5880/GFZ.GRAVIS_06_L3_TWS (CC BY 4.0)
- `GRAVIS-3_COSTG_0200_TWS_europe_subset.nc` – COST-G GravIS RL02, version 0001,
  doi:10.5880/COST-G.GRAVIS_02_L3_TWS (CC BY 4.0)

Original global files (~510 MB each) from https://isdc-data.gfz.de/grace/GravIS/ (downloaded
14 Nov 2026). They were cut to 12°W–32°E, 34–62°N with `subset_gravis.py`; longitudes were
converted from 0–360 to −180–180.

## ERA5-Land monthly means (Copernicus C3S)

- `era5_land_monthly_alps_2002_2026.nc` – snow depth water equivalent (`sd`, m w.e.) and
  volumetric soil water layers 1–4 (`swvl1`–`swvl4`), monthly averaged reanalysis,
  2002-01 to 2026-10, area 43–49°N, 4–18°E, 0.1° grid.
- Copernicus Climate Change Service (2019): ERA5-Land monthly averaged data from 1950 to
  present. C3S Climate Data Store, doi:10.24381/cds.68d2bb30. Model description: Muñoz-Sabater
  et al. (2021), ESSD 13, 4349–4383.
- Request reproduced in `download_era5_land.py` (needs a free CDS account and API key).

## Natural Earth

`raw/natural_earth/` – 1:10m countries and lakes (public domain), clipped to Europe.
