# lcd-data

**lcd-data** is a Python toolkit for downloading and processing [Local Climatological Data version 2 (LCDv2) ](https://www.ncei.noaa.gov/products/land-based-station/local-climatological-data) data.

It provides:

- A top level command-line tool that
  
  - automates the download of LCD v2 station observations for
    - individual stations
    - U.S. states and territories
    - Regional Transmission Organization (RTO) / Independent System Operator (ISO) regions
    
  - constructs full-hourly UTC time series of
    - temperature at 2 m
    - dew point temperature at 2 m
    - relative humidity at 2 m
    - wind speed at 10 m
  
    from the irregularly spaced, local time LCD v2 station observation time series, for a selected station or for stations in the selected US state/territory or RTO/ISO region, and a user-specified time range. The time series are saved in a netCDF file.  
    
- Modules for processing LCD v2 station observations.

LCD v2 is provided by the [National Centers for Environmental Information (NCEI)](https://www.ncei.noaa.gov/).

## Installation (Linux / macOS)

```bash
pip install git+https://github.com/jankazil/lcd-data
```

## Overview

The package provides a command-line tool that selects stations by geography (a single station by GHCNh identifier, a U.S. state or territory, RTO/ISO regions, and the special region CONUS representing the contiguous U.S.), checks data availability, downloads LCD v2 observation files for a given year range, constructs full-hourly UTC time series for the observables, and saves them in a NetCDF file. It optionally generate plots showing the original and the interpolated time series.

Geospatial region selection is based on U.S. Energy Information Administration definitions of RTO/ISO footprints, and U.S. Census Bureau state/territory boundaries, included with the package.

The list of GHCNh station identifiers is available [here](https://www.ncei.noaa.gov/oa/global-historical-climatology-network/hourly/doc/ghcnh-station-list.txt). LCD v2 contains only U.S. stations.

## Workflow

The following describes the internal workflow performed by the command-line tool:

1. Load the region geometry (RTO/ISO polygons or U.S. state/territory boundaries) if a region is specified; skip this step if a station ID is provided.
2. Retrieve the station list from NCEI and either filter it spatially by region or select the specified station.
3. Filter the stations by data availability for the requested year range, either online by probing NCEI or offline by checking local files.
4. Save the filtered station list for reference.
5. Download LCD v2 observation files from NCEI for the selected stations and years, skipping files already present that match by ETag.
6. Create full-hourly UTC time series for temperature (T), dew point temperature (Td), relative humidity (RH), and wind speed by converting local observation time to UTC and interpolating the data to full hours. Remove temperatures above 60 °C. Perform interpolation only across gaps of up to 2 hours. Derive RH from T and Td.
7. Optionally create comparison plots for the original and interpolated series.
8. Save the full-hourly UTC time series in a NetCDF file, for the given station or the stations in the state/region.

**Notes:** Interpolation of station observation time series across many years and/or many stations can be slow due to inherent limitations of Python. Creating plots is very slow and recommended only for individual stations (as opposed to regions).

## Command-line interface (CLI)

The CLI is exposed as `"build-lcd-dataset"` when installed:

```bash
# Provide usage information, a list of two-letter U.S. state/territory codes, and a list of RTO/ISO region names:
build-lcd-dataset --help
  
build-lcd-dataset 2020 2025 USW00003017 /path/to/data -p /path/to/plots
  
build-lcd-dataset 2022 2022 ERCOT /path/to/data -n 32
  
build-lcd-dataset 2021 2021 CO /path/to/data --offline
```

**Positional arguments**  

- `start_year` and `end_year`: Inclusive range of years.  
- `region_name`: Two-letter U.S. state or territory code (e.g., `CA`, `PR`), the special region `CONUS`, one of the RTO/ISO codes (`ERCOT`, `CAISO`, `ISONE`, `NYISO`, `MISO`, `SPP`, `PJM`), or a U.S. station identifier as in [the GHCNh station list](https://www.ncei.noaa.gov/oa/global-historical-climatology-network/hourly/doc/ghcnh-station-list.txt).  
- `data_dir`: Destination directory for station lists, downloads, and outputs.

**Options**  

- `-n, --n INT`: Maximum number of parallel downloads.  
- `-o, --offline`: Work offline; expect required files to be present in `data_dir`.  
- `-p, --plotdir PATH`: Directory to write comparison plots; plotting is slow.

## Sample results

Original and interpolated full-hourly UTC time series in November 2024, Twentynine Palms, CA:

![LCD station USL000ANVC1 time series, December 2024](plots/USW00093121.Nov-2024.png)  

## Public API

### Modules

#### `lcd_data.ncei`
Utilities for station metadata and LCD v2 downloads.

- `download_stations_meta_files(local_dir)`: Download GHCNh and LCD v2 station meta documents.
- `lcd_data_file_name(year, station_id)`: Construct the canonical LCD v2 observation file name.
- `lcd_data_file_paths(start_year, end_year, station_ids, local_dir)`: Build local paths for all expected files.
- `lcd_data_url(year, station_id)`: Build the absolute URL to a LCD v2 observation file.
- `lcd_data_urls(station_ids, start_year, end_year, n_jobs)`: Probe NCEI server to list existing files.
- `download_many(...)` and `download_threaded(...)`: Concurrent file downloads with optional refresh behavior.
- `download_file(url, local_dir, refresh=False, verbose=False)`: Robust download with ETag checking and retries.

#### `lcd_data.rto_iso`
Helpers to work with RTO/ISO region polygons.

- `REGION_NAMES`: `['CAISO', 'ERCOT', 'ISONE', 'NYISO', 'MISO', 'PJM', 'SPP']`.
- `regions(rto_iso_geojson)`: Read GeoJSON and return a GeoDataFrame with merged geometries for each region.
- `region(rto_iso_geojson, region_name)`: Return a GeoDataFrame for the requested region.

#### `lcd_data.saturation`
Saturation vapor pressure and relative humidity utilities.

- `esatw(T)`: Saturation vapor pressure over liquid water (hPa) using an 8th‑order polynomial fit.
- `rh(T, Td)`: Relative humidity (%) computed from temperature and dew point.

#### `lcd_data.region_codes`
Provides

  - `lcd_data.region_codes.countries`: Three-letter ISO 3166-1 alpha-3 country codes
  - `lcd_data.region_codes.us_states_territories`: Two-letter U.S. state or territory codes
  - `lcd_data.region_codes.conus`: The special `CONUS` region code
  - `lcd_data.region_codes.rto_iso_regions`: RTO/ISO region codes

#### `lcd_data.stations`
Station catalog handling, filtering, reading, interpolation, and writing.

- `Stations.from_url()` / `Stations.from_file(path)`: Build the station catalog from GHCNh-format metadata.
- Spatial selection by region geometry with `filter_by_region(region_gdf)` and by bounding box with `filter_by_coordinates(...)`.
- Availability filters: `filter_by_data_availability_online(start_time, end_time, n_jobs, verbose)` and `filter_by_data_availability_offline(data_dir, start_time, end_time, verbose)`.
- Station utilities: `filter_by_id(station_id)`, `ids()`, `save_station_list(path)`.
- `read_station_observations(...)`: Read and clean per‑station LCD v2 observation files; convert times to UTC, coerce numeric columns, correct Celsius-with-18.3° base fields, drop non-observational report types, limit unrealistic temperatures, and compute hourly RH.
- `construct_hourly(...)`: Build full-hourly UTC series for `T`, `Td`, `RH`, and `windspeed`, with optional plotting and gap-limited interpolation.
- `write_utc_hourly_netcdf(path)`: Save the hourly dataset to NetCDF with safe encodings.

## Development

### Code Quality and Testing Commands

- `make fmt` - Runs ruff format, which automatically reformats Python files according to the style rules in `pyproject.toml`
- `make lint` - Runs ruff check - -fix, which lints the code (checks for style errors, bugs, outdated patterns, etc.) and auto-fixes what it can.
- `make check` - Runs fmt and lint.
- `make type` - Currently disabled. Runs mypy, the static type checker, using the strictness settings from `pyproject.toml`. Mypy is a static type checker for Python, a dynamically typed language. Because static analysis cannot account for all dynamic runtime behaviors, mypy may report false positives which do no reflect actual runtime issues.
- `make test` - Runs pytest with reporting (configured in `pyproject.toml`).

## Disclaimers

The LCD v2 data accessed by this software are publicly available from NOAA's National Centers for Environmental Information (NCEI) and are subject to their terms of use. This project is not affiliated with or endorsed by NOAA.

This software uses U.S. Census Bureau and U.S. Energy Information Administration data, but is neither endorsed nor certified by the U.S. Census Bureau or the U.S. Energy Information Administration.

## Author

Jan Kazil - jan.kazil.dev@gmail.com - [jankazil.com](https://jankazil.com)

## License

BSD-3-Clause
