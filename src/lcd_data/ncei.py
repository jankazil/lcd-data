'''
Tools for download of Local Climatological Data (LCD) version 2
(https://www.ncei.noaa.gov/products/land-based-station/local-climatological-data)
from National Centers for Environmental Information (NCEI).
'''

import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from pathlib import Path

import requests

from lcd_data import web

# Notes:
#
# - the GHCNh station list file provides more metadata than the LCD
#   station list file even for stations in the US,
# - the LCD station list file is not limited to US stations, but gives
#   nearly the same stations as the GHCNh station list file
#
# Therefore, the GHCNh station list is better and is used in lcd-data.
# This means that one needs to determine later if LCD data are actually
# available for the US stations in the GHCNh station list.

# GHCNh station meta information file URLs

ghcnh_url = 'https://www.ncei.noaa.gov/oa/global-historical-climatology-network'

ghcnh_station_list_url = ghcnh_url + '/hourly/doc/ghcnh-station-list.txt'
ghcnh_station_doc_url = ghcnh_url + '/hourly/doc/ghcnh_DOCUMENTATION.pdf'

# LCD URLs

lcd_url = 'https://www.ncei.noaa.gov/oa/local-climatological-data'

lcd_station_list_url = lcd_url + '/v2/doc/lcdv2-station-list.txt'
lcd_station_doc_url = lcd_url + '/v2/doc/lcdv2_DOCUMENTATION.pdf'


def download_stations_meta_files(local_dir: Path):
    '''
    Downloads the LCD files with station meta information

    Args:
        local_dir (Path): Local directory where the downloaded files will be saved.
    '''

    # Create parent directory if needed

    local_dir.mkdir(parents=True, exist_ok=True)

    # Download files - we download both the GHCNh and the LCD files.

    web.download_file(ghcnh_station_list_url, local_dir, verbose=True)
    web.download_file(ghcnh_station_doc_url, local_dir, verbose=True)
    web.download_file(lcd_station_list_url, local_dir, verbose=True)
    web.download_file(lcd_station_doc_url, local_dir, verbose=True)

    return


def lcd_data_file_name(year: int, station_id: str) -> str:
    """
    Constructs the name of a NCEI LCD data file.

    Args:
        year (int): Gregorian year of the data
        station_id (str): Station ID (GHCNh station identifier)

    Returns:
        str: Local file name.
    """

    file_name = 'LCD_' + station_id + '_' + str(year) + '.csv'

    return file_name


def lcd_data_file_paths(
    start_year: int,
    end_year: int,
    station_ids: list[str],
    local_dir: Path,
) -> list[Path]:
    """
    For a given year range (inclusive) and given station IDs (GHCNh station identifiers),
    and a local directory path, returns the paths to LCD data files in that directory path.

    Args:
        start_year (int): Gregorian year of the first data file to be downloaded
        end_year (int): Gregorian year of the last data file to be downloaded
        station_ids (list[str]): A list containing LCD station IDs (GHCNh station identifiers)
        local_dir (Path): Local directory where the downloaded files will be saved.

    Returns:
        list[Path]: List of local paths of the downloaded files.
    """

    # Construct local file paths

    all_local_file_paths = []

    for year in range(start_year, end_year + 1):
        local_file_paths = []

        for station_id in station_ids:
            local_file_paths.append(local_dir / lcd_data_file_name(year, station_id))

        all_local_file_paths = all_local_file_paths + local_file_paths

    return all_local_file_paths


def lcd_data_url(year: int, station_id: str) -> str:
    """
    Constructs the URL of a NCEI LCD station data file.

    Args:
        year (int): Gregorian year of the data
        station_id (str): Station ID (GHCNh station identifier)

    Returns:
        str: URL of a NCEI LCD station data file
    """

    url = lcd_url.rstrip('/') + '/v2/access/' + str(year) + '/' + lcd_data_file_name(year, station_id)

    return url


def lcd_data_urls(station_ids: list[str], start_year: int, end_year: int, n_jobs: int = 16, verbose: bool = False) -> list[str]:
    """
    Collects and returns URLs of all available LCD data file from the NCEI server
    for a list of stations and for an inclusive range of years.

    Args:
        station_ids (list[str]) : List of station IDs (GHCNh station identifiers).
        start_year (int): First year to include (inclusive).
        end_year   (int): Last year to include (inclusive). Must be >= start_year.
        n_jobs (int): Maximum number of parallel web access operations
        verbose (bool): If True, print information. Defaults to False.

    Returns:
        list[str]: Absolute URLs of files under the specified year directories.
    """

    print()

    urls = []

    for year in range(start_year, end_year + 1):
        year_dir = lcd_url + '/index.html#v2/access/' + str(year)

        if verbose:
            print('Collecting data file URLs from NCEI server directory', year_dir, flush=True)

        candidates = [f"{lcd_url}/v2/access/{year}/LCD_{sid}_{year}.csv" for sid in station_ids]

        with ThreadPoolExecutor(max_workers=n_jobs) as exe:
            for url, ok in zip(candidates, exe.map(web.head_ok, candidates), strict=False):
                if ok:
                    urls.append(url)

    return urls


def download_many(
    start_year: int,
    end_year: int,
    station_ids: list[str],
    local_dir: Path,
    n_jobs: int = 1,
    refresh: bool = False,
    verbose: bool = False,
) -> list[Path]:
    """
    Downloads LCD data files for a given year range (inclusive) and given station IDs
    (GHCNh station identifiers), even for files that already exists locally, if requested.
    A given number of parallel threads is used to accelerate download. The routine is
    parallelized over stations, but not over the year range.

    Why an inclusive listing? Because that is the canonical definition of a range of years -
    it is not meant to exclude the last one.

    Args:
        start_year (int): Gregorian year of the first data file to be downloaded
        end_year (int): Gregorian year of the last data file to be downloaded
        station_ids (list[str]): A list containing LCD station IDs (GHCNh station identifiers)
        local_dir (Path): Local directory where the downloaded files will be saved.
        n_jobs (int): Maximum number of parallel downloads
        refresh (bool, optional): If True, download even if the file already exists. Defaults to False.
        verbose (bool): If True, print information. Defaults to False.

    Returns:
        list[Path]: List of local paths of the downloaded files.
    """

    if verbose:
        print()
        print('Downloading observations from NCEI server')
        print()

    # Construct URLs and local file paths

    all_local_file_paths = []

    for year in range(start_year, end_year + 1):
        urls = []

        local_file_paths = []

        for station_id in station_ids:
            urls.append(lcd_data_url(year, station_id))

        local_file_paths = web.download_threaded(urls, local_dir, n_jobs=n_jobs, refresh=refresh, verbose=verbose)

        all_local_file_paths = all_local_file_paths + local_file_paths

    return all_local_file_paths


def get_period_of_record(station_id: str, token: str, timeout: int = 60, max_retries: int = 1200, retry_delay: int = 3) -> dict:
    '''
    Retrieve the Period of Record (minimum and maximum available dates) for a
    given LCD (Local Climatological Data) station from NOAA's CDO (Climate Data
    Online) v2 API, with retry logic.

    Args:
        station_id (str): Station ID (GHCNh station identifier).
        token (str): A valid NOAA CDO API token for authentication.
        timeout (int): (Optional) Seconds to wait for a response before raising
                       a requests.exceptions.Timeout error.
        max_retries (int): (Optional) Number of times to retry the request
                           after a failure.
        retry_delay (int): (Optional) Number of seconds to wait between retries.

    Returns:
        dict: A dictionary containing:
            - station_id (str)
            - name (str or None)
            - mindate (datetime or None)
            - maxdate (datetime or None)
            - source_url (str)

    Raises:
        requests.exceptions.HTTPError: If the API request fails after all retries.
        requests.exceptions.RequestException: For other network-related errors.
    '''
    wban_id = station_id[-5:]
    station_id = f'WBAN:{wban_id}'

    url = f'https://www.ncei.noaa.gov/cdo-web/api/v2/stations/{station_id}'
    headers = {'token': token}
    params = {'datasetid': 'LCD'}

    last_exception = None

    for attempt in range(1, max_retries + 1):
        try:
            r = requests.get(url, headers=headers, params=params, timeout=timeout)
            r.raise_for_status()
            rec = r.json()

            mindate = rec.get('mindate')
            if mindate is not None:
                mindate = datetime.strptime(mindate, '%Y-%m-%d')

            maxdate = rec.get('maxdate')
            if maxdate is not None:
                maxdate = datetime.strptime(maxdate, '%Y-%m-%d')
                maxdate = maxdate.replace(hour=23, minute=59, second=59)

            return {
                'station_id': station_id,
                'name': rec.get('name'),
                'mindate': mindate,
                'maxdate': maxdate,
                'source_url': r.url,
            }

        except (requests.exceptions.RequestException, ValueError) as e:
            last_exception = e
            if attempt < max_retries:
                time.sleep(retry_delay)
            else:
                raise last_exception from None
