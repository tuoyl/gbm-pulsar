#!/usr/bin/env python
# coding: utf-8

import numpy as np
from gbmpulsar.GBM_data_filter import UTC_TO_MET, filter
from gbmpulsar.barycor.barycor import barycor
import argparse
import datetime
import glob
import os
import wget
from tqdm import tqdm
from astropy.io import fits
from astropy.table import Table

"""
A pipeline tool to do the following procesures:
    1. Download the GBM data (if required)
    2. filter the data for the given time range (UTC)
    3. save the data to FITS file (store the data daily) with PHA and detector number (if required)
    4. barycentric correction

The input paramters should be:
    1. Time range of interests in UTC (e.g., '2022-07-27T00:00:00')
    2. GBM data directory (keep the hierarchy of the GBM file directory), if data exists then skip the downloading
    3. output directory to save the filtered data
    4. stem output
    5. RA, DEC in degree (J2000)
    optionals:
    1. Energy (keV) range to filtering (default is to keep all PHA)
    2. whether to store the PHA information (default is No for the sake of storage saving)
    3. whether to store the Detector information (default is No for the sake of storage saving)
    4. whether to Download data
    5. whether to execute Barycentric correction
"""


NOTICE = """

/---------------------------------------------\\
|   GBM data analysis tool for X-ray pulsar   |
\---------------------------------------------/

Usage:

        python GBM_pulsar_pipeline.py --gbm_dir="/path/to/GRM/data" \\
                --tstart='2022-07-27T00:00:00' --tstop='2022-07-28T00:00:00' \\
                --output_dir="/path/to/save/data" \\
                --stem="gbmCrab" --ra=83.883225 --dec=22.014458333333334 \\
                --barycor --jplephem="./barycor/de421.bsp"

    some optional:

        python GBM_pulsar_pipeline.py --gbm_dir="/path/to/GRM/data" \\
                --tstart='2022-07-27T00:00:00' --tstop='2022-07-28T00:00:00' \\
                --output_dir="/path/to/save/data" \\
                --stem="gbmCrab" --ra=83.883225 --dec=22.014458333333334 \\
                --barycor --jplephem="./barycor/de421.bsp" \\
                --energylow=8 --energyhigh=25 \\
                --store_pha --store_det --accelerate
    """

def parse_args():

    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=NOTICE)

    parser.add_argument(
            "--gbm_dir",
            type=str,
            required=True,
            help="GBM data directory (several folders named after years in this directory)")

    parser.add_argument(
            "--tstart",
            type=str,
            required=True,
            help="start time of time range (in UTC)")


    parser.add_argument(
            "--tstop",
            type=str,
            required=True,
            help="stop time of time range (in UTC)")

    parser.add_argument(
            "--output_dir",
            type=str,
            required=False,
            help="directory to save the product file")

    parser.add_argument(
            "--stem",
            type=str,
            required=True,
            help="STEM for FITS output file")

    parser.add_argument(
            "--ra",
            type=float,
            required=True,
            help="RA for the source (degree, J200)")

    parser.add_argument(
            "--dec",
            type=float,
            required=True,
            help="Declination for the source (degree, J200)")

    # optional
    parser.add_argument(
            "--download",
            action="store_true",
            help="Retrieve data from GBM FTP, if data does not exist")

    parser.add_argument(
            "--barycor",
            action="store_true",
            help="To carry out barycentric correction")

    parser.add_argument(
            "--jplephem",
            type=str,
            help="specify the JPLEPH file")

    parser.add_argument(
            "--energylow",
            type=float,
            required=False,
            help="lower limits to filter the photons")

    parser.add_argument(
            "--energyhigh",
            type=float,
            required=False,
            help="upper limits to filter the photons")

    parser.add_argument(
            "--store_pha",
            action="store_true",
            help="to store the PHA column to output file")

    parser.add_argument(
            "--store_det",
            action="store_true",
            help="to store the Detector column to output file")

    parser.add_argument(
            "--accelerate",
            action="store_true",
            required=False,
            default=False,
            help="to accelerate the barycentric correction processes by barycor sampling photons, and interpolate\
                    the barycored arrival times for whole time series")

    return parser.parse_args()

def get_one_day_files(gbm_data_dir, year, month, day, hour='all', direction='left'):
    """
    return all available data in specific day
    if hour='all', all hourly data in that day are used,

    Parameters
    ----------
    gbm_data_dir: string
        GBM data directory

    hour: string or int
        if 'all': retrieve all hourly data
        if is `int`: then the hour is the creteria to select date before/after
        the hour (include) in the day, use the parameter `direction` to select
        the data {'left':before the hour; 'right':after the hour}.

    direction: string
        'left' or 'right'

    Returns
    -------
    event_list: list
        the list of event file (absolute path) in specific day

    poshist: string
        the poshist file
        glg_poshist_all_220315_v00
    """

    year = str(year)
    month= str(month).zfill(2)
    day  = str(day).zfill(2)
    if hour == 'all':
        hour = '*'
    elif direction == 'left':
        hour = int(hour)
        hour = np.linspace(0, hour, hour+1, dtype=int)
        #hour = '[{}]'.format( ','.join([str(x) for x in hour]) )
        hour = '*'
        # FIXME: now whole day is used, hour does not working
    elif direction == 'right':
        hour = int(hour)
        hour = np.linspace(hour, 23, 24-hour, dtype=int)
        #hour = '[{}]'.format( ','.join([str(x) for x in hour]) )
        hour = '*'
        # FIXME: now whole day is used, hour does not working

    event_list = glob.glob(
            os.path.join(gbm_data_dir, year, month, day, 'current',
                "glg_tte_*_{}{}{}_{}z_v*.fit.gz".format(year[2:], month, day, hour)))
    poshist = glob.glob(
            os.path.join(gbm_data_dir, year, month, day, 'current',
                "glg_poshist_all_{}{}{}_v*.fit".format(year[2:], month, day)))
    if len(poshist) != 0:
        poshist = poshist[-1]

    return event_list, poshist

def retrieve_data(data_dir, year, month, day,
        ftp_link="https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily/",
        detectors = ['n0','n1','n2','n3','n4','n5','n6','n7','n8','n9','na','nb', 'b0', 'b1']):

    """
    Retrieve data for specfic day from GBM FTP
    e.g.:
    glg_tte_n0_220201_23z_v00.fit.gz
    glg_poshist_all_220301_v00.fit

    data_dir: string
        the directory to save ALL GRB data (parent dir)
    """

    year = str(year)
    month= str(month).zfill(2)
    day  = str(day).zfill(2)

    # Creating folder
    out_dir = os.path.join(data_dir, year, month, day, 'current')
    if not os.path.exists(out_dir):
        print(f"...Making directory {out_dir}...")
        os.system("mkdir -p {}".format(
            out_dir))

    # Download Orbital file
    version = [str(x).zfill(2) for x in np.arange(2)]
    for ver in version:
        print(f"try to download different version of data ")
        try:
            poshist  = f"glg_poshist_all_{str(year)[2:]}{month}{day}_v{ver}.fit"
            if os.path.exists(os.path.join(out_dir, poshist)):
                print(f"...{poshist} already exists...")
            else:
                url = os.path.join(ftp_link, year, month, day, 'current', poshist)
                print(f"Downloading {url}")
                wget.download(url, out_dir+'/')
        except:
            continue

    # Download NaI Data
    if (not 'b0' in detectors) & (not 'b1' in detectors):
        for ver in version:
            print(f"try to download different version of data ")
            try:
            # Download Events
                for det in detectors:
                    for hour in np.linspace(0, 23, 24, dtype=int):
                        hstr = str(hour).zfill(2)
                        evtfile = f"glg_tte_{det}_{str(year)[2:]}{month}{day}_{hstr}z_v{ver}.fit.gz"
                        if os.path.exists(os.path.join(out_dir, evtfile)):
                            continue
                        else:
                            url = os.path.join(ftp_link, year, month, day, 'current', evtfile)
                            print(f"Downloading {url}")
                            wget.download(url, out_dir)
            except:continue

    # Download NaI/BGO Data
    if ('b0' in detectors) or ('b1' in detectors):
        for ver in version:
            print(f"try to download different version of data ")
            try:
            # Download Events
                for det in ['b0', 'b1']:
                    for hour in np.linspace(0, 23, 24, dtype=int):
                        hstr = str(hour).zfill(2)
                        evtfile = f"glg_tte_{det}_{str(year)[2:]}{month}{day}_{hstr}z_v{ver}.fit.gz"
                        if os.path.exists(os.path.join(out_dir, evtfile)):
                            continue
                        else:
                            url = os.path.join(ftp_link, year, month, day, 'current', evtfile)
                            print(f"Downloading {url}")
                            wget.download(url, out_dir)
            except:continue

def _resolve_evtname(evtfile):
    """
    resolve the detector name and the met array from an
    input event file

    det_shortname:
        the short name in the filename, e.g., 'n1', 'na'

    det_headname:
        the name in the header, e.g., "NAI_11"
    """
    filename = os.path.basename(evtfile)
    hdulist = fits.open(evtfile)

    det_shortname = filename.split('_')[2]
#    if len(hdulist) < 4:
#        ## File Error
#        raise IOError("File {} extension error".format(filename)) 
    det_headname = hdulist["EVENTS"].header['DETNAM']
    met = hdulist["EVENTS"].data.field("TIME")
    pha = hdulist["EVENTS"].data.field("PHA")
    return det_shortname, det_headname, met, pha

def save_to_fits(outname, data, names):
    """
    save data to FITS file

    Paramters
    ---------
    data: list
        a list of arrays to store

    names: list
        a list of column names
    """
    t = Table(data, names=names)
    t.write(outname, overwrite=True)
    os.system(f"gzip -f {outname}")

def main():
    args = parse_args()

    tstart_date = datetime.datetime.strptime(args.tstart, "%Y-%m-%dT%H:%M:%S")
    tstop_date  = datetime.datetime.strptime(args.tstop, "%Y-%m-%dT%H:%M:%S")

    begin_day = datetime.date(tstart_date.year, tstart_date.month, tstart_date.day)
    end_day   = datetime.date(tstop_date.year, tstop_date.month, tstop_date.day)
    one_day = datetime.timedelta(days=1)
    next_day  = begin_day

    detectors = ['n0','n1','n2','n3','n4','n5','n6','n7','n8','n9','na','nb']
    ## Loop for each day between Start and Stop
    while next_day <= end_day:
        year = next_day.year
        month= next_day.month
        day  = next_day.day

        if args.download:
            retrieve_data(
                    args.gbm_dir,
                    year, month, day,
                    detectors=detectors)

        print(f"...Analyzing data in {next_day}...")

        if next_day == begin_day:
            # first or last day, deal with the hour
            evtfiles, poshist = get_one_day_files(args.gbm_dir,
                    year, month, day,
                    hour=args.tstart.split(':')[1],
                    direction='right')
        elif next_day == end_day:
            # first or last day, deal with the hour
            evtfiles, poshist = get_one_day_files(args.gbm_dir,
                    year, month, day,
                    hour=args.tstop.split(':')[1],
                    direction='left')
        else:
            evtfiles, poshist = get_one_day_files(args.gbm_dir,
                    year, month, day,
                    hour='all')

        if (evtfiles == []) or (poshist == ''):
            print(f">>>file does not exists...")
            next_day += one_day
            continue

        # Save the data hourly! Dayly will cost a lot of RAM
        for hour in np.linspace(0, 23, 24, dtype=int):
            met_one_hour = np.array([])
            pha_one_hour = np.array([])
            det_one_hour = []

            for evtfile in tqdm(sorted(evtfiles), desc=f"Analyzing data in {hour+1}/24 hour of one day"):
                if os.path.basename(evtfile).split('_')[4] != str(hour).zfill(2)+'z':
                    continue
                try:
                    det_shortname, det_headname, met, pha = _resolve_evtname(evtfile)
                except IOError:
                    print("file {} Extension error".format(evtfile))
                    continue
                met_filtered, mask = filter(poshist,
                        detector=det_shortname,
                        radec=np.array([args.ra, args.dec]),
                        met=met,
                        angle_incident=70,
                        retrieve_mask=True)

                # store
                met_one_hour = np.append(met_one_hour,
                        met_filtered)
                pha_one_hour = np.append(pha_one_hour,
                        pha[mask])
                det_one_hour += [f"{det_headname}"]*met_filtered.size

            # empty
            if met_one_hour.size == 0:
                print("No photons...")
                continue

            # Barycentric Correction
            if args.barycor:
                print("...Barycentre correction...")
                orbit = fits.open(poshist)
                mjdreff = orbit[1].header['mjdreff']
                mjdrefi = orbit[1].header['mjdrefi']
                mjd = (met_one_hour/86400)+mjdreff+mjdrefi
                delta_t = barycor(mjd,
                        ra=args.ra,
                        dec=args.dec,
                        orbit=poshist,
                        jplephem=args.jplephem,
                        accelerate=args.accelerate)
                if delta_t is None:
                    continue
                TDB_one_hour = met_one_hour + delta_t
                OUT_DATA = [met_one_hour, pha_one_hour, TDB_one_hour, det_one_hour]
                OUT_COLN = ['TIME', "PHA", "TDB", "DET"]
            else:
                OUT_DATA = [met_one_hour, pha_one_hour]
                OUT_COLN = ['TIME', "PHA"]

            # Save data
            outname = os.path.join(args.output_dir, args.stem+"_{}_{}z.fits".format(next_day, str(hour).zfill(2)))
            print(f"...Saving to FITS {outname}...")
            #TODO store_pha, store_det, and energy selection not available
            save_to_fits(outname, OUT_DATA, OUT_COLN)
        next_day += one_day






if __name__ == "__main__":
    main()
