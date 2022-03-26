#!/usr/bin/env python
# coding: utf-8

import numpy as np
from GBM_data_filter import UTC_TO_MET, filter
import argparse
import datetime
import glob
import os
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
"""


NOTICE = """

/---------------------------------------------\\
|   GBM data analysis tool for X-ray pulsar   |
\---------------------------------------------/

Usage:

        python GBM_pusar_pipeline.py --gbm_dir="/path/to/GRM/data" \\
                -tstart='2022-07-27T00:00:00' --tstop=''2022-07-28T00:00:00' \\
                --output_dir="/path/to/save/data" \\
                --stem="gbmCrab" --ra=83.883225 --dec=22.014458333333334

    some optional:

        python GBM_pusar_pipeline.py --gbm_dir="/path/to/GRM/data" \\
                -tstart='2022-07-27T00:00:00' --tstop=''2022-07-28T00:00:00' \\
                --output_dir="/path/to/save/data" \\
                --stem="gbmCrab" --ra=83.883225 --dec=22.014458333333334 \\
                --energylow=8 --energyhigh=25 \\
                --store_pha --store_det
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

    wildcast = os.path.join(gbm_data_dir, year, month, day,
                "glg_tte_n*_{}{}{}_{}z_v00.fit.gz".format(year[2:], month, day, hour))
    event_list = glob.glob(
            os.path.join(gbm_data_dir, year, month, day,
                "glg_tte_n*_{}{}{}_{}z_v00.fit.gz".format(year[2:], month, day, hour)))
    poshist = glob.glob(
            os.path.join(gbm_data_dir, year, month, day,
                "glg_poshist_all_{}{}{}_v00.fit".format(year[2:], month, day)))
    if len(poshist) is not 0:
        poshist = poshist[-1]

    return event_list, poshist

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
    det_headname = hdulist[2].header['DETNAM']
    met = hdulist[2].data.field("TIME")
    pha = hdulist[2].data.field("PHA")
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
            det_one_hour = np.array([])

            for evtfile in tqdm(evtfiles):
                if os.path.basename(evtfile).split('_')[4] != str(hour).zfill(2)+'z':
                    continue
                det_shortname, det_headname, met, pha = _resolve_evtname(evtfile)
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
                det_one_hour = np.append(det_one_hour,
                        np.ones(mask.size)*int(det_headname.split('_')[1]))

            # Save data
            outname = os.path.join(args.output_dir, args.stem+"_{}_{}z.fits".format(next_day, str(hour).zfill(2)))
            #TODO store_pha, store_det, and energy selection not available
            save_to_fits(outname, [met_one_hour, pha_one_hour], ['TIME', "PHA"])

        next_day += one_day






if __name__ == "__main__":
    main()
