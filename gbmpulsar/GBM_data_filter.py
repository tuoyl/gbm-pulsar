#!/usr/bin/env python
# coding: utf-8

"""
The Tool is first development by Xiao Shuo, and then optimize by YL
GBM data selection tool for specific persistent X-ray source,
the selection cretira are the event list when the incident angle between
the source and the each GRD detetor are less than 70 degrees.

The input parameter should be:
    1. The coordinates of the source
    2. event file of specific GBM detector (one of 12 GRDs)
    3. the output file to save the filtered event file (if empty then no output)
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import argparse
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates.representation import CartesianRepresentation
from astropy.coordinates import SkyCoord
from astropy import units as u
from gbm.data.poshist import PosHist

def src_j2000_to_elv(ra, dec, x_sc, y_sc, z_sc):
    """
    calculate the ELV for given source with RA, DEC (J2000)

    Parameter
    ---------
    ra: float
        Right ascension
    dec: float
        Declination
    x_sc: array-like
        The x position in cartesian coordinates of spacecraft
    y_sc: array-like
        The y position in cartesian coordinates of spacecraft
    z_sc: array-like
        The z position in cartesian coordinates of spacecraft

    Return
    ------
    elv: array-like
        The elivation angle (<0 is blocked by Earth)
    """

    #x, y, z of source
    x_src, y_src, z_src = radec2xyz(ra, dec)

    #vector dot
    x_sc = -x_sc
    y_sc = -y_sc
    z_sc = -z_sc

    dot_value = x_src*x_sc + y_src*y_sc + z_src*z_sc

    #mol of vector_src and vector_sc
    mol_src = np.sqrt(x_src**2 + y_src**2 + z_src**2)
    mol_sc  = np.sqrt(x_sc**2  + y_sc**2  + z_sc**2)

    elv = np.rad2deg( np.arccos( dot_value/(mol_src * mol_sc) )) - 67.01571812044855

    return elv

def xyz2radec(xyz):
    """
    provided by Xiao Shuo
    convert from J2000 coordicate in Cartesian to RA, Dec
    #将J2000的xyz转换到radec

    Parameter
    ---------
    xyz: array
        the x,y,z value in Cartesian coordinates

    Return
    ------
    ra_dec: array
        the RA and DEC

    """
    radec=CartesianRepresentation(xyz[0],xyz[1],xyz[2]).represent_as(UnitSphericalRepresentation)
    ra=radec.lon.degree
    dec=radec.lat.degree
    ra_dec=np.array([ra,dec])
    return ra_dec

def radec2xyz(ra, dec):
    """
    convert RA, DEC to x, y z
    """
    c = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
    return (c.cartesian.x.value,
            c.cartesian.y.value,
            c.cartesian.z.value)

def UTC_TO_MET(utc,
        ref_utc='2001-01-01T00:00:00.000',
        leap_sec=5):
    """
    convert UTC to GBM MET

    Parameters
    ----------
    utc: string
        the format is '2021-08-12T16:00:00'

    ref_utc: string
        the reference time of mission (GBM is '2001-01-01T00:00:00.000')

    leap_sec: int
        the leap second, Only valid for data afeter 2017-01-01, and the leap second is 5 seconds for them.
    """
    ref_mjd = Time(ref_utc, scale='utc', format='isot').mjd
    met_gbm = (Time(utc).mjd-ref_mjd)*86400.+leap_sec   ###只对2017-1-1日后数据有效!!!+5是闰秒。统一换算到了GBM时间系统
    return round(met_gbm,6)                            #https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/xTime

def angle(x1,x2):
    """
    angle (degree) of two vector

    Parameters
    ----------
    x1, x2: array
        two vector. e.g., x1 = np.array([1,2,4])
    """
    return math.degrees(np.arccos(
        np.dot(x1, x2)/(np.linalg.norm(x1)*np.linalg.norm(x2))))

def filter(GBM_poshist, detector, radec, met, angle_incident=70, retrieve_mask=False):
    """
    select the events which has incident angle with the source less than
    creteria (default is 70 degrees)

    Parameters
    ----------
    GBM_poshist: string
        file of GBM poshist file (absolute path)

    detector: string
        one of the GBM detetor, optionals are
        ['n0','n1','n2','n3','n4','n5','n6','n7','n8','n9','na','nb']

    radec: array
        array for RA and Dec

    met: array
        The MET time of intersted events (GBM MET)

    angle_incident: float
        incident angle (degree) to select the data

    retrieve_mask: bool
        if True, then return the mask to filter the mask as well

    Returns
    -------
        met: array
            filtered MET array

        mask: array (bool)
            if retrieve_mask=True
    """

    poshist = PosHist.open(GBM_poshist)
    hdulist = fits.open(GBM_poshist)
    tstart = hdulist[1].header['TSTART']
    tstop  = hdulist[1].header['TSTOP']
    sclk_utc = hdulist[1].data['sclk_utc']
    pos_x  = hdulist[1].data['POS_X']
    pos_y  = hdulist[1].data['POS_Y']
    pos_z  = hdulist[1].data['POS_Z']


    mask = np.ones(met.size, dtype=bool)
    met_raw = met
    # make sure photons are with poshist range
    mask *= (met>tstart)&(met<tstop)
    met = met[mask]

    # NOTE: To reduce the consumption of computational expense,
    # the hourly GBM data is divided into 60 segments here,
    # and n each segment calculate one pointing angle in calculated,
    # if meets the requirements, whole segment is selected. It would be very slow to calculate each photon
    segments = 60
    met_edges = np.linspace(met.min(), met.max(), segments)
    radec_edges = poshist.detector_pointing(detector, met_edges)

    met_filtered = np.array([])
    pi_filtered  = np.array([])

    #mask for segments to merge
    mask_gtis = np.zeros(met_raw.size, dtype=bool)

    for i in range(segments-1):
        ## calculate the source pointing
        source_pointing = radec2xyz(radec_edges[0][i], radec_edges[1][i])

        ## calculate the incident angle between source and detector
        if angle(
                radec2xyz(radec[0],radec[1]),
                source_pointing) <= angle_incident:
            mask_gtis = mask_gtis | ((met_raw>=met_edges[i])&(met_raw<met_edges[i+1]))
    #        met_filtered = np.append(met_filtered,
    #                met_raw[mask])
    mask = (mask&mask_gtis)

    ## Filter ELV for mask
    mask_elv = np.zeros(met_raw.size, dtype=bool)
    for i in range(segments-1):
        idx_seg = np.searchsorted(sclk_utc, met_edges[i])
        elv = src_j2000_to_elv(radec[0], radec[1],
                               pos_x[idx_seg],
                               pos_y[idx_seg],
                               pos_z[idx_seg])
        if elv > 0:
            mask_elv = mask_elv | ((met_raw>=met_edges[i])&(met_raw<=met_edges[i+1]))
    mask = (mask&mask_elv)

    met_filtered = met_raw[mask]
    if retrieve_mask:
        return met_filtered, mask
    else:
        return met_filtered


if __name__ == "__main__":

    pos_file = "test/glg_poshist_all_220301_v00.fit"
    evt_file = "test/glg_tte_n0_220301_00z_v10.fit.gz"
    evt_file = "test/glg_tte_b0_220301_10z_v00.fit.gz"
    hdulist = fits.open(evt_file)
    met = hdulist[2].data.field("TIME")
    PI  = hdulist[2].data.field("PHA")
    bins = np.arange(met.min(), met.max(), 1)
    y, x = np.histogram(met, bins)
    plt.figure()
    plt.errorbar(x[:-1], y)
    clean_time, mask = filter(pos_file,
            'b0',
            np.array([83.22, 22.01]),
            met=met,
            angle_incident=70,
            retrieve_mask=True)
    print(clean_time, mask)
    y, x = np.histogram(clean_time, bins)
    plt.errorbar(x[:-1], y)
    plt.show()

