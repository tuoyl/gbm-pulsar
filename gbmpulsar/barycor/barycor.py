from astropy.io import fits
import numpy as np
from jplephem.spk import SPK
from gbmpulsar.barycor.tdb2tdt import tdb2tdt
from scipy.interpolate import interp1d
import argparse
import os

# To run the code one needs to have access to one of the ephemeride files provided by jpl:
# https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/
# here I use the DE421 file


c = 299792.458 # km/s  (kilometers per second)

def barycor(date,ra,dec, orbit=False, return_correction=True, approx_einstein=10, jplephem=None, accelerate=False):
    """Apply barycentric correction to vector
       date (in MJD) assuming coordinates ra,dec given in degrees
       optionally use orbit of a satellinte in heasarc format
       This includes geometric and shapiro corrections, but einstein
       correction is not implemented for now. Accurate to about 3e-4s.
       If return_correction=True (default), correction in seconds is returned.
       approx einstein is used to define number of points to calculate interpolated einstein correction (0 means do for each data point).
    see also https://asd.gsfc.nasa.gov/Craig.Markwardt/bary/ for detail and https://pages.physics.wisc.edu/~craigm/idl/ephem.html for the codebase
    which at least inspired the code below.
    v0.2: 21.03.2022 by V. Doroshenko
       """
    ra = np.radians(np.double(ra))
    dec = np.radians(np.double(dec))


    if accelerate:
        date_raw = date
        print("Accelerating barycor")
        ## Draw photons for accelerate
        N_segment = int(date_raw.size/60) # divide time into 60 segnents
        date = np.linspace(date_raw.min(), date_raw.max(), N_segment)

    jd =  np.array(date,dtype=np.float64) + 2400000.5
    if not os.path.exists(jplephem):
        raise FileNotFoundError("File {} not found".format(jplephem))
    kernel = SPK.open(jplephem)


    msol = 0.0002959122082855911

    # get positions and velocities of the sun, earth-moon barycenter and earth
    (x_sun,y_sun,z_sun),(vx_sun,vy_sun,vz_sun) = kernel[0,10].compute_and_differentiate(jd)
    (x_em,y_em,z_em),(vx_em,vy_em,vz_em) = kernel[0,3].compute_and_differentiate(jd)
    (x_e,y_e,z_e),(vx_e,vy_e,vz_e) = kernel[3,399].compute_and_differentiate(jd)


    x_earth = x_em + x_e
    y_earth = y_em + y_e
    z_earth = z_em + z_e

    vx_earth = (vx_em + vx_e)/86400 # velocities are always distance units per day as per jplephem documentation, so need to divide to get km/s
    vy_earth = (vy_em + vy_e)/86400
    vz_earth = (vz_em + vz_e)/86400


    if orbit:
        orbit = fits.open(orbit)
        mjdref = orbit[1].header['mjdreff']+orbit[1].header['mjdrefi']
        #try:
        #    minmet = (np.min(date)-1.-mjdref)*86400
        #    maxmet =(np.max(date)+1.-mjdref)*86400
        #except:
        #    minmet = (date-1-mjdref)*86400
        #    maxmet = (date+1-mjdref)*86400
        minmet = (np.min(date)-mjdref)*86400
        maxmet =(np.max(date)-mjdref)*86400

        #try:
        #    t = orbit[1].data.field('time')
        #except:
        #    t = orbit[1].data.field('sclk_utc')
        try:
            t = orbit[1].data.field('sclk_utc')
        except:
            t = orbit[1].data.field("TIME")

        #mi, ma = t.searchsorted([minmet,maxmet])
        #t = t[mi:ma]/86400.+mjdref
        mask = (t>minmet-1)&(t<maxmet+1)
        t = t[mask]/86400 + mjdref
        if t.size == 0:
            print(f">>> Orbit Time not cover data")
            return None

        # interpolate orbit to observed time and convert to km and km/s
        if 'pos_x' in [x.lower() for x in orbit[1].data.names]:
            x_s = np.interp(date, t, orbit[1].data.field('pos_x')[mask]/1000.)
            y_s = np.interp(date, t, orbit[1].data.field('pos_y')[mask]/1000.)
            z_s = np.interp(date, t, orbit[1].data.field('pos_z')[mask]/1000.)

            vx_s = np.interp(date, t, orbit[1].data.field('vel_x')[mask]/1000.)
            vy_s = np.interp(date, t, orbit[1].data.field('vel_y')[mask]/1000.)
            vz_s = np.interp(date, t, orbit[1].data.field('vel_z')[mask]/1000.)
        elif 'x' in [x.lower() for x in orbit[1].data.names]:
            x_s = np.interp(date, t, orbit[1].data.field('x')[mask]/1000.)
            y_s = np.interp(date, t, orbit[1].data.field('y')[mask]/1000.)
            z_s = np.interp(date, t, orbit[1].data.field('z')[mask]/1000.)

            vx_s = np.interp(date, t, orbit[1].data.field('vx')[mask]/1000.)
            vy_s = np.interp(date, t, orbit[1].data.field('vy')[mask]/1000.)
            vz_s = np.interp(date, t, orbit[1].data.field('vz')[mask]/1000.)



        x_obs, y_obs, z_obs = x_earth + x_s, y_earth + y_s, z_earth + z_s
        vx_obs, vy_obs, vz_obs = vx_earth + vx_s, vy_earth + vy_s, vz_earth + vz_s
        # orbital correction
        ocor = (vx_earth*x_s+vy_earth*y_s+vz_earth*z_s)/c**2
    else:
        x_obs, y_obs, z_obs = x_earth, y_earth , z_earth
        vx_obs, vy_obs, vz_obs = vx_earth , vy_earth , vz_earth
        ocor = 0.

    # #components of the object unit vector:
    x_obj = np.cos(dec)*np.cos(ra)
    y_obj = np.cos(dec)*np.sin(ra)
    z_obj = np.sin(dec)

    #geometric correction
    geo_corr = (x_obs*x_obj + y_obs*y_obj + z_obs*z_obj)/c

    #einstein correction
    if approx_einstein == 0:
        einstein_corr = tdb2tdt(jd)
    else:
        xx = np.linspace(jd.min(),jd.max(),approx_einstein)
        einstein_corr = tdb2tdt(xx)
        einstein_corr = np.interp(jd,xx,einstein_corr)
    #shapiro correction ("Shapiro") = - (2 G Msun/c^3) log(1 + cos th)
    sun_dist = np.sqrt((x_sun-x_obs)**2+(y_sun-y_obs)**2+(z_sun-z_obs)**2)
    costh = ((x_obs-x_sun)*x_obj+(y_obs-y_sun)*y_obj + (z_obs-z_sun)*z_obj)/sun_dist
    shapiro_corr = - 9.8509819e-06*np.log(1.+costh)
    corr = geo_corr + ocor + einstein_corr - shapiro_corr

    ## Interpolate for accelerated photons
    if accelerate:
        ## interpolate
        corr_fun = interp1d(date, corr, kind='quadratic')
        corr = corr_fun(date_raw)

    if return_correction:
        return corr
    else:
        return date + corr/86400.

def cp_table(old_table):
    col_names = old_table.names
    col_type = old_table.formats
    cp_col = []
    for i in range(len(col_names)):
        cp_col.append( fits.Column(name=col_names[i], array=old_table.field(col_names[i]), format=col_type[i]) )
    new_table = fits.BinTableHDU.from_columns(cp_col)
    return new_table
def write_file(datafile, tdb):
    #read old
    hdulist_old = fits.open(datafile)
    prim_hdr_old = hdulist_old[0].header
    prim_hdr_new = fits.PrimaryHDU(header=prim_hdr_old)
    hdr_all = []
    for i in range(len(hdulist_old)):
        hdr_all.append(hdulist_old[i].header)
    #modify main table
    table1 = hdulist_old[1].data
    col_names = table1.names
    col_type = table1.formats
    cp_col = []
    if 'TDB' not in col_names:
        print("...adding a column to event file...")
        new_col = fits.Column(name='TDB', array=tdb, format='1D')
        for i in range(len(col_names)):
            cp_col.append( fits.Column(name=col_names[i], array=table1.field(col_names[i]), format=col_type[i]) )
        cp_col.append(new_col)
        new_table1 = fits.BinTableHDU.from_columns(cp_col)
        if len(hdulist_old) <=2:
            hdulist_new = fits.HDUList([prim_hdr_new, new_table1])
        else:
            table_rest = []
            for i in np.arange(2, len(hdulist_old), 1):
                table_rest.append(cp_table(hdulist_old[i].data))
            hdulist_new = fits.HDUList([prim_hdr_new, new_table1]+table_rest)
        #update header and info
        hdr_all[1].append('TTYPE' + str(hdr_all[1]['TFIELDS']), 'TDB', 'label for field')
        hdr_all[1].append('TFORM' + str(hdr_all[1]['TFIELDS']), '1D', 'format of field')
        for i in range(len(hdulist_new)):
            hdulist_new[i].header = hdr_all[i]
        hdulist_new.writeto(datafile, overwrite=True)
        hdulist_new.close()
        hdulist_old.close()
    else:
        print("...column TDB exists, overwrite the column...")
        table1['TDB'][:] = tdb
        hdulist_old.writeto(datafile, overwrite=True)
        hdulist_old.close()
    return


def parse_args():

    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description='barycor --evtfile="/path/to/eventfile" --extnum=1 --colnum=1')

    parser.add_argument(
            "--evtfile",
            type=str,
            required=True,
            help="The filename of L1 Event file")

    parser.add_argument(
            "--orbitfile",
            type=str,
            required=True,
            help="The filename of orbit file")

    parser.add_argument(
            "--extnum",
            type=int,
            required=True,
            help="The extension number to specify the Time (Start from 0)")

    parser.add_argument(
            "--colnum",
            type=int,
            required=True,
            help="The column number to specify the Time (Start from 0)")

    parser.add_argument(
            "--ra",
            type=float,
            required=True,
            help="RA of source")

    parser.add_argument(
            "--dec",
            type=float,
            required=True,
            help="Declination of the source")

    parser.add_argument(
            "--jplephem",
            type=str,
            required=True,
            help="the JPL Ephemeris file, files could be downloaded from \
                    https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/")

    parser.add_argument(
            "--accelerate",
            action="store_true",
            required=False,
            default=False,
            help="to accelerate the barycentric correction processes by barycor sampling photons, and interpolate\
                    the barycored arrival times for whole time series")

    return parser.parse_args()

def main():
    args = parse_args()

    hdulist = fits.open(args.evtfile)
    met = hdulist[args.extnum].data.field(args.colnum)
    orbit = fits.open(args.orbitfile)
    mjdref = orbit[1].header['mjdrefi'] + orbit[1].header['mjdreff']
    mjd = met/86400. + mjdref

    dtime = barycor(mjd,
            ra=args.ra,
            dec=args.dec,
            orbit=args.orbitfile,
            jplephem=args.jplephem,
            accelerate=args.accelerate)

    print("Finish bary correction...")
    print("Writing file...")
    write_file(args.evtfile, met+dtime)

if __name__ == "__main__":
    main()
