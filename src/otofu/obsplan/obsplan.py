
__all__ = ["query_simbad", "plot_FC", "plot_radius_hw", "check_field", "check_altitudes", "gaia_to_jcr"]

#%%
import numpy as np
import pandas as pd
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.mast import Catalogs
from astropy.time import Time
import matplotlib.pyplot as plt

#%%
def query_simbad(name):
    print ("# querying %s in simbad..."%name)
    result_table = Simbad.query_object(name)
    print (result_table)
    c = SkyCoord(ra=result_table["RA"], dec=result_table["DEC"], frame='icrs', unit=(u.hourangle, u.deg), obstime=Time("J2000.0"))
    return c

def search_gaia_cone(ra, dec, radius):
    from astroquery.gaia import Gaia
    coord = SkyCoord(ra=float(ra)*u.degree, dec=float(dec)*u.degree, frame="icrs", obstime=Time("J2016.0"))
    j = Gaia.cone_search_async(coord, u.Quantity(radius, u.arcsec), table_name='gaiaedr3.gaia_source')
    r = j.get_results().to_pandas()
    return r

def search_gaia_square(coord, fov, table="gaiadr2.gaia_source", nmax=100000):
    from astroquery.gaia import Gaia
    Gaia.MAIN_GAIA_TABLE = table
    Gaia.ROW_LIMIT = nmax
    gaia_sources = Gaia.query_object_async(coordinate=coord, width=fov, height=fov).to_pandas()
    if len(gaia_sources)==nmax:
        print ("# more sources found than nmax=%d."%nmax)
    print ("RA range: %.2fdeg" % (np.max(gaia_sources.ra)-np.min(gaia_sources.ra)))
    print ("Dec range: %.2fdeg" % (np.max(gaia_sources.dec)-np.min(gaia_sources.dec)))
    return gaia_sources

def serach_tic_using_gaiadr2(source_ids):
    tic_sources = Catalogs.query_criteria(catalog="TIC", GAIA=source_ids)
    return tic_sources.to_pandas()

def plot_FC(coord, fov, name, survey="DSS"):
    from astroplan.plots import plot_finder_image
    plt.figure(figsize=(8,6))
    ax, hdu=plot_finder_image(coord, survey=survey, reticle=True, fov_radius=fov)
    plt.title(r'%s %s'%(name, survey))
    plt.savefig("%s_fov.png"%name, dpi=200, bbox_inches="tight")

def plot_radius_hw(d, name):
    plt.figure()
    plt.xlabel("$H_w$ mag")
    plt.ylabel("stellar radius ($R_\odot$)")
    plt.yscale("log")
    plt.plot(d.Hwmag, d.rad, '.')
    plt.title("stars around %s"%name)
    plt.savefig("%s_rad_hw.png"%name, dpi=200, bbox_inches="tight");

def check_field(name, fov, parallax_cut=2., plot=False):
    """
    if type(target)=='str':
        coord = query_simbad(name)[0]
    elif isinstance(target, SkyCoord):
        coord = target
    else:
        print ("# input should be a simbad name or SkyCoord object.")
        return None
    """
    coord = query_simbad(name)[0]

    gaia_sources = search_gaia_square(coord, fov)
    gaia_sources = gaia_sources[gaia_sources.parallax > parallax_cut].reset_index(drop=True)
    tic_sources = serach_tic_using_gaiadr2(gaia_sources.source_id)
    tic_sources["source_id"] = np.array(tic_sources.GAIA).astype(np.int)
    tic_sources = tic_sources.rename(columns={'ra': 'ra_tic', 'dec': 'dec_tic'})

    d = pd.merge(gaia_sources, tic_sources, on='source_id')
    d["mg"] = d.phot_g_mean_mag + 5 * np.log10(d.parallax) - 10
    d["Hwmag"] = 0.7 * d.Jmag + 0.3 * d.Hmag

    if plot:
        plot_FC(coord, fov, name)
        plot_radius_hw(d, name)

    return coord, d

#%%
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.coordinates import get_sun, get_moon
from matplotlib.dates import DateFormatter
# altitude from Keck
def check_altitudes(name, ra, dec, times, maxalt=False):
    maxalts = []

    keck = EarthLocation.of_site('Keck Observatory')
    utcoffset = -10*u.hour
    c = SkyCoord(ra=ra*u.degree, dec=dec*u.degree)

    fig = plt.figure(figsize=(28,12))

    for i,time in enumerate(times):
            dt = np.linspace(-8, 8, 50)
            obstimes_local = Time(time) + dt*u.hour
            obstimes = obstimes_local - utcoffset
            obsframes = AltAz(obstime=obstimes, location=keck)
            alts = c.transform_to(obsframes)
            moonalts = get_moon(obstimes).transform_to(obsframes)
            sunalts = get_sun(obstimes).transform_to(obsframes)

            tarr = np.array([np.datetime64(t) for t in obstimes_local.value])
            #fig = plt.figure(figsize=(12,8)) # single
            fig.add_subplot(2, 4, i+1)
            plt.ylim(0, 90)
            plt.ylabel('altitude (deg)')
            plt.xlabel('local time (HST: UTC$%d$)'%utcoffset.value)
            plt.xlim(tarr[0], tarr[-1])
            plt.plot(tarr, alts.alt, label='%s'%name, lw=3, alpha=0.8)
            plt.plot(tarr, sunalts.alt, color='k', lw=1, label='sun', zorder=-100)
            plt.plot(tarr, moonalts.alt, color='gray', ls='dashed', lw=1, label='moon', zorder=-100)
            plt.fill_between(tarr, 0, 90, sunalts.alt>-18*u.deg, color='gray', zorder=-1000, alpha=0.4, lw=0.5)
            plt.fill_between(tarr, 0, 90, sunalts.alt>-0*u.deg, color='k', zorder=-1000, alpha=0.2, lw=0.5)
            plt.gca().xaxis.set_major_formatter(DateFormatter('%m/%d %H:%M'))
            for label in plt.gca().get_xticklabels():
                    label.set_rotation(30)
                    label.set_horizontalalignment('right')
            #if i==6:
            #        plt.legend(loc='upper left', bbox_to_anchor=(1,1))

            if maxalt:
                maxalts.append(np.max(alts.alt[sunalts.alt<-18*u.deg].value))
            else:
                alts_specified = c.transform_to(AltAz(obstime= Time(time) - utcoffset, location=keck))
                maxalts.append(np.max(alts_specified.alt.value))

    handler, label = plt.gca().get_legend_handles_labels()
    fig.legend(handler, label, loc='upper left', bbox_to_anchor=(0.76,0.47), prop={'size':20})
    fig.tight_layout(h_pad=2., w_pad=0.05)

    return np.array(maxalts), fig

#%%
def check_altitude(ra, dec, time):
    keck = EarthLocation.of_site('Keck Observatory')
    utcoffset = -10*u.hour
    c = SkyCoord(ra=ra*u.degree, dec=dec*u.degree)

    obstimes_local = Time(time)
    obstimes = obstimes_local - utcoffset
    obsframes = AltAz(obstime=obstimes, location=keck)
    alts = c.transform_to(obsframes)
    moonalts = get_moon(obstimes).transform_to(obsframes)
    sunalts = get_sun(obstimes).transform_to(obsframes)

    return alts.alt, moonalts.alt, sunalts.alt

#%% conversion from Gaia mag to JC R mag
def gaia_to_jcr(g, bprp):
    gminr = -0.003226 + 0.3833*bprp - 0.1345*bprp**2
    return g - gminr


#%%
"""
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style='ticks', font_scale=1.6, font='sans-serif')
plt.rcParams["figure.figsize"] = (14,7)
from matplotlib import rc
rc('text', usetex=False)

#%%
fov = 0.45*u.deg

#%%
name = "TOI-700"
name = "LTT 1445A"
name = "TRAPPIST-1"

#%%
coord, d = check_field(name, fov, plot=True)

#%%
#from astroquery import skyview
#sk = skyview.SkyViewClass()
#sk.list_surveys()
"""
