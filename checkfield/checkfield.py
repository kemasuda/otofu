#%%
import numpy as np
import pandas as pd
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.mast import Catalogs
from astroquery.gaia import Gaia
from astropy.time import Time
from astroplan.plots import plot_finder_image

import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style='ticks', font_scale=1.6*0.5, font='sans-serif')
plt.rcParams["figure.figsize"] = (14*0.5,7*0.5)
from matplotlib import rc
rc('text', usetex=False)


#%%
def query_simbad(name):
    print ("# querying %s in simbad..."%name)
    result_table = Simbad.query_object(name)
    print (result_table)
    c = SkyCoord(ra=result_table["RA"], dec=result_table["DEC"], frame='icrs', unit=(u.hourangle, u.deg), obstime=Time("J2000.0"))
    return c

def search_gaia_cone(ra, dec, radius):
    coord = SkyCoord(ra=float(ra)*u.degree, dec=float(dec)*u.degree, frame="icrs", obstime=Time("J2016.0"))
    j = Gaia.cone_search_async(coord, u.Quantity(radius, u.arcsec), table_name='gaiaedr3.gaia_source')
    r = j.get_results().to_pandas()
    return r

def search_gaia_square(coord, fov, table="gaiadr2.gaia_source", nmax=100000):
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
