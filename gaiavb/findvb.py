#%%
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
from otofu.obsplan import check_altitude

#%%
def add_vbinfo(data):
    data['mg1'] = data['phot_g_mean_mag1'] + 5*np.log10(data['parallax1']) - 10
    data['mg2'] = data['phot_g_mean_mag2'] + 5*np.log10(data['parallax2']) - 10
    data['delG'] = data.phot_g_mean_mag2 - data.phot_g_mean_mag1
    data['delC'] = np.abs(data.bp_rp2 - data.bp_rp1)
    data['darcsec'] = data.pairdistance*3600.
    data['draproj'] = (data.ra2 - data.ra1)*np.cos(data.dec1*np.pi/180.)
    data['ddecproj'] = (data.dec2 - data.dec1)
    data['sep'] = np.sqrt(data.draproj**2 + data.ddecproj**2)*3600.
    data['PA'] = np.arctan2(data.draproj, data.ddecproj)*180./np.pi
    data['poserror'] = np.sqrt(data.ra_error1**2 + data.dec_error1**2 + data.ra_error2**2 + data.dec_error2**2 \
            + (10*data.pmra_error1)**2 + (10*data.pmdec_error1)**2 + (10*data.pmra_error2)**2 + (10*data.pmdec_error2)**2)
    return data

def search_twomass(ra, dec, source_id1, radius=5.):
    d2m = pd.DataFrame({})
    for i in range(len(ra)):
        result = Vizier.query_region(SkyCoord(ra=ra[i], dec=dec[i], unit=(u.deg, u.deg), frame='icrs'),
                radius=radius*u.arcsec, catalog=["II/246/out"])
        try:
            _d = result["II/246/out"].to_pandas()
        except:
            continue
        _d['source_id1'] = source_id1[i]
        d2m = d2m.append(_d)
    d2m = d2m.rename({"_2MASS": "2MASS"}, axis='columns')
    return d2m

def select_targets(time, alt_min=45, gmag_max=8, deltagmag_min=1, deltagmag_max=4, sep_max=30., twomass_rad=5., R_max=0.1):
    data = pd.read_pickle("eb21binary_g11.pkl")
    data = add_vbinfo(data)

    idx = (data.R_chance_align<R_max) & (deltagmag_min<data.delG) & (data.delG<deltagmag_max) & (data.darcsec<sep_max) & (data.phot_g_mean_mag1<gmag_max)
    dtarg = data[idx].reset_index(drop=True)
    print ("# %d potential targets found."%len(dtarg))

    print ("checking altitudes at %s..."%time)
    alts = []
    for ra, dec in zip(dtarg.ra1, dtarg.dec1):
        _alt, _, _  = check_altitude(ra, dec, time)
        alts.append(_alt.to_value(u.deg))
    dtarg['alt'] = alts

    idxtarg = (dtarg.alt > alt_min)
    dvis = dtarg[idxtarg].reset_index(drop=True)
    print ("# %d visible targets found."%len(dvis))

    print ("retrieving 2MASS info...")
    d2m = search_twomass(dvis.ra1, dvis.dec1, dvis.source_id1, radius=twomass_rad)
    d2m_nodup = d2m.drop_duplicates("source_id1", keep='first')
    print ("# %d 2MASS sources found within %d arcsec (%d without duplication)."%(len(d2m), twomass_rad, len(d2m_nodup)))

    d = pd.merge(dvis, d2m.drop_duplicates("source_id1", keep='first')[['source_id1', '2MASS', "Jmag", "Hmag", "Kmag"]], on='source_id1', how='left')

    return d

#%%
time = '2021-8-20 02:30:00'

#%%
d = select_targets(time)

#%%
keys = ['2MASS', 'delG', 'sep', 'PA', 'phot_g_mean_mag1', 'Jmag', 'Hmag', 'Kmag', 'source_id1', 'ra1', 'dec1', 'alt', 'poserror']
dict = {"delG": "delta G mag", 'sep': 'separation (arcsec)', 'PA': 'PA (deg)', 'poserror': 'poserror (mas)', 'alt': 'altitude at %s'%time}

#%%
dout = d.sort_values("sep")[keys].rename(dict, axis='columns')
dout.to_csv("vbs.csv", index=False)
