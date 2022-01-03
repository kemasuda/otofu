#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from astropy.constants import R_sun, R_earth
from astroquery.mast import Catalogs

#%% confirmed planets
dconf = pd.read_csv('data/confirmed_20210427.csv', comment='#')
dconf['disc_telescope'].iloc[np.array(dconf.disc_telescope!=dconf.disc_telescope)] = "N/A"
dconf_tess = dconf[dconf.disc_telescope.str.contains("TESS")].reset_index(drop=True) # TESS confirmed
dconf_tra = dconf[(~dconf.disc_telescope.str.contains("TESS"))&(dconf.tran_flag==1.)].reset_index(drop=True) # other transits

#%% TID missing for some of the confirmed-tranet hosts
func = lambda x: int(x[4:]) if x==x else -1
dconf_tra['tid'] = list(map(func, dconf_tra.tic_id))
dconf_tess['tid'] = list(map(func, dconf_tess.tic_id))
print ('# %d stars with TIC IDs.'%len(dconf_tra[dconf_tra.tid>0].drop_duplicates('tid', keep='first')))

#%% TIC information for the confirmed-tranet hosts
catalog_data = Catalogs.query_criteria(catalog="TIC", ID=dconf_tra.tid).to_pandas()
catalog_data['tid'] = np.array(catalog_data.ID).astype(int)
print ("# %d TIC stars found."%len(catalog_data))
dconf_tra = pd.merge(dconf_tra, catalog_data[["tid", "mass", "rad", "Teff", "Jmag", "Hmag"]], on='tid').reset_index(drop=True)

#%% known TOIs
dtoi = pd.read_csv('data/TOI_20210427.csv', comment='#')
print ('# %d stars with TOIs.'%len(dtoi.drop_duplicates("tid", keep='first')))

#%% TIC information for the TOIs
catalog_data = Catalogs.query_criteria(catalog="TIC", ID=dtoi.tid).to_pandas()
#catalog_data2 = Catalogs.query_criteria(catalog="ctl", ID=dtoi.tid).to_pandas()
catalog_data['tid'] = np.array(catalog_data.ID).astype(int)
#catalog_data2['tid'] = np.array(catalog_data2.ID).astype(int)
print ("# %d TIC stars found."%len(catalog_data))
dtoi = pd.merge(dtoi, catalog_data[["tid", "mass", "rad", "Teff", "Jmag", "Hmag"]], on='tid').reset_index(drop=True)

#%% fill nan in stellar/planetary radii
dconf_tra['rad'] = dconf_tra['rad'].fillna(dconf_tra['st_rad'])
dtoi['rad'] = dtoi['rad'].fillna(dtoi['st_rad'])
dconf_tra['pl_rade'] = dconf_tra['pl_rade'].fillna(dconf_tra.rad*R_sun*np.sqrt(dconf_tra.pl_trandep*1e-2)/R_earth)
dtoi['pl_rade'] = dtoi['pl_rade'].fillna(dtoi.rad*R_sun*np.sqrt(dtoi.pl_trandep*1e-6)/R_earth)

#%% missing Teff for TRAPPIST-1
print (dconf_tra.Teff[dconf_tra.hostname=="TRAPPIST-1"])
dconf_tra.Teff[dconf_tra.hostname=="TRAPPIST-1"] = 2566
print (dconf_tra.Teff[dconf_tra.hostname=="TRAPPIST-1"])

#%%
dconf_tess['mass'] = dconf_tess.st_mass
dconf_tess['rad'] = dconf_tess.st_rad
dconf_tess['Hmag'] = dconf_tess.sy_hmag
dconf_tess.mass[dconf_tess.tid==259377017] = 0.362274

#%% remove TOIs duplicated in confirmed TESS/tranet catalogs
dtoi['pl_name'] = dtoi.toi
dtoi['hostname'] = ["TOI-%d"%dtoi.toipfx[i] for i in range(len(dtoi))]
dtoinodup = dtoi[(~dtoi.tid.isin(dconf_tess.tid))&(~dtoi.tid.isin(dconf_tra.tid))].reset_index(drop=True)

#%% confirmed tranets from Kepler/K2 and others
kepidx = dconf_tra.disc_facility.str.contains("Kepler")+dconf_tra.disc_facility.str.contains("K2")
dconf_kepler = dconf_tra[kepidx].reset_index(drop=True)
dconf_other = dconf_tra[~kepidx].reset_index(drop=True)

#%%
# dtoinodup: TESS candidates (TOIs) - TESS confimed planets
# dconf_tess: TESS confirmed planets
# dconf_kepler: confirmed planets from Kepler/K2
# dconf_other: other confirmed planets
datasets = dict(zip(['TESS candidate', 'TESS confirmed', "Kepler/K2", "ground"], [dtoinodup, dconf_tess, dconf_kepler, dconf_other]))
datasets.keys()

#%%
if os.path.exists("compiled_catalog.h5"):
    os.system("rm compiled_catalog.h5")
store = pd.HDFStore('compiled_catalog.h5')
for k in datasets.keys():
    store[k] = datasets[k]
store.close()

#%% read
#d = pd.HDFStore('compiled_catalog.h5')
#d.keys()
