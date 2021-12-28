#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astroquery.gaia import Gaia

#%%
datadir = "data/"

#%%
d = pd.read_csv(datadir + "UltracoolSheet - Main.csv")

#%% non-M dwarfs with higher mass companions
idxcomp = (d.has_higher_mass_companion!="N") & (d.spt_opt.str.contains("M")!=True) & (d.spt_ir.str.contains("M")!=True) #& ((d.spt_opt==d.spt_opt)|(d.spt_ir==d.spt_ir))
d = d[idxcomp].reset_index(drop=True)

#%% <30'' separation
isep = d.sep_companion < 30
#list(d.keys())

#%%
d = d[isep].reset_index(drop=True)
d['num'] = d.index
len(d)

#%%
d.sort_values("sep_companion")[["has_higher_mass_companion", "sep_companion", "spt_opt", "spt_ir", "J_2MASS", "H_2MASS", "Ks_2MASS", "GaiaDetect_DR2", "sourceID_DR2", "G_DR2", 'z_P1']]

#%% Gaia EDR3 info
results = pd.DataFrame({})
for i in range(len(d)):
        print (i)
        s = d.iloc[i]
        coord = SkyCoord(ra=float(s.ra_j2000_formula)*u.degree, dec=float(s.dec_j2000_formula)*u.degree, pm_ra_cosdec=float(s.pmra_formula)*u.mas/u.yr, pm_dec=float(s.pmdec_formula)*u.mas/u.yr, frame="icrs", obstime=Time("J2000"), distance=s.dist_plx_formula*u.pc)
        new_coord = coord.apply_space_motion(new_obstime=Time("J2016.0"))
        radius = u.Quantity(35., u.arcsec)
        j = Gaia.cone_search_async(new_coord, radius, table_name='gaiaedr3.gaia_source')
        r = j.get_results()
        r = r.to_pandas()
        r['num'] = s.name
        results = results.append(r)
results = results.reset_index(drop=True)

#%%
d = pd.merge(d, results, on='num', how='left')

#%%
d.to_csv("targets_30arcsec_merged.csv", index=False)
