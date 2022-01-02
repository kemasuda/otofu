#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os, sys
#sys.path.append(os.path.join(os.path.dirname(os.path.abspath(os.getcwd())), 'checkfield'))
from otofu.obsplan import *
from texio import add_target

#%%
import seaborn as sns
sns.set(style='ticks', font_scale=1.6, font='Times New Roman')
from matplotlib import rc
rc('text', usetex=True)

#%%
d = pd.read_csv("targets_30arcsec_merged.csv")
d["Rmag"] = gaia_to_jcr(d.phot_g_mean_mag, d.bp_rp)

#%%
#times = ['2021-8-1 00:00:00', '2021-10-1 00:00:00', '2021-12-1 00:00:00', '2022-2-1 00:00:00', '2022-4-1 00:00:00', '2022-6-1 00:00:00', '2022-8-1 00:00:00'] # Subaru S21A
times = ['2021-3-27 00:00:00']
times = ['2021-5-27 22:00:00']
times = ['2021-7-16 02:30:00']

#%%
keys1 = ['name', 'spt_opt', 'spt_ir'] + ["sep_companion", "exoplanet", 'ra_j2000_formula', 'dec_j2000_formula', 'pmra_formula', 'pmdec_formula', 'J_2MASS', 'H_2MASS', 'Ks_2MASS', 'z_P1'] #['RAJ2000', 'DEJ2000']
keys2 = ['source_id', 'phot_g_mean_mag', 'Rmag',  'pmra', 'pmdec', 'parallax']

#%%
pdic = {'name': "name", 'spt_opt': "spectral type (optical)", 'spt_ir': "spectral type (IR)",
'sep_companion': "separation of companion (arcsec)", 'exoplanet': "planetary mass?",
'ra_j2000_formula': "RA (J2000)", 'dec_j2000_formula': "Dec (J2000)",
'pmra_formula': "RA proper motion (mas/yr)", 'pmdec_formula': "Dec proper motion (mas/yr)", 'J_2MASS': "2MASS $J$ (mag)", 'H_2MASS': "2MASS $H$ (mag)", 'Ks_2MASS': "2MASS $K$ (mag)", 'z_P1': "Pan-STARRS $z$ (mag)", 'source_id': "Gaia EDR3 source ID",  'phot_g_mean_mag': "Gaia $G$ (mag)", 'Rmag': "estimated $R$ (mag)", 'pmra': "Gaia RA proper motion (mas/yr)", 'pmdec': "Gaia Dec proper motion (mas/yr)", 'parallax': "Gaia parallax (mas)", 'dist': "J2016.0 distance from target (arcsec)"}#,  'dist_'+new_ref_time: "%s distance from target (arcsec)"%new_ref_time, 'ra_'+new_ref_time: "RA (%s)"%new_ref_time, 'dec_'+new_ref_time: "Dec (%s)"%new_ref_time}

#%%
dsort = d.sort_values("sep_companion").reset_index(drop=True)

#%%
names = np.array(dsort.name.drop_duplicates(keep='first'))
len(names)

#%%
jmag_cut = 14.+999
hmag_cut = 13.+999
rmag_cut = 16.
alt_cut = 40.-5.

#%%
alt_cut = 40
hmag_cut = 14.5

#%%
count = 0
tfile = open('targets_Jul2021.tex', 'w')
figdir = 'alts_Jul2021/'
if not os.path.exists(figdir):
    os.system("mkdir %s"%figdir)

#%%
tnames = []
for name in names:
        s = d[d.name==name].reset_index(drop=True)

        if s.sep_companion[0]==0.:
                print (s)
                continue

        #continue

        if s.H_2MASS[0] > hmag_cut:
            continue

        #if not (s.J_2MASS[0]<jmag_cut or s.H_2MASS[0]<hmag_cut):
        #        continue

        #if not (np.nanmin(s.Rmag)<rmag_cut):
        #        continue

        name = s.name[0].strip(" ")
        #ra, dec = s[["target_ra_"+new_ref_time, "target_dec_"+new_ref_time]].iloc[0]
        ra, dec = s[["ra_j2000_formula", "dec_j2000_formula"]].iloc[0]
        maxalts, fig = check_altitudes(name, ra, dec, times)
        if np.max(maxalts) < alt_cut:
                continue
        fig.savefig(figdir+"%s.png"%name, dpi=300, bbox_inches="tight")

        if len(s)>1:
                s = s.sort_values("phot_g_mean_mag").reset_index(drop=True)
                targ = np.argmin(s.dist)
                idx = [targ] + [i for i in range(len(s)) if i != targ]
                s = s.iloc[idx].reset_index(drop=True)
                nsources = len(s)
                if nsources > 4:
                        nsources = 4
        else:
                nsources = 1
        if tfile is not None:
                add_target(s[:nsources], name, tfile, keys1, keys2, pdic)
        count += 1
        print (s.name[0])
        tnames.append(name)

if tfile is not None:
        tfile.close()
