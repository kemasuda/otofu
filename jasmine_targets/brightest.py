#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#%%
import seaborn as sns
sns.set(style='ticks', font_scale=1.6, font='times')
#sns.set(style='whitegrid', font_scale=1.6, font='sans-serif')
from matplotlib import rc
rc('text', usetex=True)

#%%
from astropy.constants import R_sun, R_earth
from otofu.jasmine import noise_model, HZperiod, m2t
def assign_snHZ(d, mass_key, rad_key, teff_key, jmag_key, hmag_key):
    m, r, t = d[mass_key], d[rad_key], d[teff_key]
    d['Hwmag'] = d[jmag_key]*0.7 + d[hmag_key]*0.3
    d['pHZ'] = HZperiod(m, radius=r, teff=t, s=0.85) # inner edge
    d['durHZ'] = (13./24.) * (d.pHZ/365.25)**(1./3.) * r/m**(1./3.) * np.pi/4.
    d['ptraHZ'] = 1 / (3.7528*d.pHZ**(2./3.)*m**(1./3.)/r)
    d['depth_earth'] = r**(-2) * (R_earth/R_sun).value**2
    d['snHZ'] = d.depth_earth*1e6 / noise_model(d.Hwmag) * np.sqrt(d.durHZ*1440/5./2.) # last factor 2 from Earth occultation
    #d['snHZ_tess'] = d.depth_earth*1e6/tessnoise_5min(d.Tmag)*np.sqrt(d.durHZ*1440/5.)

#%%
# dtoinodup: TESS candidates (TOIs) - TESS confimed planets
# dconf_tess: TESS confirmed planets
# dconf_kepler: confirmed planets from Kepler/K2
# dconf_other: other confirmed planets
d = pd.HDFStore('compiled_catalog.h5')
datakeys = [k[1:] for k in d.keys()]
datakeys = list(np.array(datakeys)[[0,1,3,2]]) # change order
datadict = dict(zip(datakeys, [d[key] for key in datakeys]))
datalist = [d[key] for key in datakeys]
print (datakeys)

#%% CTL dwarfs
d = pd.read_csv('data/CTLv8_t4000_h11.csv', comment='#')
d = d[d.rad<0.6].reset_index(drop=True)
_ = assign_snHZ(d, 'mass', 'rad', 'Teff', 'Jmag', 'Hmag')

#%%
def simulate(ds):
    rands = np.random.rand(len(ds))
    idx = rands < ds.ptraHZ
    return ds[idx].reset_index(drop=True)

#%%
Nsim = 10000
sn_threshold = 7
res = pd.DataFrame({})
for i in range(Nsim):
    dsim = simulate(d)
    res = res.append(dsim[dsim.snHZ>sn_threshold][["Hmag", "mass"]])
    #res = res.append(dsim[["Hmag", "mass", "snHZ"]])

#%%
#x, y = np.array(res.Hmag[res.snHZ>10]), np.array(res.mass[res.snHZ>10])
x, y = np.array(res.Hmag), np.array(res.mass)
medges = np.arange(0.05, 0.55, 0.07)
medges = np.arange(0.08, 0.60, 0.08)
mvals = 0.5 * (medges[1:] + medges[:-1])

#%%
hmag1, hmag3, hmag10 = [], [], []
for mu, ml in zip(medges[1:], medges[:-1]):
    _idx = (ml<=y) & (y<mu)
    plt.figure()
    plt.title("$%.2f$-$%.2f\,M_\odot$"%(ml, mu))
    plt.xlabel("$H_w$ mag")
    plt.ylabel("$N_\mathrm{exp}$ of transiting HZ Earths")
    #plt.hist(x[_idx], cumulative=True, weights=w[_idx])
    sortx = np.sort(x[_idx])
    yval = np.cumsum(np.ones_like(sortx)) / Nsim
    if np.max(yval)>1:
        hmag1.append(sortx[np.argmin(np.abs(yval-1.))])
    else:
        hmag1.append(np.nan)
    if np.max(yval)>3:
        hmag3.append(sortx[np.argmin(np.abs(yval-3.))])
    else:
        hmag3.append(np.nan)
    if np.max(yval)>10:
        hmag10.append(sortx[np.argmin(np.abs(yval-10.))])
    else:
        hmag10.append(np.nan)
    plt.yscale('log')
    plt.axvline(x=hmag1[-1])
    plt.axvline(x=hmag3[-1])
    plt.axvline(x=hmag10[-1])
    plt.plot(sortx, yval)
hmag1 = np.array(hmag1)
hmag3 = np.array(hmag3)
hmag10 = np.array(hmag10)

#%% add ticks for Teff
marr = np.linspace(0.1, 0.6, 100000)
tarr = m2t(marr)
mticks = []
#tticks = [2750, 3000, 3250, 3500, 3750]
tticks = [2800, 3100, 3250, 3400, 3600] # M1-5
for t in tticks:
    mticks.append(marr[np.argmin(np.abs(tarr-t))])

#%%
fig = plt.figure(figsize=(17,8))
ax = fig.gca()
ax.set_ylabel("stellar mass ($M_\odot$)")
ax.set_xlabel("$H_\mathrm{w}$ (mag)")
#plt.yscale('log')
#plt.xscale('log')
plt.xlim(5, 12)
plt.ylim(0.05, 0.6)

#"""
for i, (key, marker, size) in enumerate(zip(datakeys, ['^', 's', 'o', '*'], [8, 8, 8, 16])):
    _d = datadict[key]
    _idx = (_d.pl_rade < 1.5)
    plt.plot(_d["Hmag"][_idx], _d["mass"][_idx], '.', marker=marker, mfc='none', label=key, color='C%d'%i, markersize=size)
    _idx = (_d.pl_rade < 1.5) & (HZperiod(_d.mass, s=0.85)<_d.pl_orbper) & (_d.pl_orbper<HZperiod(_d.mass, s=0.25))
    plt.plot(_d["Hmag"][_idx], _d["mass"][_idx], '.', marker=marker, color='C%d'%i, markersize=size)
#"""

plt.plot(8, 0.1, '.', markersize=0, label='(filled: classical HZ)')

""" line
cont1 = np.poly1d(np.polyfit(hmag1[hmag1==hmag1], mvals[hmag1==hmag1], deg=3))
cont3 = np.poly1d(np.polyfit(hmag3[hmag3==hmag3], mvals[hmag3==hmag3], deg=3))
cont10 = np.poly1d(np.polyfit(hmag10[hmag10==hmag10], mvals[hmag10==hmag10], deg=3))

h0 = np.linspace(7, 10, 100)
plt.plot(h0, cont1(h0), ls='--', lw=1, color='k')

h0 = np.linspace(7.5, 10.5, 100)
plt.plot(h0, cont3(h0), ls='--', lw=1, color='k')

h0 = np.linspace(8.6, 11, 100)
plt.plot(h0, cont10(h0), ls='dotted', lw=1, color='k')
"""

# points
#plt.plot(hmag1, mvals, 'o', color='tan')
#plt.plot(hmag3, mvals, 'o', color='tan')
#plt.plot(hmag10, mvals, 'o', color='tan')

# step
plt.plot(np.sort(np.tile(hmag1, 2))[::-1], np.sort(np.r_[medges[1:], medges[0:-1]]), ls='dashed', lw=1, color='k')
#plt.plot(np.sort(np.tile(hmag3, 2))[::-1], np.sort(np.r_[medges[1:], medges[:-1]]), ls='dotted', lw=1, color='k')
plt.plot(np.sort(np.tile(hmag10, 2))[::-1], np.sort(np.r_[medges[1:], medges[0:-1]]), ls='dotted', lw=1, color='k')

#plt.text(7.5-0.5, 0.54, "brightest transiting HZ Earth\n(100\% occurrence)", ha='right', va='bottom', fontsize=16)
plt.text(7.5-0.3, 0.54, r"brightest hosts of transiting HZ Earths\\\phantom{aa}detectable w/ JASMINE (100\% occurrence)", ha='right', va='bottom', fontsize=16)
plt.text(7.5+1.7-0.5, 0.54+0.02, "(10\% occurrence)", ha='right', va='bottom', fontsize=16)

plt.title("transiting planets smaller than $1.5\,R_\oplus$")
plt.legend(loc='lower left')

ax2 = ax.twinx()
ax2.set_xscale('linear')#, subsx=[-2])
ax2.set_yticks(mticks)
ax2.set_yticklabels(['${0:d}$'.format(int(t)) for t in tticks])
ax2.set_ylabel("effective temperature (K)")
ax2.set_ylim(0.05, 0.6)

plt.savefig("brightest.png", dpi=200, bbox_inches="tight");
