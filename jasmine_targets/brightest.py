#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.constants import R_sun, R_earth
from astroquery.mast import Catalogs
from otofu.jasmine import noise_model, HZperiod

#%%
import seaborn as sns
sns.set(style='ticks', font_scale=1.6, font='sans-serif')
#sns.set(style='whitegrid', font_scale=1.6, font='sans-serif')
from matplotlib import rc
rc('text', usetex=True)

#%%
def assign_snHZ(d, period_key, mass_key, rad_key, teff_key, jmag_key, hmag_key):
    m, r, t, per = d[mass_key], d[rad_key], d[teff_key], d[period_key]
    d['Hwmag'] = d[jmag_key]*0.7 + d[hmag_key]*0.3
    try:
        d['pHZ'] = 365.25*(m**(-0.5))*(r**1.5)*((t/5777)**3)
    except:
        d['pHZ'] = HZperiod(d[mass_key])
    d['durHZ'] = (13./24.)*(d.pHZ/365.25)**(1./3.)*r/m**(1./3.)*np.pi/4.
    d['ptraHZ'] = 1/(3.7528*d.pHZ**(2./3.)*m**(1./3.)/r)
    d['depth_earth'] = r**(-2)*(R_earth/R_sun).value**2
    d['snHZ'] = d.depth_earth*1e6 / noise_model(d.Hwmag) * np.sqrt(d.durHZ*1440/5.)
    #d['snHZ_tess'] = d.depth_earth*1e6/tessnoise_5min(d.Tmag)*np.sqrt(d.durHZ*1440/5.)

#%%
d = pd.HDFStore('compiled_catalog.h5')
datakeys = [k[1:] for k in d.keys()]

#%%
# dtoinodup: TESS candidates (TOIs) - TESS confimed planets
# dconf_tess: TESS confirmed planets
# dconf_kepler: confirmed planets from Kepler/K2
# dconf_other: other confirmed planets
datadict = dict(zip(datakeys, [d[key] for key in datakeys]))
datalist = [d[key] for key in datakeys]

#%%
d = pd.read_csv('data/CTLv8_t4000_h11.csv', comment='#')
d = d[d.rad<0.6].reset_index(drop=True)
d['pHZ'] = HZperiod(d.mass)
_ = assign_snHZ(d, 'pHZ', 'mass', 'rad', 'Teff', 'Jmag', 'Hmag')

#%%
_idx = (d.Teff==d.Teff) & (d.mass==d.mass)
m2t = np.poly1d(np.polyfit(d.mass[_idx], d.Teff[_idx], deg=5))
t2m = np.poly1d(np.polyfit(d.Teff[_idx], d.mass[_idx],  deg=5))

#%%
marr = np.linspace(0.1, 0.6, 100000)
tarr = m2t(marr)
mticks = []
#tticks = [2750, 3000, 3250, 3500, 3750]
tticks = [2800, 3100, 3250, 3400, 3600] # M1-5
for t in tticks:
    mticks.append(marr[np.argmin(np.abs(tarr-t))])

#%%
def simulate(ds):
    rands = np.random.rand(len(ds))
    idx = rands < ds.ptraHZ
    return ds[idx].reset_index(drop=True)

#%%
Nsim = 5000
res = pd.DataFrame({})
for i in range(Nsim):
    dsim = simulate(d)
    res = res.append(dsim[dsim.snHZ>10][["Hmag", "mass"]])
    #res = res.append(dsim[["Hmag", "mass", "snHZ"]])

#%%
#x, y = np.array(res.Hmag[res.snHZ>10]), np.array(res.mass[res.snHZ>10])
x, y = np.array(res.Hmag), np.array(res.mass)
medges = np.arange(0.05, 0.55, 0.07)
mvals = 0.5 * (medges[1:] + medges[:-1])

#%%
hmag1, hmag3, hmag10 = [], [], []
for mu, ml in zip(medges[1:], medges[:-1]):
    _idx = (ml<=y) & (y<mu)
    plt.figure()
    #plt.hist(x[_idx], cumulative=True, weights=w[_idx])
    sortx = np.sort(x[_idx])
    yval = np.cumsum(np.ones_like(sortx))/Nsim
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
    plt.axvline(x=hmag10[-1])
    plt.plot(sortx, yval)
hmag1 = np.array(hmag1)
hmag3 = np.array(hmag3)
hmag10 = np.array(hmag10)

#%%
fig = plt.figure(figsize=(18,8))
ax = fig.gca()
ax.set_ylabel("stellar mass ($M_\odot$)")
ax.set_xlabel("$H_\mathrm{w}$ (mag)")
#plt.yscale('log')
#plt.xscale('log')
plt.xlim(5, 12)
plt.ylim(0.05, 0.6)

#"""
for i, (key, marker, size) in enumerate(zip(datakeys, ['^', 's', '*', 'o'], [8, 8, 16, 8])):
    _d = datadict[key]
    _idx = (_d.pl_rade < 1.5)
    plt.plot(_d["Hmag"][_idx], _d["mass"][_idx], '.', marker=marker, mfc='none', label=key, color='C%d'%i, markersize=size)
    _idx = (_d.pl_rade < 1.5) & (HZperiod(_d.mass, s=0.85)<_d.pl_orbper) & (_d.pl_orbper<HZperiod(_d.mass, s=0.25))
    plt.plot(_d["Hmag"][_idx], _d["mass"][_idx], '.', marker=marker, color='C%d'%i, markersize=size)
#"""

h0 = np.linspace(7, 10, 100)
#""" line
#plt.plot(hmag1, mvals, 'o', color='tan')
#plt.plot(hmag10, mvals, 'o', color='tan')
cont1 = np.poly1d(np.polyfit(hmag1, mvals, deg=3))
cont3 = np.poly1d(np.polyfit(hmag3[hmag3==hmag3], mvals[hmag3==hmag3], deg=3))
cont10 = np.poly1d(np.polyfit(hmag10[hmag10==hmag10], mvals[hmag10==hmag10], deg=3))
plt.plot(h0, cont1(h0), ls='--', lw=1, color='k')#, label='brightest tHZ Earth detectable with exo-Jasmine\n(100\% occurrence)')
#h0 = np.linspace(7.5, 10.5, 100)
#plt.plot(h0, cont3(h0), ls='--', lw=1, color='k')
h0 = np.linspace(8.6, 11, 100)
plt.plot(h0, cont10(h0), ls='dotted', lw=1, color='k')#, label='brightest tHZ Earth detectable with exo-Jasmine\n(10\% occurrence)')
#"""

# step
#plt.plot(np.sort(np.tile(hmag1, 2))[::-1], np.sort(np.r_[medges[1:], medges[:-1]]), ls='dashed', lw=1, color='k')
#plt.plot(np.sort(np.tile(hmag10, 2))[::-1], np.sort(np.r_[medges[1:], medges[:-1]]), ls='dotted', lw=1, color='k')

plt.text(7.5-0.5, 0.54, "brightest transiting HZ Earth\n(100\% occurrence)", ha='right', va='bottom', fontsize=16)
plt.text(7.5+1.7-0.5, 0.54+0.02, "(10\% occurrence)", ha='right', va='bottom', fontsize=16)

plt.title("known transiting planets smaller than $1.5\,R_\oplus$ (filled: classical HZ)")
plt.legend(loc='lower left')

ax2 = ax.twinx()
ax2.set_xscale('linear')#, subsx=[-2])
ax2.set_yticks(mticks)
ax2.set_yticklabels(['${0:d}$'.format(int(t)) for t in tticks])
ax2.set_ylabel("effective temperature (K)")
ax2.set_ylim(0.05, 0.6)

#plt.savefig("brightest_sim.png", dpi=200, bbox_inches="tight")
#plt.savefig("brightest_known.png", dpi=200, bbox_inches="tight")
plt.savefig("brightest_known_with_best.png", dpi=200, bbox_inches="tight");
