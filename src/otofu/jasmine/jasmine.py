import numpy as np

#%% jasmine noise model in ppm as a function of Hw=0.7J+0.3H for 5min
def noise_model(hw, hwlim=[7.756708395177219, 12.756708395177219], coeff=[2.73659218e+00, -9.95974284e+01,  1.37688923e+03, -8.47717270e+03,1.95861269e+04]):
    m = np.poly1d(coeff)
    if np.sum((hw<hwlim[0])|(hwlim[1]<hw)):
        print ('Hw outside the valid range:', hwlim)
    return m(hw)

#%% HZ period (classical HZ: s=0.25-0.85)
def HZperiod(mass, radius=None, teff=None, s=1):
    if radius is None or teff is None:
        return HZperiod_mass(mass, s=s)
    else:
        return HZperiod_mrt(mass, radius, teff, s=s)

#%% Lbol derived from empirical Lbol-mass relation (0.1<mass<0.4)
def HZperiod_mass(mass, s=1, z=np.poly1d([ 0.40127856,  1.54449068,  4.0068301 , -1.17924394])):
    return 12 / np.sqrt(mass/0.15) * (np.exp(z(np.log(mass)))/3e-3)**(3./4.) * s**(-3./4.)

#%% Lbol calculated from radius and temperature
def HZperiod_mrt(m, r, t, s=1):
    return 365.25 * (m**(-0.5)) * (r**1.5) * ((t/5777)**3) * s**(-3./4.)

#%% teff from mass for M dwarfs with 0.1<mass<0.6
def m2t(m, z=np.poly1d([203300.24449536, -402362.46252036, 317996.63455074, -123679.76418643, 24608.54436212,    1264.4338542])):
    return z(m)

#%% mass-teff ticks for 0.1<mass<0.5
def get_mtticks(tticks=[2800, 3100, 3250, 3400, 3600]):
    marr = np.linspace(0.1, 0.6, 100000)
    tarr = m2t(marr)
    mticks = [marr[np.argmin(np.abs(tarr-t))] for t in tticks]
    return mticks, tticks

#%% JASMINE S/N for HZ Earth
from astropy.constants import R_sun, R_earth
def assign_snHZ(d, mass_key, rad_key, teff_key, jmag_key, hmag_key, duty_cycle=0.5):
    m, r, t = d[mass_key], d[rad_key], d[teff_key]

    d['Hwmag'] = d[jmag_key]*0.7 + d[hmag_key]*0.3
    d['pHZ'] = HZperiod(m, radius=r, teff=t, s=0.85) # inner edge
    d['durHZ'] = (13./24.) * (d.pHZ/365.25)**(1./3.) * r/m**(1./3.) * np.pi/4.
    d['ptraHZ'] = 1 / (3.7528*d.pHZ**(2./3.)*m**(1./3.)/r)
    d['depth_earth'] = r**(-2) * (R_earth/R_sun).value**2
    d['snHZ'] = d.depth_earth*1e6 / noise_model(d.Hwmag) * np.sqrt(d.durHZ*1440/5.*duty_cycle)
    #d['snHZ_tess'] = d.depth_earth*1e6/tessnoise_5min(d.Tmag)*np.sqrt(d.durHZ*1440/5.)

    return d

"""
#%% Lbol-mass relation for 0.1 < mass < 0.4
import pandas as pd
import matplotlib.pyplot as plt
d = pd.read_csv('../../../jasmine_targets/data/CTLv8_t4000_h11.csv', comment='#')
d = d[d.rad<0.6].reset_index(drop=True)

#%%
logm, logl = np.log(d.mass), np.log(d.lum)
idx = (logm==logm) & (logl==logl) & (d.mass<0.4)
logm, logl = logm[idx], logl[idx]
z = np.poly1d(np.polyfit(logm, logl, deg=3))

#%%
logmarr = np.linspace(-2.4, -0.9, 100)
plt.plot(logm, logl, '.')
plt.plot(logmarr, z(logmarr), '-')

#%%
m = 0.2
print (HZperiod(m, s=0.85))
print (HZperiod(m, radius=0.22, teff=3200, s=0.85))

#%%
linm, lint = d.mass, d.Teff
logm, logl = np.log(d.mass), np.log(d.lum)
idx = (linm==linm) & (lint==lint)
linm, lint = linm[idx], lint[idx]

#%%
m2t = np.poly1d(np.polyfit(linm, lint, deg=5))
t2m = np.poly1d(np.polyfit(lint, linm,  deg=5))

#%%
plt.plot(linm, lint, '.')
plt.plot(linm, m2t(linm), '.');

#%%
plt.plot(lint, linm, '.')
plt.plot(lint, t2m(lint), '.');

#%%
z = np.poly1d(np.polyfit(logm, logl, deg=3))
"""
