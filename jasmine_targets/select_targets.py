#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.constants import R_jup, R_sun, R_earth, au, M_sun, M_earth, M_jup
from astroquery.mast import Catalogs

#%%
import seaborn as sns
sns.set(style='ticks', font_scale=1.6, font='sans-serif')
#sns.set(style='whitegrid', font_scale=1.6, font='sans-serif')
from matplotlib import rc
rc('text', usetex=True)

#%% ejasmine noise in ppm as a function of Hw=0.7J+0.3H for 5min
def noise_model(hw, hwlim=[7.756708395177219, 12.756708395177219], coeff=[2.73659218e+00, -9.95974284e+01,  1.37688923e+03, -8.47717270e+03,1.95861269e+04]):
    m = np.poly1d(coeff)
    if np.sum((hw<hwlim[0])|(hwlim[1]<hw)):
        print ('Hw outside the valid range:', hwlim)
    return m(hw)

def HZperiod(mass, s=1, z=np.poly1d([ 0.40127856,  1.54449068,  4.0068301 , -1.17924394])):
    return 12*s**(-3./4.)/np.sqrt(mass/0.15)*(np.exp(z(np.log(mass)))/3e-3)**(3./4.)

def assign_snHZ(d, period_key, mass_key, rad_key, teff_key, jmag_key, hmag_key):
    m, r, t, per = d[mass_key], d[rad_key], d[teff_key], d[period_key]
    d['pHZ'] = 365.25*(m**(-0.5))*(r**1.5)*((t/5777)**3)
    #pHZ = HZperiod(d[mass_key])
    d['depth_earth'] = r**(-2)*(R_earth/R_sun).value**2
    d['durHZ'] = (13./24.)*(d.pHZ/365.25)**(1./3.)*r/m**(1./3.)*np.pi/4.
    d['ptraHZ'] = 1/(3.7528*d.pHZ**(2./3.)*m**(1./3.)/r)
    d['Hwmag'] = d[jmag_key]*0.7 + d[hmag_key]*0.3
    d['snHZ'] = d.depth_earth*1e6/noise_model(d.Hwmag)*np.sqrt(d.durHZ*1440/5.)
    #d['snHZ_tess'] = d.depth_earth*1e6/tessnoise_5min(d.Tmag)*np.sqrt(d.durHZ*1440/5.)

#%% confirmed planets
dconf = pd.read_csv('data/confirmed_20210427.csv', comment='#')
dconf['disc_telescope'].iloc[np.array(dconf.disc_telescope!=dconf.disc_telescope)] = "N/A"
dconf_tess = dconf[dconf.disc_telescope.str.contains("TESS")].reset_index(drop=True)
dconf_tra = dconf[(~dconf.disc_telescope.str.contains("TESS"))&(dconf.tran_flag==1.)].reset_index(drop=True)

#%% known TOIs
#dconf_tess = pd.read_csv("data/TESS_confirmed_20210427.csv", comment='#') this is a subset
dtoi = pd.read_csv('data/TOI_20210427.csv', comment='#')
print ('# %d stars with TOIs.'%len(dtoi.drop_duplicates("tid", keep='first')))

#%%
func = lambda x: int(x[4:]) if x==x else -1
dconf_tra['tid'] = list(map(func, dconf_tra.tic_id))
dconf_tess['tid'] = list(map(func, dconf_tess.tic_id))
print ('# %d stars with TIC IDs.'%len(dconf_tra[dconf_tra.tid>0].drop_duplicates('tid', keep='first')))

#%%
catalog_data = Catalogs.query_criteria(catalog="TIC", ID=dconf_tra.tid).to_pandas()
catalog_data['tid'] = np.array(catalog_data.ID).astype(int)
print ("# %d TIC stars found."%len(catalog_data))
dconf_tra = pd.merge(dconf_tra, catalog_data[["tid", "mass", "rad", "Teff", "Jmag", "Hmag"]], on='tid').reset_index(drop=True)

#%%
catalog_data = Catalogs.query_criteria(catalog="TIC", ID=dtoi.tid).to_pandas()
#catalog_data2 = Catalogs.query_criteria(catalog="ctl", ID=dtoi.tid).to_pandas()
catalog_data['tid'] = np.array(catalog_data.ID).astype(int)
#catalog_data2['tid'] = np.array(catalog_data2.ID).astype(int)
print ("# %d TIC stars found."%len(catalog_data))
dtoi = pd.merge(dtoi, catalog_data[["tid", "mass", "rad", "Teff", "Jmag", "Hmag"]], on='tid').reset_index(drop=True)

#%%
dconf_tra['rad'] = dconf_tra['rad'].fillna(dconf_tra['st_rad'])
dtoi['rad'] = dtoi['rad'].fillna(dtoi['st_rad'])

#%%
dconf_tra['pl_rade'] = dconf_tra['pl_rade'].fillna(dconf_tra.rad*R_sun*np.sqrt(dconf_tra.pl_trandep*1e-2)/R_earth)
dtoi['pl_rade'] = dtoi['pl_rade'].fillna(dtoi.rad*R_sun*np.sqrt(dtoi.pl_trandep*1e-6)/R_earth)

#%% missing Teff
dconf_tra.Teff[dconf_tra.hostname=="TRAPPIST-1"] = 2566

#%%
assign_snHZ(dtoi, 'pl_orbper', 'mass', 'rad', 'Teff', 'Jmag', 'Hmag')
assign_snHZ(dconf_tra, 'pl_orbper', 'mass', 'rad', 'Teff', 'Jmag', 'Hmag')
assign_snHZ(dconf_tess, 'pl_orbper', 'st_mass', 'st_rad', 'st_teff', 'sy_jmag', 'sy_hmag')
dconf_tess['mass'] = dconf_tess.st_mass
dconf_tess['rad'] = dconf_tess.st_rad

#%%
dtoi['pl_name'] = dtoi.toi
dtoi['hostname'] = ["TOI-%d"%dtoi.toipfx[i] for i in range(len(dtoi))]
dtoinodup = dtoi[(~dtoi.tid.isin(dconf_tess.tid))&(~dtoi.tid.isin(dconf_tra.tid))].reset_index(drop=True)
dconf_tess.mass[dconf_tess.tid==259377017] = 0.362274

#%%
kepidx = dconf_tra.disc_facility.str.contains("Kepler")+dconf_tra.disc_facility.str.contains("K2")
dconf_kepler = dconf_tra[kepidx].reset_index(drop=True)
dconf_other = dconf_tra[~kepidx].reset_index(drop=True)

#%%
m0 = np.linspace(0.03, 0.65, 100)
phz0 = HZperiod(m0, s=0.85)
phz0_upp = HZperiod(m0, s=0.25)

#dtoidup[dtoidup.toipfx==1468][['tid', 'toi', 'mass', 'pl_orbper']]
#dconf_tess[dconf_tess.hostname=="LTT 3780"][['tid', 'mass', 'pl_orbper']]

#%%
dconf_tess['Hmag'] = dconf_tess.sy_hmag

#%%
# dtoinodup: TESS candidates (TOIs) - TESS confimed planets
# dconf_tess: TESS confirmed planets
# dconf_kepler: confirmed planets from Kepler/K2
# dconf_other: other confirmed planets
datasets = dict(zip(['TESS candidate', 'TESS confirmed', "Kepler/K2", "ground"], [dtoinodup, dconf_tess, dconf_kepler, dconf_other]))

#%%
d = pd.read_csv('~/Dropbox/research_notes/data/CTLv8_t4000_h11.csv', comment='#')
d = d[d.rad<0.6].reset_index(drop=True)
d['pHZ'] = HZperiod(d.mass)
assign_snHZ(d, 'pHZ', 'mass', 'rad', 'Teff', 'Jmag', 'Hmag')

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
snthreshold = 10.
fig = plt.figure(figsize=(18,8))
ax = fig.gca()
ax.set_ylabel("stellar mass ($M_\odot$)")
ax.set_xlabel("orbital period (days)")
#plt.yscale('log')
plt.xscale('log')
plt.ylim(0.05, 0.57)
xkey, ykey = 'pl_orbper', 'pl_rade'
xkey, ykey = 'pl_orbper', 'mass'
ax.fill_betweenx(m0, phz0, phz0_upp, color='gray', alpha=0.2)
#ax.text(48, 0.48, "classical\nhabitable zone", color='k', zorder=1000, alpha=0.8)
ax.text(30, 0.3, "habitable zone", color='k', zorder=1000, alpha=0.8)
for df, mk, lab, ms in zip([dtoinodup, dconf_tess, dconf_kepler, dconf_other], ['^', 's', 'o', '*'], ['TESS candidate', 'TESS confirmed', "Kepler/K2", "ground"], [40, 40, 40, 80]):
    #_df = df[['tid', 'pl_rade']].groupby('tid', as_index=False).max()
    #df = pd.merge(df, _df.rename({"pl_rade": "maxrade"}, axis='columns'), on='tid')
    _idx = (df['snHZ'] > snthreshold) & (df.pl_rade < 1.5)
    plt.scatter(df[_idx][xkey], df[_idx][ykey], alpha=1, #facecolors='none',
    c=df[_idx].snHZ, vmin=10, vmax=40, cmap='coolwarm',
    edgecolors='k', lw=0.5, marker=mk, label=lab, s=ms)
    if mk=='^':
        plt.colorbar(pad=0.08, label='JASMINE S/N for habitable Earth')
    _idx = (df['snHZ'] > snthreshold) & (df.pl_rade >= 1.5)
    plt.scatter(df[_idx][xkey], df[_idx][ykey], alpha=0.3, #facecolors='none',
    c=df[_idx].snHZ, vmin=10, vmax=40, cmap='coolwarm',
    edgecolors='k', lw=0.5, marker=mk, s=ms)
    _idx = (df['snHZ'] > snthreshold)
    _tids = list(set(df[_idx].tid))
    for _tid in _tids:
        __idx = _idx & (df.tid==_tid)
        plt.plot(df[__idx][xkey], df[__idx][ykey], ls='dotted', lw=0.5, color='gray')
        #""" multi
        if (np.sum(__idx&(df.mass<0.5)) > 1):
            #print (df[__idx][['pl_name', 'hostname', 'mass', 'pl_orbper', 'pl_rade', 'pHZ', 'snHZ', "Hwmag"]])
            #ax.text(np.min(df[__idx]['pl_orbper'])*0.05+0.1, df[__idx]['mass'].iloc[0], df[__idx]['hostname'].iloc[0], color='k', zorder=1000, alpha=0.8, va='center', fontsize=10)
            ax.text(np.min(df[__idx]['pl_orbper'])*0.95, df[__idx]['mass'].iloc[0], df[__idx]['hostname'].iloc[0], zorder=1000, alpha=0.8, va='center', ha='right', fontsize=10)
        """
        if (np.sum(__idx&(df.mass<0.5)) == 1 and df[__idx].snHZ.iloc[0] > 20):
            ax.text(np.min(df[__idx]['pl_orbper'])*1.02, df[__idx]['mass'].iloc[0], df[__idx]['hostname'].iloc[0], zorder=1000, alpha=0.8, va='bottom', ha='left', fontsize=10)
        """
plt.legend(loc='lower right')
plt.xlim(0.2, 200)

ax2 = ax.twinx()
ax2.set_yticks(mticks)
ax2.set_yticklabels(['${0:d}$'.format(int(t)) for t in tticks])
ax2.set_ylabel("effective temperature (K)")
ax2.set_ylim(0.05, 0.55)
plt.savefig("select_targets_known_multi2.png", dpi=200, bbox_inches="tight")
#plt.savefig("select_targets_known_single2.png", dpi=200, bbox_inches="tight")

#%%
plt.plot(d.Teff, d.Tmag, '.')

#%%
d.mass[np.nanargmin(m2t(d.mass)-2750)]
plt.plot(d.mass, d.Teff, '.')
plt.plot(d.mass, m2t(d.mass), '.')

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
for i, (key, marker, size) in enumerate(zip(datasets.keys(), ['^', 's', 'o', '*'], [8, 8, 8, 16])):
    _d = datasets[key]
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
ax2.set_xscale('linear', subsx=[-2])
ax2.set_yticks(mticks)
ax2.set_yticklabels(['${0:d}$'.format(int(t)) for t in tticks])
ax2.set_ylabel("effective temperature (K)")
ax2.set_ylim(0.05, 0.6)

#plt.savefig("brightest_sim.png", dpi=200, bbox_inches="tight")
#plt.savefig("brightest_known.png", dpi=200, bbox_inches="tight")
plt.savefig("brightest_known_with_best.png", dpi=200, bbox_inches="tight")

#%% extended: https://figshare.com/articles/dataset/TESS_Extended_Mission_Yield_Simulations/11775081 (Barclay+2018)
dtesssim = pd.read_csv("data/detected_planets-CTL8.csv", comment='#')
dtesssim["mass"] = dtesssim["MASS"]
dtesssim["prad"] = dtesssim["planetRadius"]
dtesssim["period"] = dtesssim["planetPeriod"]
dtesssim.keys()

#%% extended: (Barclay+2018) original
dtesssim = pd.read_csv("data/barclay18.tsv", comment='#', delimiter='|')
dtesssim['period'] = dtesssim.Perp
dtesssim['prad'] = dtesssim.Rp
dtesssim['mass'] = dtesssim.Mass

#%% extended:
dtesssim = pd.read_csv("data/detected_planets.csv", comment='#')
dtesssim["mass"] = dtesssim["Star-mass"]
dtesssim["prad"] = dtesssim["Planet-radius"]
dtesssim["period"] = dtesssim["Planet-period"]
dtesssim.keys()

#%%
snthreshold = 10.
fig = plt.figure(figsize=(18,8))
ax = fig.gca()
ax.set_ylabel("stellar mass ($M_\odot$)")
ax.set_xlabel("orbital period (days)")
#plt.yscale('log')
plt.xscale('log')
plt.ylim(0.05, 0.57)
xkey, ykey = 'pl_orbper', 'pl_rade'
xkey, ykey = 'pl_orbper', 'mass'
"""
simidx = dtesssim.prad<1.75
ax.plot(dtesssim.period[simidx], dtesssim.mass[simidx], 'o', color='C0', label='$<1.5$ Earth radii')
ax.plot(dtesssim.period[~simidx], dtesssim.mass[~simidx], '.', alpha=0.2, color='C0', label='$>1.5$ Earth radii')
"""
plt.scatter(dtesssim.period, dtesssim.mass, c=dtesssim.prad, vmin=1, vmax=1.5, s=12, alpha=1)
plt.colorbar(pad=0.08, label='planet radius (Earth radii)')
ax.fill_betweenx(m0, phz0, phz0_upp, color='gray', alpha=0.2)
#ax.text(48, 0.48, "classical\nhabitable zone", color='k', zorder=1000, alpha=0.8)
ax.text(30, 0.3, "classical\nhabitable zone", color='k', zorder=1000, alpha=0.8)
#plt.legend(loc='lower right')
plt.xlim(0.2, 200)

ax2 = ax.twinx()
ax2.set_yticks(mticks)
ax2.set_yticklabels(['${0:d}$'.format(int(t)) for t in tticks])
ax2.set_ylabel("effective temperature (K)")
ax2.set_ylim(0.05, 0.55)
plt.title("TESS Yield Simulations from Barclay et al. (2018)")
plt.savefig("tess_barclay2018.png", dpi=200, bbox_inches="tight")
#plt.savefig("select_targets_known_single2.png", dpi=200, bbox_inches="tight")

#%%
"""
fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 1, figsize=(14*1.3,10*1.3), sharex=True)
ax5.set_xlabel("stellar mass ($M_\odot$)")
ax5.set_xlim(0.1, 0.6)
ax1.plot(d.mass, d.pHZ, '.')
ax1.set_ylabel("orbital period (days)")
ax1.set_yscale("log")
ax2.plot(d.mass, d.durHZ*24., '.')
ax2.set_ylabel("transit duration (hr)")
#ax2.set_yscale("log")
ax3.plot(d.mass, d.ptraHZ*100, '.')
ax3.set_ylabel("transit probability (\%)")
ax4.plot(d.mass, d.depth_earth*1e2, '.')
ax4.set_ylabel("transit depth (\%)")
ax4.set_yscale("log")
ax5.plot(d.mass, d.snHZ, '.')
ax5.set_ylabel("transit S/N")
ax5.set_yscale("log")
ax5.axhline(y=10, color='gray')
ax1.yaxis.set_label_coords(-0.05, 0.5)
ax2.yaxis.set_label_coords(-0.05, 0.5)
ax3.yaxis.set_label_coords(-0.05, 0.5)
ax4.yaxis.set_label_coords(-0.05, 0.5)
ax5.yaxis.set_label_coords(-0.05, 0.5)
#plt.savefig("HZearth_mass.png", dpi=200, bbox_inches="tight")

#%%　ds: searchable by ejasmine
idxs = (d.snHZ>10) & (d.pHZ<30)
ds = d[idxs].reset_index(drop=True)
dall = d

#%%
dknown[dknown.st_teff!=dknown.st_teff]

#%%
#_dtoi = pd.read_csv('data/TOI_20200421.csv', comment='#')
#_dtoi = pd.read_csv('data/TOI_20201115.csv', comment='#')
#_dtoi = pd.read_csv('data/toi+_20210427.csv', comment='#')
#dtoi = pd.merge(_dtoi, d, left_on='TIC', right_on='ID', how='left') # supercedes TOI catalog values; some still miss radius etc
#dtoi.rad[dtoi.rad!=dtoi.rad] = dtoi['Star Radius Value'][dtoi.rad!=dtoi.rad]
_dtoi = pd.read_csv('data/TOI_20210427.csv', comment='#')
dtoi = pd.merge(_dtoi, d, left_on='tid', right_on='ID', how='left')
dtoi.rad[dtoi.rad!=dtoi.rad] = dtoi['st_rad'][dtoi.rad!=dtoi.rad]
dtoi['prad'] = dtoi.rad*R_sun*np.sqrt(dtoi.pl_trandep*1e-6)/R_earth
dtoi['idxHZe'] = (0.5<dtoi.pl_orbper/dtoi.pHZ)&(dtoi.pl_orbper/dtoi.pHZ<2)&(dtoi.prad<radmax)
nhze = dtoi[['tid', 'idxHZe']].groupby('tid', as_index=False).sum()
dtoi = pd.merge(dtoi, nhze.rename(columns={'idxHZe':'nHZe'}), on='tid', how='left')


#%%
#dconf = pd.read_csv('~/Dropbox/research_notes/data/confirmed_20200326.csv', comment='#')
dconf = pd.read_csv('data/confirmed_20210427.csv', comment='#')
dconf['plname'] = dconf.hostname+' '+dconf.pl_letter
dconftra = dconf[dconf.tran_flag==1]
dconftra['Hw'] = dconftra.sy_jmag*0.7+dconftra.sy_hmag*0.3
ids = []
for _tic in dconftra.tic_id:
    try:
        ids.append(int(_tic[4:]))
    except:
        ids.append(-1)
dconftra["ID"] = ids

#%%
dconftra = pd.merge(dconftra, d[['ID', 'rad', 'Teff', 'mass']], on='ID', how='left')
dconftra.st_teff[dconftra.st_teff!=dconftra.st_teff] = dconftra.Teff[dconftra.st_teff!=dconftra.st_teff]
dconftra.st_mass[dconftra.st_mass!=dconftra.st_mass] = dconftra.mass[dconftra.st_mass!=dconftra.st_mass]
#dconftra.st_teff[dconftra.hostname=='GJ 436']# = 3318
dconftra.st_teff[dconftra.pl_name=='TRAPPIST-1 h'] = 2559
dconftra.st_mass[dconftra.pl_name=='TRAPPIST-1 h'] = 0.08

#%%
radmax = 1.5
dknown = dconftra[(dconftra.st_mass<0.5)&(dconftra.Hw<11.5)].reset_index(drop=True)
dknown['pHZ'] = 365.25*(dknown.st_mass**(-0.5))*(dknown.st_rad**1.5)*((dknown.st_teff/5777)**3)
dknown['idxHZe'] = (0.5<dknown.pl_orbper/dknown.pHZ)&(dknown.pl_orbper/dknown.pHZ<2)&(dknown.pl_rade<radmax)
nhze = dknown[['hostname', 'idxHZe']].groupby('hostname', as_index=False).sum()
dknown = pd.merge(dknown, nhze.rename(columns={'idxHZe':'nHZe'}), on='hostname', how='left')

#%%
dconftra[(dconftra.st_teff!=dconftra.st_teff)&(dconftra.pl_rade<2)&(dconftra.st_mass<0.6)][['sy_hmag']]

#%%
"""
_p, _m = d.pHZ, d.mass
idx = (_p==_p) & (_m==_m)
_p, _m = _p[idx], _m[idx]
z = np.poly1d(np.polyfit(np.log(_p), np.log(_m), deg=3))
p0 = np.logspace(0.5, 2, 100)
plt.xscale("log")
plt.plot(_p, _m, '.')
plt.plot(p0, np.exp(z(np.log(p0))))
"""

#%%
#dtess = pd.read_csv("data/TESS_confirmed_20201121.csv", comment='#')
dtess = pd.read_csv("data/TESS_confirmed_20210427.csv", comment='#')

#%%
from sklearn.neighbors import KDTree
X = np.array(dconftra[['ra', 'dec']])
tree = KDTree(X)
dist, ind = tree.query(np.array(dtoi[['ra_x', 'dec_x']]), k=1)
dist = dist.reshape(1,-1)[0]
dtoi["dist"] = dist # toi to confirmed transits

X = np.array(dtess[['ra', 'dec']])
tree = KDTree(X)
dist, ind = tree.query(np.array(dtoi[['ra_x', 'dec_x']]), k=1)
dist = dist.reshape(1,-1)[0]
dtoi["dist2"] = dist # toi to tess

X = np.array(dconftra[['ra', 'dec']])
tree = KDTree(X)
dist, ind = tree.query(np.array(dtess[['ra', 'dec']]), k=1)
dist = dist.reshape(1,-1)[0]
dtess["dist"] = dist # tess conf to other conf


"""

#%% HZ period from stellar mass
# based on empirical luminosity - mass relation
def pHZ_mass(m, s=1):
    logm, logl = np.log(d.mass), np.log(d.lum)
    idx = (logm==logm)&(logl==logl)&(d.mass<0.4)
    logm, logl = logm[idx], logl[idx]
    z = np.poly1d(np.polyfit(logm, logl, deg=3))
    return 12*s**(-3./4.)/np.sqrt(m/0.15)*(np.exp(z(np.log(m)))/3e-3)**(3./4.)

dconftra[dconftra['hostname']=="TRAPPIST-1"][['st_mass', 'pl_orbper', 'pl_rade', 'sy_dist', 'st_mass']]

#%%
rmax = 1.5#+0.5
m0 = np.linspace(0.03, 0.65, 100)
phz0 = pHZ_mass(m0, s=0.85)
phz0_upp = pHZ_mass(m0, s=0.25)
#phz0_low = pHZ_mass(m0, s=1)
fig = plt.figure(figsize=(16,7))
ax = fig.gca()
ax.set_ylabel("stellar mass ($M_\odot$)")
ax.set_xlabel("orbital period (days)")
ax.set_xscale("log")
ax.set_ylim(0.05, 0.55)
ax.set_xlim(0.25, 250)
#ax.set_yscale("log")
ax.fill_betweenx(m0, phz0, phz0_upp, color='gray', alpha=0.2)
#ax.text(48, 0.48, "classical\nhabitable zone", color='k', zorder=1000, alpha=0.8)
ax.text(30, 0.3, "classical\nhabitable zone", color='k', zorder=1000, alpha=0.8)
idxe = (dconftra.pl_rade<rmax)
idxe &= (~dconftra.disc_telescope.str.contains("TESS"))
#idxe &= (dconftra.sy_dist<50)
#idxe &= (dconftra.sy_hmag<11)
ax.plot(dconftra.pl_orbper[idxe], dconftra.st_mass[idxe], '*', #mfc='none',
    mew=1, markersize=10, color='C1', label="ground/Kepler/K2")
#idxe = (dconftra.pl_rade<1.5)&(dconftra.st_j>11.5)
#ax.plot(dconftra.pl_orbper[idxe], dconftra.st_mass[idxe], '.', #mfc='none',
#    mew=1, markersize=4, color='gray', label="")
idxetoi = (dtess.pl_rade<rmax) #& (dtess.dist>1)
ax.plot(dtess.pl_orbper[idxetoi], dtess.st_mass[idxetoi], 's', markersize=6, mew=1., label='TESS confirmed', color="C0")
idxetoi = (dtoi.prad<rmax) & (dtoi.dist2>1) & (dtoi.tfopwg_disp!="FP") & (dtoi.dist2>1)
ax.plot(dtoi.pl_orbper[idxetoi], dtoi.mass[idxetoi], 's', markersize=6, mfc='none', mew=1., label='TESS candidate', color="C0")
#idxsim = dtesssim.prad<1.5
#ax.plot(dtesssim.period[idxsim], dtesssim.mass[idxsim], '.', color='C3', label='TESS simulation (Barclay+2018)')
#ax.legend(loc='lower right')
ax.legend(loc='upper right')
ax.set_title("transiting planets smaller than 1.5 Earth radii within 50pc")
#plt.savefig("mass_period2.png", dpi=200, bbox_inches="tight")
#plt.savefig("mass_period2_tess.png", dpi=200, bbox_inches="tight")
#set(dconftra[idxe&(dconftra.st_mass<0.6)].pl_hostname)
#for n,x,y in zip(dknown.plname[idxe], dknown.pl_orbper[idxe], dknown.st_mass[idxe]):
#    plt.text(x,y,n)

#%%
dkoi = pd.read_csv("~/Dropbox/research_notes/data/koi_20170927.csv")
#dkoi = dkoi[dkoi["koi_smass"]<0.6].reset_index(drop=True)
dkoi = dkoi[(dkoi.koi_prad>0.5)&(dkoi.koi_prad<4)&(dkoi.koi_smass<0.6)]
dkoi["hzin"] = 20#pHZ_mass(dkoi.koi_smass, s=0.85)
dkoi["hzout"] = 50#pHZ_mass(dkoi.koi_smass, s=0.25)

#%%
disps = ["CONFIRMED", "CANDIDATE"]
disps = ["CONFIRMED"]
dkoi["10d"] = (dkoi.koi_disposition.isin(disps)) & (dkoi.koi_period<7.)
dkoi = pd.merge(dkoi, dkoi.groupby("kepid", as_index=False).sum()[["kepid", "10d"]], on='kepid')

#%%
dkoi["10d_hz"] = dkoi.koi_disposition.isin(disps) & (dkoi["10d_y"]>0) & (dkoi.hzin < dkoi.koi_period) & (dkoi.koi_period < dkoi.hzout)
dkoi = pd.merge(dkoi, dkoi.groupby("kepid", as_index=False).sum()[["kepid", "10d_hz"]], on='kepid')

#%%
koistars = dkoi.drop_duplicates("kepid", keep='first')

#%%
np.sum(koistars["10d_hz_y"]>0)/np.sum(koistars["10d_y"]>0)

#%%
plt.yscale("log")
plt.xscale("log")
plt.plot(dkoi.koi_period[dkoi["10d_x"]], dkoi.koi_prad[dkoi["10d_x"]], '.')
plt.plot(dkoi.koi_period[dkoi["10d_hz_x"]], dkoi.koi_prad[dkoi["10d_hz_x"]], '.')

#%%
list(dconftra.keys())
dconftra[idxe][["pl_hostname", "pl_facility", "pl_telescope", "pl_instrument"]]

#%%
dtess_b18 = pd.read_csv("data/v2.0_primary_1_trial_0.csv")
dtess_b18["m1"] = np.poly1d(np.array([ 1.84588214, -1.83692212,  1.47262963, -0.00589032]))(dtess_b18["r1"])

#%%
gcol, scol = "C1", "C0"
gmk, smk = "*", "s"
rmax = 1.5#+0.5*2
m0 = np.linspace(0.03, 0.65, 100)
phz0 = pHZ_mass(m0, s=0.85)
phz0_upp = pHZ_mass(m0, s=0.25)
#phz0_low = pHZ_mass(m0, s=1)
fig = plt.figure(figsize=(16,7))
ax = fig.gca()
ax.set_ylabel("stellar mass ($M_\odot$)")
ax.set_xlabel("orbital period (days)")
ax.set_xscale("log")
ax.set_ylim(0.05, 0.55)
ax.set_xlim(0.25, 250)
#ax.set_yscale("log")
ax.fill_betweenx(m0, phz0, phz0_upp, color='gray', alpha=0.2)
#ax.text(48, 0.48, "classical\nhabitable zone", color='k', zorder=1000, alpha=0.8)
ax.text(30, 0.3, "classical\nhabitable zone", color='k', zorder=1000, alpha=0.8)
idxe = (dconftra.pl_rade<rmax)
idxe &= (dconftra.sy_dist<50)
#idxe &= (dconftra.st_h<12)
idxs = dconftra.pl_facility.str.contains("Kepler")+dconftra.pl_facility.str.contains("K2")+dconftra.pl_facility.str.contains("TESS")+dconftra.pl_facility.str.contains("CoRoT")
idxeg = idxe & (idxs==0)
idxes = idxe & idxs
ax.plot(dconftra.pl_orbper[idxeg], dconftra.st_mass[idxeg], gmk,
    mew=1, markersize=10, color=gcol, label="ground")
ax.plot(dconftra.pl_orbper[idxes], dconftra.st_mass[idxes], smk,
    mew=1, markersize=6, color=scol, label="space")
#idxe = (dconftra.pl_rade<1.5)&(dconftra.st_j>11.5)
#ax.plot(dconftra.pl_orbper[idxe], dconftra.st_mass[idxe], '.', #mfc='none',
#    mew=1, markersize=4, color='gray', label="")
idxetoi = (dtess.pl_rade<rmax)&(dtess.dist>1)
ax.plot(dtess.pl_orbper[idxetoi], dtess.st_mass[idxetoi], 's', markersize=6, mew=1., label='', color="C0")
"""
idxetoi = (dtoi.prad<rmax)&(dtoi.dist>1)&(dtoi.tfopwg_disp!="FP")&(dtoi.dist2>1)
ax.plot(dtoi.pl_orbper[idxetoi], dtoi.mass[idxetoi], 's', markersize=6, mfc='none', mew=0.5, label='TESS candidate', color="C0")
"""
#idxsim = dtesssim.prad<rmax
#ax.plot(dtesssim.period[idxsim], dtesssim.mass[idxsim], '.', color='C3', label='TESS simulation (Barclay+2018)')
#ax.legend(loc='lower right')
idxb18 = dtess_b18.r2<rmax
ax.plot(dtess_b18.p[idxb18], dtess_b18.m1[idxb18], '.', color="C3", label='TESS simulation (Bouma+2017)')
ax.legend(loc='upper right')
ax.set_title("transiting planets smaller than 1.5 Earth radii within 50pc")
plt.savefig("mass_period_gs_tess2.png", dpi=200, bbox_inches="tight")
#set(dconftra[idxe&(dconftra.st_mass<0.6)].pl_hostname)
#for n,x,y in zip(dknown.plname[idxe], dknown.pl_orbper[idxe], dknown.st_mass[idxe]):
#    plt.text(x,y,n)

#%%
gcol, scol = "C1", "C0"
gmk, smk = "*", "s"
rmax = 1.5#+0.5*2
m0 = np.linspace(0.03, 0.65, 100)
phz0 = pHZ_mass(m0, s=0.85)
phz0_upp = pHZ_mass(m0, s=0.25)
#phz0_low = pHZ_mass(m0, s=1)
fig = plt.figure(figsize=(16,7))
ax = fig.gca()
ax.set_ylabel("stellar mass ($M_\odot$)")
ax.set_xlabel("orbital period (days)")
ax.set_xscale("log")
ax.set_ylim(0.05, 0.55)
ax.set_xlim(0.25, 250)
#ax.set_yscale("log")
ax.fill_betweenx(m0, phz0, phz0_upp, color='gray', alpha=0.2)
#ax.text(48, 0.48, "classical\nhabitable zone", color='k', zorder=1000, alpha=0.8)
ax.text(28, 0.3, "habitable zone", color='k', zorder=1000, alpha=0.8)
idxe = (dconftra.pl_rade<rmax)
idxe &= (dconftra.st_dist<50)
#idxe &= (dconftra.st_h<12)
idxs = dconftra.pl_facility.str.contains("Kepler")+dconftra.pl_facility.str.contains("K2")+dconftra.pl_facility.str.contains("TESS")+dconftra.pl_facility.str.contains("CoRoT")
idxeg = idxe & (idxs==0)
idxes = idxe & idxs
ax.plot(dconftra.pl_orbper[idxeg], dconftra.st_mass[idxeg], gmk,
    mew=1, markersize=10, color=gcol, label="ground")
ax.plot(dconftra.pl_orbper[idxes], dconftra.st_mass[idxes], smk,
    mew=1, markersize=6, color=scol, label="space")
idxetoi = (dtess.pl_rade<rmax)&(dtess.dist>1)
ax.plot(dtess.pl_orbper[idxetoi], dtess.st_mass[idxetoi], smk, markersize=6, mew=1., label='', color=scol)
ax.legend(loc='upper right')
ax.set_title("transiting planets smaller than 1.5 Earth radii within 50pc")
plt.savefig("mass_period_simple.png", dpi=200, bbox_inches="tight")
#set(dconftra[idxe&(dconftra.st_mass<0.6)].pl_hostname)
#for n,x,y in zip(dknown.plname[idxe], dknown.pl_orbper[idxe], dknown.st_mass[idxe]):
#    plt.text(x,y,n)

#%%
db19 = pd.read_csv("data/Ballard19.txt", delim_whitespace=True, comment='#')
db19

#%%
plt.xscale("log")
idxt2 = (db19.detected==1)&(db19.radius<1.5)
stars_known = db19[idxt2].drop_duplicates("starid", keep='first')
print (len(stars_known))
plt.plot(db19.period[idxt2], db19.teff[idxt2], '.')
idxt2 = (db19.detected==0)&(db19.transit==1)&(db19.radius<1.5)&(db19.period>10.)
stars_potential = db19[idxt2].drop_duplicates("starid", keep='first')
print (len(stars_potential))
print (np.sum(stars_potential.starid.isin(stars_known.starid)))
plt.plot(db19.period[idxt2], db19.teff[idxt2], '.')

#%%
stars_known
stars_potential

#%%
dtoi.tfopwg_disp
list(dtoi.keys())
dtoi[idxetoi&(dtoi.mass<0.6)&(dtoi.tfopwg_disp!="FP")]#[['pl_orbper']]

#%%
bins = np.linspace(0, 15, 100)
plt.hist(d.snHZ/d.snHZ_tess, bins=bins);

#%%
plt.hist(d.Jmag-d.Hmag)

#%%
list(d.keys())
plt.xscale("log")
plt.plot(1e3/d.plx[(d.mass<0.6)&(d.rad<0.6)], d.Hmag[(d.mass<0.6)&(d.rad<0.6)], '.')

#%%
p_obs = 365.
searchables = d[d.snHZ>10]
searchables_tess = d[d.snHZ_tess*np.sqrt(p_obs/d.pHZ)>10]

#%%
b = np.linspace(0.05, 0.7, 500)
plt.figure(figsize=(14,6))
plt.yscale("log")
plt.ylabel("cumulative count")
plt.xlabel("stellar mass ($M_\odot$)")
plt.xlim(0.06, 0.6)
plt.ylim(0.5, 5e4)
plt.hist(searchables.mass, cumulative=True, bins=b, histtype='step', lw=1.5)
plt.hist(searchables.mass, cumulative=True, bins=b, histtype='step', lw=1.5, weights=searchables.ptraHZ);

#%%
#b = np.linspace(0.05, 0.7, 500)
b = np.logspace(0, 3, 500)
plt.figure(figsize=(14,6))
plt.yscale("log")
plt.ylabel("cumulative count")
plt.xlabel("distance (pc)")
#plt.xlim(0.06, 0.6)
plt.xscale("log")
plt.ylim(0.5, 50)
plt.xlim(10, 100)
for m0 in [0.05, 0.15, 0.25, 0.35, 0.45]:
    idx = (m0<searchables.mass)&(searchables.mass<m0+0.1)
    plt.hist(1e3/searchables.plx[idx], cumulative=True, bins=b, histtype='step', lw=1.5, weights=searchables.ptraHZ[idx]*0.3, label='%.2f-%.2f'%(m0, m0+0.1))
plt.legend()
#plt.hist(searchables.mass, cumulative=True, bins=b, histtype='step', lw=1.5, weights=searchables.ptraHZ);

#%%
dknown[["pl_hostname", "st_mass"]]
dtess[dtess.pl_hostname=="TOI-700"][["pl_hostname", "st_mass"]]

#%% 20% occurrence
focc = 0.3
b = np.arange(0.075, 0.7, 0.05)
plt.figure(figsize=(14,5))
plt.yscale("log")
plt.ylabel("number in each bin")
plt.xlabel("stellar mass ($M_\odot$)")
plt.xlim(0.05, 0.65)
plt.ylim(0.3, 30)
#plt.hist(d.mass, bins=b, histtype='step', lw=1.5, weights=d.ptraHZ, ls='solid', color='gray')
plt.hist(searchables.mass, bins=b, histtype='step', lw=2, weights=searchables.ptraHZ*focc, ls='solid',
        label='%d%s occurrence'%(focc*100, "\%"))
plt.hist(searchables_tess.mass, bins=b, histtype='step', lw=1.5, weights=searchables_tess.ptraHZ*focc, ls='dotted',
        label='TESS 1yr, %d%s occurrence'%(focc*100, "\%"))
#for x in [0.08, 0.41]:
#    plt.axvline(x=x, color='gray')
plt.legend(loc='best')
plt.title("HZ Earth-sized planets searchable with exo-JASMINE")
plt.savefig("yield_comparison.png", dpi=200, bbox_inches="tight")
#plt.savefig("yield_jasmine_known.png", dpi=200, bbox_inches="tight")
#plt.savefig("yield_jasmine.png", dpi=200, bbox_inches="tight")

#%%
b = np.arange(0.075, 0.7, 0.05)
plt.figure(figsize=(14,5))
plt.yscale("log")
plt.ylabel("number in each bin")
plt.xlabel("stellar mass ($M_\odot$)")
plt.xlim(0.05, 0.65)
plt.ylim(0.05, 500)
#plt.hist(d.mass, bins=b, histtype='step', lw=1.5, weights=d.ptraHZ, ls='solid', color='gray')
plt.hist(searchables.mass, bins=b, histtype='step', lw=1.5, weights=searchables.ptraHZ, ls='solid',
        label='100\% occurrence')
plt.hist(searchables.mass, bins=b, histtype='step', lw=1.5, weights=searchables.ptraHZ*0.1, ls='dashed',
        label='10\% occurrence')
plt.hist(searchables_tess.mass, bins=b, histtype='step', lw=1.5, weights=searchables_tess.ptraHZ*0.1, ls='dotted',
        label='TESS 1yr, 10\% occurrence')
plt.legend(loc='best')
plt.title("HZ Earth-sized planets searchable with exo-JASMINE")
#plt.savefig("yield_comparison.png", dpi=200, bbox_inches="tight")
#plt.savefig("yield_jasmine.png", dpi=200, bbox_inches="tight")

#%%
from astroquery.mast import Catalogs
catalog_data = Catalogs.query_criteria(catalog="TIC", Hmag=[11.,12.], Teff=[2000,4000], rad=[0,0.6])
