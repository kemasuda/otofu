#%%
import astropy.units as u
from otofu.obsplan import check_field

#%%
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style='ticks', font_scale=1.6*0.5, font='sans-serif')
plt.rcParams["figure.figsize"] = (14*0.5,7*0.5)
from matplotlib import rc
rc('text', usetex=False)

#%%
fov = 0.45*u.deg

#%%
names = ["TOI-700", "LTT 1445A", "TRAPPIST-1"]

#%%
for name in names[:1]:
    coord, d = check_field(name, fov, plot=True)
