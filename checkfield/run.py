#%%
import astropy.units as u
from checkfield import check_field

#%%
fov = 0.45*u.deg

#%%
names = ["TOI-700", "LTT 1445A", "TRAPPIST-1"]

#%%
for name in names:
    coord, d = check_field(name, fov, plot=True)
