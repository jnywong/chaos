import chaosmagpy
import shtns
import pyshtools as pysh
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from paropy.plot_utils import rad_to_deg, get_Z_lim
import scipy.special as sp
import numpy as np

lmax = 20
nphi = 480
ntheta = 200
fig_aspect = 0.8
n_levels = 30

def get_grid(nphi=256, ntheta=128):

    phi = np.linspace(0., 2*np.pi, nphi)
    x, w = sp.roots_legendre(ntheta)
    theta = np.sort(np.arccos(x))

    p2D = np.zeros([nphi, ntheta])
    th2D = np.zeros([nphi, ntheta])

    for i in range(nphi):
        p2D[i, :] = phi[i]

    for j in range(ntheta):
        th2D[:, j] = theta[j]

    return p2D, th2D

# Load CHAOS data
model = chaosmagpy.load_CHAOS_matfile('CHAOS-7.7.mat')
time = 2020.
timeCHAOS = chaosmagpy.data_utils.dyear_to_mjd(time)

# Gauss, max. degree 20
# g10, g11, h11, g20, g21, h21, g22, h22
ylm = model.synth_coeffs_tdep(timeCHAOS)
dtglm = model.synth_coeffs_tdep(timeCHAOS,deriv=1) # Secular Variation

mmax = lmax
sh = shtns.sht(lmax, mmax)
nlat, nphi = sh.set_grid(nlat = ntheta, nphi = nphi)

# Construct complex Gauss coefficients
k=0
j=0
m=0
glm = np.zeros(int(lmax*(lmax+3)/2))
hlm = np.zeros(int(lmax*(lmax+3)/2))
for l in range(1,lmax+1):
    for i in range(2*l+1):
        if i==0:            
            glm[j] = ylm[k]
            hlm[j]= 0.0
            j +=1
            m+=1
        else:
            if np.mod(i,2) == 1:
                glm[j] = ylm[k]
                j+=1
            else:
                hlm[m] = ylm[k]
                m+=1

        k+=1
# Set g00 = 0 (no monopoles) 
glm = np.insert(glm,0,0)
hlm = np.insert(hlm, 0, 0)*1j
Ylm = np.zeros(int(lmax*(lmax+3)/2+1),dtype=np.cdouble)
Ylm = glm + hlm

X,Y = get_grid(nphi=nphi, ntheta = ntheta)
Br = np.zeros_like(X)

cilm = pysh.shio.SHCindexToCilm(Ylm)
clm = pysh.SHCoeffs.from_array(cilm, normalization='ortho',csphase=-1)

# Z = sh.synth(Ylm)
Z = Br

w, h = plt.figaspect(fig_aspect)
fig, ax = plt.subplots(1, 1, figsize=(1.5*w, h),
                        subplot_kw={'projection': ccrs.Mollweide()})

# phi = np.linspace(0,2*np.pi,nphi)
# theta = np.linspace(0, np.pi, ntheta)
# X, Y = rad_to_deg(phi, theta)
# Z_lim = get_Z_lim(Z)
Z_lim = 20000
levels = np.linspace(-Z_lim, Z_lim, n_levels)
c = ax.contourf(X, Y, Z, levels, transform=ccrs.PlateCarree(), cmap='PuOr_r',
                      extend='both')

fig.show()
