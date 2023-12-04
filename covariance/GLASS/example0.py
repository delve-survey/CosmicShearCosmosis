import numpy as np
import camb
from cosmology import Cosmology

import glass.shells
import glass.ext.camb


# cosmology for the simulation
h = 0.7
Oc = 0.25
Ob = 0.05

# basic parameters of the simulation
lmax = 1000

# set up CAMB parameters for matter angular power spectrum
pars = camb.set_params(H0=100*h, omch2=Oc*h**2, ombh2=Ob*h**2,
                       NonLinear=camb.model.NonLinear_both)

# get the cosmology from CAMB
cosmo = Cosmology.from_camb(pars)

# shells of 200 Mpc in comoving distance spacing
zb = glass.shells.distance_grid(cosmo, 0., 1., dx=200.)

# uniform matter weight function
# CAMB requires linear ramp for low redshifts
ws = glass.shells.tophat_windows(zb, weight=glass.ext.camb.camb_tophat_weight)

# compute angular matter power spectra with CAMB
cls = glass.ext.camb.matter_cls(pars, lmax, ws)

np.save('cls.npy', cls)
