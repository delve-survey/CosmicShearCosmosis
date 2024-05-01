import numpy as np
import sys, os
from astropy.io import fits
from astropy.table import QTable


OUTFILE = "baseline_April23rd2024.fits"

#Clean up first
os.system('rm %s' % OUTFILE)
os.system('rm %s' % OUTFILE.replace('baseline', 'baryon') )


FILES = {'n_of_z' : '/project/chihway/dhayaa/DECADE/SOMPZ/Runs/20240408/n_of_z.npy',
         'z_grid' : '/project/chihway/dhayaa/DECADE/SOMPZ/Runs/20240408/z_grid.npy',
         'Cov' : '/home/dhayaa/DECADE/CosmicShearCosmosis/covariance/analytic/April23rd2024_analytic/FinalCov.npy',
        }

#STEP 1: Open the placeholder file and put into it the Covariance and n(z) that we want
f = fits.open('placeholder.fits')

covmat, header_covmat = fits.getdata('placeholder.fits', 'COVMAT',    header=True)
source, header_source = fits.getdata('placeholder.fits', 'nz_source', header=True)
xip,    header_xip    = fits.getdata('placeholder.fits', 'xip', header=True)
xim,    header_xim    = fits.getdata('placeholder.fits', 'xim', header=True)


z_grid = np.load(FILES['z_grid'])
n_of_z = np.load(FILES['n_of_z'])


header_source['NAXIS2'] = len(z_grid)
source = QTable()
source['Z_LOW']  = z_grid - (z_grid[1] - z_grid[0])
source['Z_MID']  = z_grid
source['Z_HIGH'] = z_grid + (z_grid[1] - z_grid[0])
source['BIN1'] = n_of_z[0]
source['BIN2'] = n_of_z[1]
source['BIN3'] = n_of_z[2]
source['BIN4'] = n_of_z[3]

covmat_NEW = np.load(FILES['Cov'])

base_header = fits.open('placeholder.fits')
fits.append(OUTFILE, base_header[0].data, base_header[0].header)
fits.append(OUTFILE, covmat_NEW, header_covmat)
fits.append(OUTFILE, xip, header_xip)
fits.append(OUTFILE, xim, header_xim)
fits.append(OUTFILE, source, header_source)
base_header.close()

#Setup the ini files
INI = """[runtime]
sampler = test
verbosity = standard

; Parameters used several times in this file are
; put in the DEFAULT section and are referenced as e.g. %(2PT_FILE)s
[DEFAULT]
2PT_FILE = baseline_April23rd2024.fits
BASE_DIR = /home/dhayaa/DECADE/cosmosis-standard-library 
MODE = save_2pt


[pipeline]
modules =  consistency  bbn_consistency
           camb sigma8_rescale fast_pt
           fits_nz IA pk_to_cl
           add_intrinsic
           2pt_shear
           shear_m_bias
           %(MODE)s

timing=F
debug=F
priors = priors.ini
values = values.ini
extra_output = cosmological_parameters/sigma_8 data_vector/2pt_chi2

[baryon]
file = %(BASE_DIR)s/structure/baryon_power_scaling/baryonic_interface.py
ratio_table = %(BASE_DIR)s/structure/baryon_power_scaling/data/logPkRatio_owls_AGN.dat


[test]
save_dir=outputs/shear_nontomo
fatal_errors=T


[consistency]
file = %(BASE_DIR)s/utility/consistency/consistency_interface.py

[camb]
file = %(BASE_DIR)s/boltzmann/camb/camb_interface.py
mode = all
lmax = 2500          ;max ell to use for cmb calculation
feedback=3         ;amount of output to print
AccuracyBoost=1.1 ;CAMB accuracy boost parameter
do_tensors = T
do_lensing = T
NonLinear = pk
halofit_version = takahashi
zmin_background = 0.
zmax_background = 4.
nz_background = 401
kmin=1e-4
kmax = 50.0
kmax_extrapolate = 500.0
nk=700


[bbn_consistency]
file = %(BASE_DIR)s/utility/bbn_consistency/bbn_consistency.py

[sigma8_rescale]
file = %(BASE_DIR)s/utility/sample_sigma8/sigma8_rescale.py

[fits_nz]
file = %(BASE_DIR)s/number_density/load_nz_fits/load_nz_fits.py
nz_file = %(2PT_FILE)s
data_sets = source
prefix_section = T
prefix_extension = T

[source_photoz_bias]
file = %(BASE_DIR)s/number_density/photoz_bias/photoz_bias.py
mode = additive
sample = nz_source
bias_section = wl_photoz_errors
interpolation = linear

[fast_pt]
file = %(BASE_DIR)s/structure/fast_pt/fast_pt_interface.py
do_ia = T
k_res_fac = 0.5
verbose = F

[IA]
file = %(BASE_DIR)s/intrinsic_alignments/tatt/tatt_interface.py
sub_lowk = F
do_galaxy_intrinsic = F
ia_model = nla

[pk_to_cl]
file = %(BASE_DIR)s/structure/projection/project_2d.py
ell_min_logspaced = 0.1
ell_max_logspaced = 5.0e5
n_ell_logspaced = 100 
shear-shear = source-source
shear-intrinsic = source-source
intrinsic-intrinsic = source-source
intrinsicb-intrinsicb=source-source
verbose = F
get_kernel_peaks = F
sig_over_dchi = 20. 
shear_kernel_dchi = 10. 

[add_intrinsic]
file=%(BASE_DIR)s/shear/add_intrinsic/add_intrinsic.py
shear-shear=T
position-shear=F
perbin=F

[add_eb]
file = %(BASE_DIR)s/intrinsic_alignments/add_b_mode/add_b_mode_cl.py

[2pt_shear]
file = %(BASE_DIR)s/shear/cl_to_xi_fullsky/cl_to_xi_interface.py
ell_max = 40000
xi_type = EB
theta_file=%(2PT_FILE)s
bin_avg = T
; these get
input_section_name = shear_cl  shear_cl_bb
output_section_name = shear_xi_plus  shear_xi_minus


[shear_m_bias]
file = %(BASE_DIR)s/shear/shear_bias/shear_m_bias.py
m_per_bin = True
; Despite the parameter name, this can operate on xi as well as C_ell.
cl_section = shear_xi_plus shear_xi_minus
;cross_section = galaxy_shear_xi
verbose = F

[2pt_like]
file = %(BASE_DIR)s/likelihood/2pt/2pt_point_mass/2pt_point_mass.py
do_pm_marg = False
do_pm_sigcritinv = False
sigma_a = 10000.0
no_det_fac = False
include_norm = False ;True
data_file = %(2PT_FILE)s
data_sets = xip xim
make_covariance=F
covmat_name=COVMAT


[save_2pt]
file = %(BASE_DIR)s/likelihood/2pt/save_2pt.py
theta_min = 2.5
theta_max = 250.0
n_theta = 20
real_space = T
make_covariance = F
shear_nz_name = nz_source
position_nz_name = nz_lens
filename = %(2PT_FILE)s
overwrite = T
auto_only = galaxy_xi
spectrum_sections = shear_xi_plus shear_xi_minus 
output_extensions = xip xim 
two_thirds_midpoint = T
copy_covariance=%(2PT_FILE)s
"""

VALUES = """[cosmological_parameters]
omega_m      =  0.3
h0           =  0.69
omega_b      =  0.048
n_s          =  0.97
w            =  -1.0
A_s     = 2.1e-9 ;placeholder, doesn't matter as sigma8 will be rescaled to match the cosmolike input
sigma8_input = 0.82355

; New parametrization and names:
mnu = 0.0
num_massive_neutrinos = 0
nnu = 3.046

omega_k      =  0.0
tau          =  0.0697186

[shear_calibration_parameters]
m1 = -0.0063
m2 = -0.0198
m3 = -0.0241
m4 = -0.0369

[intrinsic_alignment_parameters]
z_piv   =  0.62
A1      = 0.19
A2      = 0.0
alpha1  = -2.6
alpha2  = 0.0
bias_ta = 0.0
"""

PRIORS = """
"""


with open('./shear.ini', 'w') as f:  f.write(INI)
with open('./values.ini', 'w') as f: f.write(VALUES)
with open('./priors.ini', 'w') as f: f.write(PRIORS)
    
with open('./shear_test.ini', 'w') as f:  f.write( INI.replace("save_2pt", "2pt_like") )

os.system('cosmosis shear.ini')
os.system('cosmosis shear_test.ini')



with open('./shear_baryon.ini', 'w') as f:   f.write( (INI
                                                       .replace("camb sigma8_rescale", "camb sigma8_rescale baryon")
                                                       .replace("filename = %(2PT_FILE)s", 
                                                                "filename = %s" % OUTFILE.replace('baseline', 'baryon'))
                                                      )
                                                    )
os.system('cosmosis shear_baryon.ini')