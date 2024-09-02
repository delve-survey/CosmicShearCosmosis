import numpy as np, h5py
import healpy as hp
import sys, os
import multiprocessing as mp

FILES = {'wide' : '/project/chihway/data/decade/metacal_gold_combined_20240209.hdf',
         'n_of_z' : '/project/chihway/dhayaa/DECADE/SOMPZ/Runs/20240831/Summary/mean_nz_Final.npy',
         'cov_exec' : '/home/dhayaa/DECADE/CosmoCov/covs/cov',
        }


#STEP ONE: get the mask of the survey
with h5py.File(FILES['wide'], 'r') as f:
        
        ra   = f['RA'][:]
        dec  = f['DEC'][:]
        mask = f['baseline_mcal_mask_noshear'][:] > 0
        
        ra   = ra[mask]
        dec  = dec[mask]
        del mask

        
print("LOADED CATALOG")

Survey_Mask = np.bincount(hp.ang2pix(4096, ra, dec, lonlat = True), minlength = hp.nside2npix(4096)) > 0
Survey_Size = np.sum(Survey_Mask) * hp.nside2pixarea(4096, degrees = True) #Survey size in deg^2
Cls = hp.anafast(Survey_Mask, use_pixel_weights = True)
Cls = Cls/np.average(Survey_Mask) #Do 1/fsky normalization

np.savetxt('./Footprint_Cls.dat', np.stack([np.arange(Cls.size), Cls], axis = 1), fmt = ['%d', '%e'], delimiter = ' ')


print("FINISHED CLS")

#STEP TWO: rescale the number densities for cosmoscov friendly units
sigma_e   = np.array([0.233, 0.259, 0.248, 0.289, ])
sigma_ref = np.average(sigma_e)
n_eff     = np.array([1.239, 1.150, 1.169, 1.153, ])

new_n_eff = n_eff * (sigma_ref/sigma_e)**2
    
    
print("FINISHED RESCALING")

#STEP THREE: write out the n(z) file in right format
zbins  = np.arange(0.01, 5.05, 0.05)
zbinsc = zbins[:-1] + (zbins[1] - zbins[0])/2.
z_grid = zbinsc
n_of_z = np.load(FILES['n_of_z']).T

output = np.zeros([len(z_grid), 5])
output[:, 0] = z_grid - (z_grid[1] - z_grid[0])/2 #This is to go from z_bin_center to z_bin_min. Cosmocov wants latter

for i in range(4): output[:, i + 1] = n_of_z[:, i]
    
np.savetxt('./Redshifts.nz', output, fmt = ['%e'] * 5, delimiter = ' ')

print("FINISHED RESDHIFTS")

#STEP FOUR: write out the config
TEXT = """#
# Cosmological parameters
#
Omega_m : 0.3
Omega_v : 0.7
sigma_8 : 0.82355
n_spec : 0.97
w0 : -1.0
wa : 0.0
omb : 0.048
h0 : 0.69
coverH0 : 2997.92458
rho_crit : 7.4775e+21
f_NL : 0.0
pdelta_runmode : halofit
#
# Survey and galaxy parameters
#
# area in degrees
# n_gal,lens_n_gal in gals/arcmin^2
area : %(SURVEY_AREA)0.2f
sourcephotoz : multihisto
lensphotoz : 
source_tomobins : 4
lens_tomobins : 
sigma_e : %(SIGMA_REF_SQRT2)0.8f
shear_REDSHIFT_FILE : %(REDSHIFT_FILE)s
clustering_REDSHIFT_FILE : 
source_n_gal : %(Ngal1)0.8f,%(Ngal2)0.8f,%(Ngal3)0.8f,%(Ngal4)0.8f
lens_n_gal : 
lens_tomogbias :
lens_tomo_bmag : 
# IA parameters
IA : 0
A_ia : 0.0
eta_ia : 0.0
oneplusz0_ia : 1.6
#
# Covariance paramters
#
# tmin,tmax in arcminutes
tmin : 2.5
tmax : 250.0
ntheta : 20
ng : 1
cng : 1
c_footprint_file : %(FOOTPRINT_CLS_FILE)s
outdir : %(OUTPUT_DIR)s
filename : %(OUT_NAME)s
ss : true
ls : false
ll : false
"""

args = {'SURVEY_AREA' : Survey_Size,
        'SIGMA_REF_SQRT2'   : sigma_ref * np.sqrt(2),
        'REDSHIFT_FILE' : './Redshifts.nz',
        'Ngal1' : new_n_eff[0], 'Ngal2' : new_n_eff[1], 'Ngal3' : new_n_eff[2], 'Ngal4' : new_n_eff[3],
        'FOOTPRINT_CLS_FILE' : './Footprint_Cls.dat',
        'OUTPUT_DIR' : './outputs/',
        'OUT_NAME'   : 'Cov',
       }


#Print args for debugging state
print('-------INPUT PARAMS----------')
for p in args.keys():
    print('%s : %s'%(p.upper(), args[p]))
print('-----------------------------')
print('-----------------------------')

with open('./config.ini', 'w') as f:
    
    f.write(TEXT % args)
    
os.makedirs('./outputs', exist_ok = True)
print("FINISHED CONFIG")


#STEP FIVE: run the covariance runner in parallel across all cores.
def cov_block(i):
    os.system('%s %d ./config.ini' % (FILES['cov_exec'], i) )

    
with mp.Pool(os.cpu_count()) as p:
    p.map(cov_block, range(210 + 1)) #Hardcoded to 210 because this is the case for us forever :p
print("FINISHED COVARIANCE BLOCKS")

os.system('cat ./outputs/Cov_* > FinalCov.cov')

C = np.loadtxt('FinalCov.cov', comments = '#')

FinalCov = np.zeros([400, 400])
for i in range(C.shape[0]):
    FinalCov[int(C[i, 0]), int(C[i, 1])] = C[i, 8] + C[i, 9] #Add gaussian and non-gauss parts
    FinalCov[int(C[i, 1]), int(C[i, 0])] = FinalCov[int(C[i, 0]), int(C[i, 1])] * 1 #To ij and ji, multiply by 1 to avoid pointer/view :p
    
np.save('FinalCov.npy', FinalCov)

print("FINISHED COVARIANCE COMPILING")
