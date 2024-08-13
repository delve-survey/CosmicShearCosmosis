import numpy as np, h5py
from scipy import interpolate
import healpy as hp
import sys, os, shutil
import multiprocessing as mp

FILES = {'wide' : '/project/chihway/data/decade/metacal_gold_combined_20240209.hdf',
         'n_of_z' : '/project/chihway/dhayaa/DECADE/SOMPZ/Runs/20240408/Summary/mean_nz_combined_Final.npy',
         'cov_exec' : '/home/dhayaa/DECADE/CosmoCov/covs/cov',
        }


#########################################
# Bunch of init datasets
#########################################

#STEP ONE: get the mask of the survey
if not os.path.isfile('./Footprint_Cls.dat'):
    
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
    np.savetxt('./Survey_Size.txt', [Survey_Size])

TEXT = """#
# Cosmological parameters
#
Omega_m : %(Omega_m)0.4f
Omega_v : %(Omega_v)0.4f
sigma_8 : %(sigma_8)0.4f
n_spec : %(n_spec)0.4f
w0 : %(w0)0.4f
wa : 0.0
omb : %(omb)0.4f
h0 : %(h0)0.4f
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
IA : %(IA)d
A_ia : %(A_ia)0.4f
eta_ia : %(eta_ia)0.4f
oneplusz0_ia : 1.6
#
# Covariance paramters
#
# tmin,tmax in arcminutes
tmin : 2.5
tmax : 250.0
ntheta : 20
ng : 0
cng : 0
c_footprint_file : %(FOOTPRINT_CLS_FILE)s
outdir : %(OUTPUT_DIR)s
filename : Cov
ss : true
ls : false
ll : false
"""

#STEP FIVE: run the covariance runner in parallel across all cores.
def cov_block(i):

    if i == 0:
        os.system('%s %d ./config.ini' % (FILES['cov_exec'], i) )
    else:
        os.system('%s %d ./config.ini > /dev/null 2>&1' % (FILES['cov_exec'], i) )


class BaseRunner:

    def __init__(self, **kwargs):

        assert "Name" in kwargs.keys(), "Please provide input `Name` for the run during init"

        self.default_args = {'Omega_m' : 0.3,
                             'Omega_v' : 0.7,
                             'sigma_8' : 0.82355,
                             'n_spec'  : 0.97,
                             'w0'      : -1,
                             'omb'     : 0.048,
                             'h0'      : 0.69,
                             'IA'      : 0,
                             'A_ia'    : 0,
                             'eta_ia'  : 0,}

        for k in kwargs.keys():
            if k in self.default_args.keys():
                self.default_args[k] = kwargs[k]
                print(f"UPDATED {k} to {kwargs[k]}")

                if k == 'Omega_m':
                    self.default_args['Omega_v'] = 1 - self.default_args[k]
                    print(f"SETTING OMEGA_V TO {self.default_args['Omega_v']}")
                    assert 'Omega_v' not in kwargs.keys(), "Cannot provide both Omega_m and Omega_v"

        self.Name = kwargs['Name']
        self.kwargs = kwargs

    def get_neff_sigmae(self):

        #STEP TWO: rescale the number densities for cosmoscov friendly units
        sigma_e   = np.array([0.245, 0.269, 0.260, 0.294, ])
        sigma_ref = np.average(sigma_e)
        n_eff     = np.array([1.383, 1.342, 1.331, 1.284, ])

        new_n_eff = n_eff * (sigma_ref/sigma_e)**2

        return new_n_eff, sigma_ref

    def get_nz(self):

        zbins  = np.arange(0.01, 5.05, 0.05)
        zbinsc = zbins[:-1] + (zbins[1] - zbins[0])/2.
        z_grid = zbinsc
        n_of_z = np.load(FILES['n_of_z']).T

        output = np.zeros([len(z_grid), 5])
        output[:, 0] = z_grid - (z_grid[1] - z_grid[0])/2 #This is to go from z_bin_center to z_bin_min. Cosmocov wants latter

        for i in range(4): output[:, i + 1] = n_of_z[:, i]
            
        np.savetxt('./Redshifts.nz', output, fmt = ['%e'] * 5, delimiter = ' ')


    def process(self):

        new_n_eff, sigma_ref = self.get_neff_sigmae()
        self.get_nz()

        args = {'SURVEY_AREA' : np.loadtxt('./Survey_Size.txt'),
                'SIGMA_REF_SQRT2'   : sigma_ref * np.sqrt(2),
                'REDSHIFT_FILE' : './Redshifts.nz',
                'Ngal1' : new_n_eff[0], 'Ngal2' : new_n_eff[1], 'Ngal3' : new_n_eff[2], 'Ngal4' : new_n_eff[3],
                'FOOTPRINT_CLS_FILE' : './Footprint_Cls.dat',
                'OUTPUT_DIR' : './outputs/',
                'OUT_NAME'   : 'Cov',
        }

        args = args | self.default_args

        #Print args for debugging state
        print('-------INPUT PARAMS----------')
        for p in args.keys():
            print('%s : %s'%(p.upper(), args[p]))
        print('-----------------------------')
        print('-----------------------------')


        with open(f'./config.ini', 'w') as f:
            f.write(TEXT % args)

        os.makedirs('./outputs', exist_ok = True)
        print("FINISHED CONFIG")

        with mp.Pool(os.cpu_count()) as p:
            p.map(cov_block, range(210 + 1)) #Hardcoded to 210 because this is the case for us forever :p
        print("FINISHED COVARIANCE BLOCKS")
        os.system(f'cat ./outputs/Cov_* > FinalCov.cov')

        C = np.loadtxt('FinalCov.cov', comments = '#')

        FinalCov = np.zeros([400, 400])
        for i in range(C.shape[0]):
            FinalCov[int(C[i, 0]), int(C[i, 1])] = C[i, 8] + C[i, 9] #Add gaussian and non-gauss parts
            FinalCov[int(C[i, 1]), int(C[i, 0])] = FinalCov[int(C[i, 0]), int(C[i, 1])] * 1 #To ij and ji, multiply by 1 to avoid pointer/view :p
            
        np.save(f'{self.Name}Cov.npy', FinalCov)

        print("FINISHED COVARIANCE COMPILING")

        os.remove('FinalCov.cov')
        os.remove('config.ini')
        os.remove('Redshifts.nz')
        shutil.rmtree('./outputs')


class NZShiftRunner(BaseRunner):
    

    def get_nz(self):

        assert 'delta_z' in self.kwargs.keys(), "Please provide `delta_z` as input into kwargs"

        zbins  = np.arange(0.01, 5.05, 0.05)
        zbinsc = zbins[:-1] + (zbins[1] - zbins[0])/2.
        z_grid = zbinsc
        n_of_z = np.load(FILES['n_of_z']).T

        n_of_z = interpolate.PchipInterpolator(z_grid, n_of_z, axis = 0, extrapolate = False)(z_grid - self.kwargs['delta_z'])
        n_of_z = np.where(np.isnan(n_of_z), 0, n_of_z)
        print(f"SHIFTED n_of_z by {self.kwargs['delta_z']}")

        output = np.zeros([len(z_grid), 5])
        output[:, 0] = z_grid - (z_grid[1] - z_grid[0])/2 #This is to go from z_bin_center to z_bin_min. Cosmocov wants latter

        for i in range(4): output[:, i + 1] = n_of_z[:, i]
            
        np.savetxt('./Redshifts.nz', output, fmt = ['%e'] * 5, delimiter = ' ')


if __name__ == '__main__':

    
    #Base run
    BaseRunner(Name = 'Fiducial').process()
    
    #Cosmology difference runs
    BaseRunner(Name = 's8_0.9', sigma8 = 0.9).process()
    BaseRunner(Name = 's8_0.7', sigma8 = 0.7).process()

    BaseRunner(Name = 'Om_0.25', Omega_m = 0.25).process()
    BaseRunner(Name = 'Om_0.35', Omega_m = 0.35).process()

    BaseRunner(Name = 'h0_0.6', h0  = 0.6).process()
    BaseRunner(Name = 'h0_0.8', h0  = 0.8).process()

    #IA difference runs
    BaseRunner(Name = 'AIA_2', IA = 1, A_ia = 2).process()
    BaseRunner(Name = 'AIA_4', IA = 1, A_ia = 4).process()

    BaseRunner(Name = 'eta_IA_1',  IA = 1, A_ia = 1, eta_ia = 1).process()
    BaseRunner(Name = 'eta_IA_m1', IA = 1, A_ia = 1, eta_ia = -1).process()

    # #n(z) varying runs
    NZShiftRunner(Name = 'delta_z_0.01',  delta_z = 0.01).process()
    NZShiftRunner(Name = 'delta_z_0.05',  delta_z = 0.05).process()
    NZShiftRunner(Name = 'delta_z_m0.05', delta_z = -0.05).process()
    