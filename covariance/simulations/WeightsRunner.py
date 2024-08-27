import numpy as np, h5py
import sys, os, gc
import joblib
from scipy import interpolate

prefix = np.random.default_rng(seed = 42).integers(2**30)
TMPDIR = os.environ['TMPDIR']
base   = TMPDIR + f'/{prefix}'
dgamma = 2*0.01


def create_reader(file_suffix):
    def reader(mask = None):
        file_path = f'{base}_{file_suffix}.npy'
        data = np.load(file_path, mmap_mode='r')
        
        if mask is None:
            return data
        else:
            return data[mask]
        
    return reader

func = {}

#STEP ONE: build a subset of the quantities to use in defining the weights
with h5py.File('/project/chihway/data/decade/metacal_gold_combined_20240209.hdf', 'r') as f:
        
    
    #Get a mask of all objects that would be used
    mask = False
    for m in ['noshear', '1p', '1m', '2p', '2m']:
        mask = mask | (f[f'baseline_mcal_mask_{m}'][:] > 0)


    for m in ['noshear', '1p', '1m', '2p', '2m']:

        np.save(f'{base}_mask_{m}.npy', f[f'baseline_mcal_mask_{m}'][:][mask])

        for q in ['mcal_g', 'mcal_s2n', 'mcal_T_ratio']:    
            np.save(f'{base}_{q}_{m}.npy', f[f'{q}_{m}'][:][mask])

            func[f'{q}_{m}'] = create_reader(f'{q}_{m}')


    # #Ran this just to determine the limits of the S2N and T_ratio of Grid
    # m = np.load(f'{base}_mask_noshear.npy') > 0
    # print("S2N [95%] :", np.percentile(np.load(f'{base}_mcal_s2n_noshear.npy')[m],     [2.5, 97.5]))
    # print("Tr  [95%] :", np.percentile(np.load(f'{base}_mcal_T_ratio_noshear.npy')[m], [2.5, 97.5]))

def one_step(SNR_range, Tratio_range, bin = None):
    
    mask_h = {}

    for m in ['noshear', '1p', '1m', '2p', '2m']:
        snr         = func[f'mcal_s2n_{m}']()
        Tr          = func[f'mcal_T_ratio_{m}']()
        mask_h[m]   = (snr > SNR_range[0]) & (snr < SNR_range[1]) & (Tr > Tratio_range[0]) & (Tr < Tratio_range[1])
    
    del Tr, snr; gc.collect()

    if bin is None:
        mask = {m : (np.load(f'{base}_mask_{m}.npy', mmap_mode = 'r') > 0) & mask_h[m] for m in mask_h.keys()}
    else:
        mask = {m : (np.load(f'{base}_mask_{m}.npy', mmap_mode = 'r') == bin) & mask_h[m] for m in mask_h.keys()}

    del mask_h; gc.collect()

    #Get sigma_e in this cell
    g1, g2   = func['mcal_g_noshear'](mask['noshear']).T
    sigma_e  = np.sqrt(np.average(g1**2 + g2**2)/2)
    Ncounts  = np.sum(mask['noshear'])
    del g1, g2; gc.collect()
    
    #Now get the responses
    R11      = (np.average(func['mcal_g_1p'](mask['noshear'])[:, 0]) - np.average(func['mcal_g_1m'](mask['noshear'])[:, 0])) / dgamma    
    R11s     = (np.average(func['mcal_g_noshear'](mask['1p'])[:, 0]) - np.average(func['mcal_g_noshear'](mask['1m'])[:, 0])) / dgamma    
    R22      = (np.average(func['mcal_g_2p'](mask['noshear'])[:, 1]) - np.average(func['mcal_g_2m'](mask['noshear'])[:, 1])) / dgamma    
    R22s     = (np.average(func['mcal_g_noshear'](mask['2p'])[:, 1]) - np.average(func['mcal_g_noshear'](mask['2m'])[:, 1])) / dgamma   

    del mask; gc.collect()

    R11_tot  = R11 + R11s
    R22_tot  = R22 + R22s    
    Rg       = (R11 + R22)/2
    weight   = (Rg/sigma_e)**2

    SNR = np.average(SNR_range)
    Tr  = np.average(Tratio_range)

    dtype = ['weight', 'sigma_e', 'Rgamma', 'Ncounts', 'R11_tot', 'R22_tot', 'R11', 'R11s', 'R22', 'R22s', 'SNR', 'T_ratio']
    dtype = [(k, float) for k in dtype]
    out   = np.array([tuple([weight, sigma_e, Rg, Ncounts, R11_tot, R22_tot, R11, R11s, R22, R22s, SNR, Tr])], dtype = dtype)
    return out


N    = 20
SNR  = np.geomspace(10,  330, N + 1)
Tr   = np.geomspace(0.5, 4.5, N + 1)
size = N * N

jobs    = [joblib.delayed(one_step)([SNR[i % N], SNR[i % N +1]], [Tr[i // N], Tr[i // N + 1]], bin = None) for i in range(size)]
results = joblib.Parallel(n_jobs = -1, verbose = 10)(jobs)
results = np.concatenate(results)

np.save('./grid_quantities_20240827.npy', results)

results = np.load('./grid_quantities_20240827.npy')
weights = results['weight']

import healpy as hp
FILES = {'wide' : '/project/chihway/data/decade/metacal_gold_combined_20240209.hdf',}


def printfunction(neff, sigma_e, R11, R22):
    print(f"{neff:.3f}  {sigma_e:.3f}  {sigma_e**2/neff:.3f}  {R11:.3f}  {R22:.3f}")


#STEP ONE: get the mask of the survey
with h5py.File(FILES['wide'], 'r') as f:
        
    ra   = f['RA'][:]
    dec  = f['DEC'][:]
    mask = f['baseline_mcal_mask_noshear'][:] > 0
    
    ra   = ra[mask]
    dec  = dec[mask]

    g1, g2 = f['mcal_g_noshear'][:][mask].T
    w      = f['mcal_g_w'][:]

    snr    = f['mcal_s2n_noshear'][:]
    Tr     = f['mcal_T_ratio_noshear'][:]

    snr = np.nan_to_num(snr)
    Tr  = np.nan_to_num(Tr)
    snr_line = results['SNR']
    Tr_line  = results['T_ratio']
    w_line   = results['weight']
    I   = interpolate.NearestNDInterpolator(np.vstack([snr_line, Tr_line]).T, w_line)
    w   = I(snr, Tr)

    Survey_Mask = np.bincount(hp.ang2pix(4096, ra, dec, lonlat = True), minlength = hp.nside2npix(4096)) > 0
    Survey_Size = np.sum(Survey_Mask) * hp.nside2pixarea(4096, degrees = True) #Survey size in deg^2
    print("SURVEY_SIZE :", Survey_Size)

    del ra, dec

    mask = f['baseline_mcal_mask_noshear'][:][mask]

    for NAME, WEIGHTS in zip(['Normal', 'Inverse', 'Nowgt'], [w, 1/w, np.ones_like(w)]):
        print(f'============= {NAME} ================')
        for i in range(4):
            bin  = mask == (i + 1)

            dgamma = 2*0.01
            R11  = (np.average(f['mcal_g_1p'][:, 0][f['baseline_mcal_mask_noshear'][:] == (i + 1)], weights = WEIGHTS[f['baseline_mcal_mask_noshear'][:] == (i + 1)]) - 
                    np.average(f['mcal_g_1m'][:, 0][f['baseline_mcal_mask_noshear'][:] == (i + 1)], weights = WEIGHTS[f['baseline_mcal_mask_noshear'][:] == (i + 1)])) / dgamma
            R11s = (np.average(f['mcal_g_noshear'][:, 0][f['baseline_mcal_mask_1p'][:] == (i + 1)], weights = WEIGHTS[f['baseline_mcal_mask_1p'][:] == (i + 1)]) - 
                    np.average(f['mcal_g_noshear'][:, 0][f['baseline_mcal_mask_1m'][:] == (i + 1)], weights = WEIGHTS[f['baseline_mcal_mask_1m'][:] == (i + 1)])) / dgamma
            R22  = (np.average(f['mcal_g_2p'][:, 1][f['baseline_mcal_mask_noshear'][:] == (i + 1)], weights = WEIGHTS[f['baseline_mcal_mask_noshear'][:] == (i + 1)]) - 
                    np.average(f['mcal_g_2m'][:, 1][f['baseline_mcal_mask_noshear'][:] == (i + 1)], weights = WEIGHTS[f['baseline_mcal_mask_noshear'][:] == (i + 1)])) / dgamma
            R22s = (np.average(f['mcal_g_noshear'][:, 1][f['baseline_mcal_mask_2p'][:] == (i + 1)], weights = WEIGHTS[f['baseline_mcal_mask_2p'][:] == (i + 1)]) - 
                    np.average(f['mcal_g_noshear'][:, 1][f['baseline_mcal_mask_2m'][:] == (i + 1)], weights = WEIGHTS[f['baseline_mcal_mask_2m'][:] == (i + 1)])) / dgamma
            
            R11  = R11 + R11s
            R22  = R22 + R22s

            wi = WEIGHTS[f['baseline_mcal_mask_noshear'][:] == (i + 1)]
            g1 = f['mcal_g_noshear'][:, 0][f['baseline_mcal_mask_noshear'][:] == (i + 1)]
            g2 = f['mcal_g_noshear'][:, 1][f['baseline_mcal_mask_noshear'][:] == (i + 1)]

            neff = np.sum(wi)**2/np.sum(wi**2) / (Survey_Size * 60 * 60)
            sigma_e = np.sqrt(np.average( (g1/R11)**2 +  (g2/R22)**2, weights = wi**2)/2)
            printfunction(neff, sigma_e, R11, R22)

        print('\n')