import numpy as np, h5py
import healpy as hp
import sys, os
import multiprocessing as mp

FILES = {'wide' : '/project/chihway/data/decade/metacal_gold_combined_20240209.hdf',}


def printfunction(neff, sigma_e, R11, R22):

    print(f"{neff:.3f}  {sigma_e:.3f}  {sigma_e**2/neff:.3f}  {R11:.3f}  {R22:.3f}")

Rb = np.array([ (0.8726000168016929, 0.8727850133159067),
                (0.7977613141932485, 0.798284252899401),
                (0.755664376239435, 0.7577162611289076),
                (0.6253332865597496, 0.6262105017778515)])

#STEP ONE: get the mask of the survey
with h5py.File(FILES['wide'], 'r') as f:
        
    ra   = f['RA'][:]
    dec  = f['DEC'][:]
    mask = f['baseline_mcal_mask_noshear'][:] > 0
    
    ra   = ra[mask]
    dec  = dec[mask]

    g1, g2 = f['mcal_g_noshear'][:][mask].T
    w      = f['mcal_g_w'][:]

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

####################################################################################################################
# NOW WE DO DES
####################################################################################################################

Rb = np.array([ (0.7636 + 0.0046, 0.7636 + 0.0046),
                (0.7182 + 0.0083, 0.7182 + 0.0083),
                (0.6887 + 0.0126, 0.6887 + 0.0126),
                (0.6154 + 0.0145, 0.6154 + 0.0145)])

# with h5py.File('/project/chihway/dhayaa/DES_Catalogs/DESY3_indexcat.h5') as f:
#     mcal_selection = f['index/select'][:]

# with h5py.File('/project/chihway/dhayaa/DES_Catalogs/DESY3_metacal_v03-004.h5') as f:

#     ra   = f['catalog/unsheared/ra'][:][mcal_selection]
#     dec  = f['catalog/unsheared/dec'][:][mcal_selection]


#     Survey_Mask = np.bincount(hp.ang2pix(4096, ra, dec, lonlat = True), minlength = hp.nside2npix(4096)) > 0
#     Survey_Size = np.sum(Survey_Mask) * hp.nside2pixarea(4096, degrees = True) #Survey size in deg^2
#     print("SURVEY_SIZE :", Survey_Size)

Survey_Size = 4143
print("FORCING SURVEY SIZE TO BE :", Survey_Size)


with h5py.File('/project/chihway/dhayaa/DES_Catalogs/DESY3_metacal_v03-004.h5') as f:
    w = f['catalog/unsheared/weight'][:]

for NAME, WEIGHTS in zip(['Normal', 'Inverse', 'Nowgt'], [w, 1/w, np.ones_like(w)]):
    print(f'============= {NAME} (DES) ================')
    for i in range(4):

        with h5py.File('/project/chihway/dhayaa/DES_Catalogs/DESY3_indexcat.h5') as s:
            with h5py.File('/project/chihway/dhayaa/DES_Catalogs/DESY3_metacal_v03-004.h5') as f:

                dgamma = 2*0.01
                R11  = (np.average(f['catalog/sheared_1p/e_1'][:][s[f'index/select_bin{i + 1}'][:]],    weights = WEIGHTS[s[f'index/select_bin{i + 1}'][:]]) - 
                        np.average(f['catalog/sheared_1m/e_1'][:][s[f'index/select_bin{i + 1}'][:]],    weights = WEIGHTS[s[f'index/select_bin{i + 1}'][:]])) / dgamma
                R11s = (np.average(f['catalog/unsheared/e_1' ][:][s[f'index/select_1p_bin{i + 1}'][:]], weights = WEIGHTS[s[f'index/select_1p_bin{i + 1}'][:]]) - 
                        np.average(f['catalog/unsheared/e_1' ][:][s[f'index/select_1m_bin{i + 1}'][:]], weights = WEIGHTS[s[f'index/select_1m_bin{i + 1}'][:]])) / dgamma
                R22  = (np.average(f['catalog/sheared_2p/e_2'][:][s[f'index/select_bin{i + 1}'][:]],    weights = WEIGHTS[s[f'index/select_bin{i + 1}'][:]]) - 
                        np.average(f['catalog/sheared_2m/e_2'][:][s[f'index/select_bin{i + 1}'][:]],    weights = WEIGHTS[s[f'index/select_bin{i + 1}'][:]])) / dgamma
                R22s = (np.average(f['catalog/unsheared/e_2' ][:][s[f'index/select_2p_bin{i + 1}'][:]], weights = WEIGHTS[s[f'index/select_2p_bin{i + 1}'][:]]) - 
                        np.average(f['catalog/unsheared/e_2' ][:][s[f'index/select_2m_bin{i + 1}'][:]], weights = WEIGHTS[s[f'index/select_2m_bin{i + 1}'][:]])) / dgamma
                
                R11  = R11 + R11s
                R22  = R22 + R22s

                wi   = WEIGHTS[s[f'index/select_bin{i + 1}'][:]]
                g1   = f['catalog/unsheared/e_1' ][:][s[f'index/select_bin{i + 1}'][:]]
                g2   = f['catalog/unsheared/e_2' ][:][s[f'index/select_bin{i + 1}'][:]]
                neff = np.sum(wi)**2/np.sum(wi**2) / (Survey_Size * 60 * 60)
                sigma_e = np.sqrt(np.average( (g1/R11)**2 +  (g2/R22)**2, weights = wi**2)/2)
                printfunction(neff, sigma_e, R11, R22)

    print('\n')