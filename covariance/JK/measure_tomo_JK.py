import h5py as h
import numpy as np
import treecorr
import os
import time

bslop = 0.001
PROCESS = int(os.environ['SLURM_PROCID']) #sbatching a job on 10 nodes will spawn 10 tasks (if you do it right) 
THREADS = int(os.environ['OMP_NUM_THREADS']) #must be exported in the batch job script, 32 for cori/haswell and 68 for cori/knl
outpath = '/home/secco/CosmicShearCosmosis/covariance/JK/output_JK/'
catlocation = '/home/secco/CosmicShearCosmosis/covariance/JK/'

#Deciding which process does which bin combination i,j:
#Will only work correctly if you start 10 processes

if PROCESS==0: 
    i,j = (1,1)
if PROCESS==1: 
    i,j = (1,2)
if PROCESS==2: 
    i,j = (1,3)
if PROCESS==3:
    i,j = (1,4)
if PROCESS==4:
    i,j = (2,2)
if PROCESS==5:
    i,j = (2,3)
if PROCESS==6:
    i,j = (2,4)
if PROCESS==7:
    i,j = (3,3)
if PROCESS==8:
    i,j = (3,4)
if PROCESS==9:
    i,j = (4,4)

time1=time.time()    
catalog_i = np.load(catlocation+'catalog_tag0613_bin%d.npy'%(i-1))
ra_i, dec_i, g1_i, g2_i, w_i = catalog_i['ra'], catalog_i['dec'], catalog_i['g1'], catalog_i['g2'], catalog_i['w']  

catalog_j = np.load(catlocation+'catalog_tag0613_bin%d.npy'%(j-1))
ra_j, dec_j, g1_j, g2_j, w_j = catalog_j['ra'], catalog_j['dec'], catalog_j['g1'], catalog_j['g2'], catalog_j['w']  
time2=time.time()

print('\n\ntreecorr cat length: zbin%d=%d, zbin%d=%d'%(i,len(g1_i),j,len(g1_j)))
print('Loading the catalog data took %1.2f seconds'%(time2-time1))

centers_file = catlocation+'tag0613_centers.txt'
cat1 = treecorr.Catalog(g1=g1_i, g2=g2_i, ra=ra_i, dec=dec_i, w=w_i, ra_units='deg',dec_units='deg',patch_centers=centers_file)
cat2 = treecorr.Catalog(g1=g1_j, g2=g2_j, ra=ra_j, dec=dec_j, w=w_j, ra_units='deg',dec_units='deg',patch_centers=centers_file)

print('ID for this job: ',PROCESS,'\nDoing combination i=',i, 'j=',j,'\nUsing',THREADS,'cores')
GG = treecorr.GGCorrelation(nbins=20,min_sep=2.5,max_sep=250.0,sep_units='arcmin',verbose=3,bin_slop=bslop)
GG.process(cat1,cat2,num_threads=THREADS)
cov_jk = GG.estimate_cov('jackknife')
np.save(outpath+'xipm_100patches_bslop'+str(bslop)+'_'+str(i)+str(j)+'_JKcov.npy',cov_jk)
cov_sample = GG.estimate_cov('sample')
np.save(outpath+'xipm_100patches_bslop'+str(bslop)+'_'+str(i)+str(j)+'_SAMPLEcov.npy',cov_sample)

GG.write(outpath+'shear2pt_tomo_bslop'+str(bslop)+'_'+str(i)+str(j)+'_JK.txt')
time3=time.time()
print('Making the measurement took %1.2f seconds'%(time3-time2))
