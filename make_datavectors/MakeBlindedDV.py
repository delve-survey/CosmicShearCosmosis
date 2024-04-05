import numpy as np
from astropy.io import fits
import sys

fid_datavec = sys.argv[1]
sim_dv = fits.open(fid_datavec)

data_dir = sys.argv[2]
output=sys.argv[3]

n = 0
for i in range(4):
    for j in range(4):
        if j>=i:
            data = np.loadtxt(data_dir+'gg_'+str(i+1)+'_'+str(j+1)+'/mean_gg')

            for i in range(len(data)):
                sim_dv['xip'].data['VALUE'][n+i] = data[i][1]
                sim_dv['xim'].data['VALUE'][n+i] = data[i][2]
                n+=1

sim_dv.writeto(output, overwrite=True)   

