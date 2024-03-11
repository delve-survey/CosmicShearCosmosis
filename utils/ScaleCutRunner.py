import numpy as np
from astropy.io import fits


class FiducialRunner:


    def __init__(self, FID, CONT):

        self.FID  = FID
        self.CONT = CONT
    
        self.Nind =  self.FID['XIP'].data['VALUE'].size
    
        #Read out all the quantities we need
        self.X = self._cosmosis_to_dvs(self.FID,   'value')
        self.Y = self._cosmosis_to_dvs(self.CONT,  'value')
        self.cov = self._cosmosis_to_cov(self.FID)
        self.ang = self._cosmosis_to_dvs(self.FID,  'ANGBIN')[:self.Nind]
        self.var = np.diagonal(self.cov)
        
    
    def _cosmosis_to_dvs(self, X, column = 'value'):
        '''
        Convert file from cosmosis format to the arrays we need
        '''

        xi = np.concatenate([X['XIP'].data[column.upper()], X['XIM'].data[column.upper()]])

        return xi
    
    
    def _cosmosis_to_cov(self, X,):

        return X['COVMAT'].data
    

    def compute_chi2(self, Mask):
        '''
        Standard chi2 calculation, of the masked datavectors
        '''

        x_masked   = self.X[Mask]
        y_masked   = self.Y[Mask]
        cov_masked = self.cov[Mask][:, Mask]

        res  = (x_masked - y_masked)
        chi2 = res @ np.linalg.inv(cov_masked) @ res 
        
        return chi2
    
    
    def cut_scale(self, Mask):

        #Load quantities
        X = self.X
        Y = self.Y
        var = self.var #Only need diagonal since we do it point by point

        #split into xi_p and xi_m
        X_p, X_m = X[:self.Nind], X[self.Nind:]
        Y_p, Y_m = Y[:self.Nind], Y[self.Nind:]
        M_p, M_m = Mask[:self.Nind], Mask[self.Nind:]

        var_p, var_m = var[:self.Nind], var[self.Nind:]

        ang = self.ang
        
        #measure residuals
        res_p = X_p[M_p] - Y_p[M_p]
        res_m = X_m[M_m] - Y_m[M_m]

        #Get the variance at the residuals
        var_p = var_p[M_p]
        var_m = var_m[M_m]
        
        #Find datapoints that are at smallest scale
        ang_mask_p = ang[M_p] == np.min(ang[M_p])
        ang_mask_m = ang[M_m] == np.min(ang[M_m])

        #I just compute the ang-bins, but don't use them here
        ang_p = ang[M_p][ang_mask_p]
        ang_m = ang[M_m][ang_mask_m]
        
        #Now compute chi2 at these points
        chi2_p = res_p[ang_mask_p]**2 / var_p[ang_mask_p]
        chi2_m = res_m[ang_mask_m]**2 / var_m[ang_mask_m]

        #If chi2 in xi_p is high then remove it. If the xip point is at higher ang scales
        #than the xim one, then remove the xim. This 2nd condition is a bit ad-hoc and I put it
        #in because there's one xi-m point at angbin==0 with very low chi2, so the algorithm
        #gets stuck on it and removes all xi_p points since those always have a higher chi2 than this onepoint
        #In reality, xim needs to have way more scales cut than xip.
        if (np.max(chi2_p) > np.max(chi2_m)) & (np.min(ang_p) <= np.min(ang_m)):

            bad_ind_in_sub  = np.argmax(chi2_p)
            bad_ind_in_mask = np.where(ang_mask_p)[0]
            bad_ind_in_dv   = np.where(M_p)[0]
            
            bad_ind = bad_ind_in_dv[bad_ind_in_mask[bad_ind_in_sub]]
            Mask[bad_ind] = False
            
        else:

            bad_ind_in_sub  = np.argmax(chi2_m)
            bad_ind_in_mask = np.where(ang_mask_m)[0]
            bad_ind_in_dv   = np.where(M_m)[0]
            
            bad_ind = self.Nind + bad_ind_in_dv[bad_ind_in_mask[bad_ind_in_sub]]
            Mask[bad_ind] = False
            
        return Mask


    def process(self, chi2_threshold):

        
        Mask = np.ones(2 * self.Nind, dtype = bool)
        chi2 = self.compute_chi2(Mask)

        while chi2 > chi2_threshold:

            Mask = self.cut_scale(Mask)
            chi2 = self.compute_chi2(Mask)


        return Mask
    
    
    def get_cuts(self, Mask):
        '''
        Takes in a mask of the full DV and then outputs scale cuts from xip and xim
        '''
        
        ang_val = self._cosmosis_to_dvs(self.FID,  'ANG')[:self.Nind]
        
        N = np.unique(self.ang).size
        
        M_p = Mask[:self.Nind]
        M_m = Mask[self.Nind:]
        
        xip_min = [np.min(ang_val[i*N:(i+1)*N][M_p[i*N:(i+1)*N]]) for i in range(self.Nind//N)]
        xim_min = [np.min(ang_val[i*N:(i+1)*N][M_m[i*N:(i+1)*N]]) for i in range(self.Nind//N)]
        
        return xip_min, xim_min
    
    
    def get_cuts_cosmosis(self, Mask):
        
        text = "[2pt_like]"
        
        ang_val = self._cosmosis_to_dvs(self.FID,  'ANG')[:self.Nind]
        
        N = np.unique(self.ang).size
        
        M_p = Mask[:self.Nind]
        M_m = Mask[self.Nind:]
        
        
        xip_min = [np.min(ang_val[i*N:(i+1)*N][M_p[i*N:(i+1)*N]]) for i in range(self.Nind//N)]
        xim_min = [np.min(ang_val[i*N:(i+1)*N][M_m[i*N:(i+1)*N]]) for i in range(self.Nind//N)]
        
        
        bins = combinations_with_replacement(range(4), 2)
        for i, b in enumerate(bins):
            text += "\nangle_range_xip_%d_%d = %0.3f 999.0" % (b[0], b[1], xip_min[i])
            
        
        text += "\n"
        bins = combinations_with_replacement(range(4), 2)
        for i, b in enumerate(bins):
            text += "\nangle_range_xim_%d_%d = %0.3f 999.0" % (b[0], b[1], xim_min[i])
            
            
        print(text)
        
        return None
    
    
    
class BinbyBinRunner(FiducialRunner):
    
    
    def cut_scale(self, Mask):

        X = self.X
        Y = self.Y
        var = self.var #Only need diagonal since we do it point by point

        #split into xi_p and xi_m

        X_p, X_m = X[:self.Nind], X[self.Nind:]
        Y_p, Y_m = Y[:self.Nind], Y[self.Nind:]

        M_p, M_m = Mask[:self.Nind], Mask[self.Nind:]

        var_p, var_m = var[:self.Nind], var[self.Nind:]

        ang = self.ang
        
        #measure worst point in xi_p and xi_m separately
        res_p = X_p[M_p] - Y_p[M_p]
        res_m = X_m[M_m] - Y_m[M_m]

        #Get their variance
        var_p = var_p[M_p]
        var_m = var_m[M_m]
 

        #Find datapoints that are at smallest scale, do it bin-by-bin
        
        ang_mask_p = np.concatenate([ang[i*20:(i+1)*20][M_p[i*20:(i+1)*20]] == np.min(ang[i*20:(i+1)*20][M_p[i*20:(i+1)*20]]) for i in range(10)])
        ang_mask_m = np.concatenate([ang[i*20:(i+1)*20][M_m[i*20:(i+1)*20]] == np.min(ang[i*20:(i+1)*20][M_m[i*20:(i+1)*20]]) for i in range(10)])


        ang_p = ang[M_p][ang_mask_p]
        ang_m = ang[M_m][ang_mask_m]
        
        #Now compute chi2 at these points
        chi2_p = res_p[ang_mask_p]**2 / var_p[ang_mask_p]
        chi2_m = res_m[ang_mask_m]**2 / var_m[ang_mask_m]

        if (np.max(chi2_p) > np.max(chi2_m)):

            bad_ind_in_sub  = np.argmax(chi2_p)
            bad_ind_in_mask = np.where(ang_mask_p)[0]
            bad_ind_in_dv   = np.where(M_p)[0]
            
            bad_ind = bad_ind_in_dv[bad_ind_in_mask[bad_ind_in_sub]]
            Mask[bad_ind] = False
            
        else:

            bad_ind_in_sub  = np.argmax(chi2_m)
            bad_ind_in_mask = np.where(ang_mask_m)[0]
            bad_ind_in_dv   = np.where(M_m)[0]
            
            bad_ind = self.Nind + bad_ind_in_dv[bad_ind_in_mask[bad_ind_in_sub]]
            Mask[bad_ind] = False
            
        return Mask
    