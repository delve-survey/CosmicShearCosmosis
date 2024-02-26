import numpy as np
import fitsio


class FiducialRunner:


    def __init__(self, FID, CONT):

        self.FID  = FID
        self.CONT = CONT
    
        self.Nind =  self.FID['XIP']['VALUE'].size
    
    def _cosmosis_to_dvs(self, X, column = 'value'):
        '''
        Convert file from cosmosis format to the arrays we need
        '''

        xi = np.concatenate([X['XIP'][column.upper()], X['XIM'][column.upper()]])

        return xi
    
    
    def _cosmosis_to_cov(self, X,):

        return X['COV']
    

    def compute_chi2(self, Mask):

        x_masked   = self._cosmosis_to_dvs(self.FID,   'value')[Mask]
        y_masked   = self._cosmosis_to_dvs(self.CONT,  'value')[Mask]
        cov_masked = self._cosmosis_to_cov(self.FID)[Mask][:, Mask]

        res  = (x_masked - y_masked)
        chi2 = res @ np.linalg.inv(cov_masked) @ res 
        
        return chi2
    
    
    def cut_scale(self, Mask):

        X = self._cosmosis_to_dvs(self.FID,  'value')
        Y = self._cosmosis_to_dvs(self.CONT, 'value')

        err = np.diagonal(self._cosmosis_to_cov(self.FID)) #Only need diagonal since we do it point by point

        #split into xi_p and xi_m

        X_p, X_m = X[:self.Nind], X[self.Nind:]
        Y_p, Y_m = Y[:self.Nind], Y[self.Nind:]

        M_p, M_m = Mask[:self.Nind], Mask[self.Nind:]

        err_p, err_m = err[:self.Nind], err[self.Nind:]

        ang = self._cosmosis_to_dvs(self.FID,  'ANG')
        
        #measure worst point in xi_p and xi_m separately
        res_p = X_p[M_p] - Y_p[M_p]
        res_m = X_m[M_m] - Y_m[M_m]

        #Get their errors
        err_p = err_p[M_p]
        err_m = err_m[M_m]


        #Find datapoints that are at smallest scale
        ang_mask_p = ang[M_p] == np.min(ang[M_p])
        ang_mask_m = ang[M_m] == np.min(ang[M_m])


        #Now compute chi2 at these points
        chi2_p = res_p[ang_mask_p]**2 / err_p[ang_mask_p]
        chi2_m = res_m[ang_mask_m]**2 / err_p[ang_mask_m]

        if np.max(chi2_p) > np.max(chi2_m):

            bad_ind_in_sub  = np.argmax(chi2_p)
            bad_ind_in_mask = np.where(ang_mask_p)[0]
            Mask[bad_ind_in_mask[bad_ind_in_sub]] = False

        else:

            bad_ind_in_sub  = np.argmax(chi2_m)
            bad_ind_in_mask = np.where(ang_mask_m)[0]
            Mask[self.Nind + bad_ind_in_mask[bad_ind_in_sub]] = False


        return Mask


    def process(self, chi2_threshold):

        
        Mask = np.ones(2 * self.Nind, dtype = bool)
        chi2 = self.compute_chi2(Mask)

        while chi2 > chi2_threshold:

            Mask = self.cut_scale(Mask)
            chi2 = self.compute_chi2(Mask)


        return Mask