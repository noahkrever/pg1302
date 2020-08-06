# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 20:46:33 2020

@author: noahk
"""

# import sys

def separateCurves(filename): #splits up curves into separate files for parallel cluster runs
    import h5py
    f=h5py.File(filename,"r")
    for i in f:
        print(i)
        f_temp1 = h5py.File('L+C+A_mcmc_rand_sig_tau_no_outlier_' + str(i) + '.h5', "a")
        # f_temp2 = h5py.File('final_curve_pg1302_LCA.h5', "a")
        
        f_temp1[str(i)+'/MJD']=f['/'+str(i)+'/MJD'].value
        # f[str(i)+'/mag']=x
        f_temp1[str(i)+'/mag']=f['/'+str(i)+'/mag'].value
        # f_temp[str(i)+'/P']= f['/'+str(i)+'/P'].value
        # f_temp[str(i)+'/A']=f[str(i)+'/A'].value
        # # f[str(i)+'/phi']=phi
        # f_temp[str(i)+'/t0']=f[str(i)+'/t0'].value
        f_temp1[str(i)+'/magErr']=f[str(i)+'/magErr'].value
        
#        for saving only the linear+crts contributions
        
        # f[str(i)+'/MJDy']=MJDy
        # f[str(i)+'/y']=y
        # f[str(i)+'/magErry']=magErry
        
        f_temp1[str(i)+'/omegaSlope']=f[str(i)+'/omegaSlope'].value
        f_temp1[str(i)+'/sigma']=f[str(i)+'/sigma'].value
        f_temp1[str(i)+'/tau']=f[str(i)+'/tau'].value
        
#         f_temp2[str(i)+'/MJD']=f['/'+str(i)+'/MJD'].value
#         # f[str(i)+'/mag']=x
#         f_temp2[str(i)+'/mag']=f['/'+str(i)+'/mag'].value
#         # f_temp[str(i)+'/P']= f['/'+str(i)+'/P'].value
#         # f_temp[str(i)+'/A']=f[str(i)+'/A'].value
#         # # f[str(i)+'/phi']=phi
#         # f_temp[str(i)+'/t0']=f[str(i)+'/t0'].value
#         f_temp2[str(i)+'/magErr']=f[str(i)+'/magErr'].value
        
# #        for saving only the linear+crts contributions
        
#         # f[str(i)+'/MJDy']=MJDy
#         # f[str(i)+'/y']=y
#         # f[str(i)+'/magErry']=magErry
        
#         f_temp2[str(i)+'/omegaSlope']=f[str(i)+'/omegaSlope'].value
#         f_temp2[str(i)+'/sigma']=f[str(i)+'/sigma'].value
#         f_temp2[str(i)+'/tau']=f[str(i)+'/tau'].value
        
        # f_temp[str(i)+'/PLS']=f[str(i)+'/PLS'].value
        # f_temp[str(i)+'/maxpeak']=f[str(i)+'/maxpeak'].value
        # f_temp[str(i)+'/logf']=f[str(i)+'/logf'].value
        
    return None
    
# import h5py
separateCurves('C:/Users/noahk/pg1302-research/functions/10_L+C+A_rand_sig_tau_no_outlier.h5')

# if __name__ == "__main__":
#     separateCurves(sys.argv)
    