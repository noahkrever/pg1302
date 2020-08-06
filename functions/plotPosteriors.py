# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 16:41:45 2020

@author: noahk
"""
# import SimulatePeriodograms

import h5py
from matplotlib import pyplot as plt
import SimulatePeriodograms

# SimulatePeriodograms.Simulate(20,'20_LC_DRWSIN_pg1302_test')

filename = '20_LCA_DRW_pg1302_test'

f=h5py.File(filename,"r")

filename2 = 'final_curve_pg1302_LCA.h5'

f2=h5py.File(filename2,"r")

for i in f:
    mag = f['/'+str(i)+'/mag'].value
    # mag=mag[0:202]
    MJD = f['/'+str(i)+'/MJD'].value
    # MJD = MJD[0:202]
    plt.plot(MJD,mag,color='gray')
    
    
mag = f2['/'+str(0)+'/mag'].value - 14.8
MJD = f2['/'+str(0)+'/MJD'].value
magErr = f2['/'+str(0)+'/magErr'].value

plt.errorbar(MJD,mag,yerr=magErr,ecolor='black',color='red',fmt='o')


