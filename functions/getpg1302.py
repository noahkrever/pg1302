# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 17:32:42 2020

@author: noahk
"""

# import SimulatePeriodograms
import numpy as np
from matplotlib import pyplot as plt
import h5py

data=np.loadtxt(r'C:\Users\noahk\pg1302-research\observation-data\pg1302_data_LINEAR+CRTS+ASASSN.txt')
   
MJD=data[:,0]
mag=data[:,1]
magErr=data[:,2]
flag = data[:,3]

SortedInd = np.argsort(MJD)
MJD = MJD[SortedInd]
mag = mag[SortedInd]
magErr = magErr[SortedInd]
flag = flag[SortedInd]

ind=np.where(flag==0)
y=mag[ind]
MJDy=MJD[ind]
magErry=magErr[ind]

MJD = np.delete(MJD,286)
mag = np.delete(mag,286)
magErr = np.delete(magErr,286)

from OptimalFrequencyGrid import OptimalFrequencyGrid
omega, omegaSlope = OptimalFrequencyGrid(MJD)
        
        
plt.errorbar(MJD, mag, yerr=magErr, fmt='o', ecolor='black')
plt.xlabel('MJD (Days)')
plt.ylabel('Magnitude')
plt.title("PG1302 L+C+A")
#       plt.legend()
plt.show()

plt.errorbar(MJDy, y, yerr=magErry, fmt='o', ecolor='black')
plt.xlabel('MJD (Days)')
plt.ylabel('Magnitude')
plt.title("PG1302 L+C")
#       plt.legend()
plt.show()

# f = h5py.File('pg1302_LINEAR+CRTS+ASASSN.h5','a')

f = h5py.File('pg1302_L+C+A_no_outlier.h5','a')

f[str(0)+'/MJD']=MJD
# f[str(i)+'/mag']=x
f[str(0)+'/mag']= mag
f[str(0)+'/magErr']=magErr

#        for saving only the linear+crts contributions

f[str(0)+'/MJDy']=MJDy
f[str(0)+'/y']=y
f[str(0)+'/magErry']=magErry

f[str(0)+'/omegaSlope']=omegaSlope
f[str(0)+'/sigma']=.11
f[str(0)+'/tau']=100

f.close()