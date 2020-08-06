# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 17:55:49 2020

@author: noahk
"""
import numpy as np
from BinLightCurve import BinLightCurveFiveDayInterval as bc
from matplotlib import pyplot as plt
import h5py

data=np.loadtxt(r"C:\Users\noahk\pg1302-research\observation-data\pg1302_data_LINEAR+CRTS+ASASSN.txt")
MJD=data[:,0]
mag =data[:,1]
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

ind2=np.where(flag==1)
a=mag[ind2] - .03489723624
MJDa=MJD[ind2]
magErra=magErr[ind2]
        
LC_MJD, LC_mag, LC_magErr = bc(MJDy,y,magErry,180)
A_MJD, A_mag, A_magErr = bc(MJDa,a,magErra,100)

bmjd = np.concatenate((LC_MJD, A_MJD))

ind = np.argsort(bmjd)

bmjd = bmjd[ind]

bmag = np.concatenate((LC_mag, A_mag))
bmag = bmag[ind]

bmagErr = np.concatenate((magErry, magErra))
bmagErr = bmag[ind]

f = h5py.File('pg1302_L+C+A_Binned', "a")
f[str(0)+'/MJD']=bmjd
f[str(0)+'/mag']= bmag
f[str(0)+'/magErr']=bmagErr
# f[str(0)+'/omegaSlope']=omegaSlope
f[str(0)+'/sigma']=.11
f[str(0)+'/tau']=100

plt.errorbar(MJD, mag, yerr=magErr, fmt='o', color='blue',zorder=1)
plt.scatter(LC_MJD, LC_mag, color='red', s=150,zorder=2, edgecolors=['black'])
plt.scatter(A_MJD, A_mag, color='lawngreen', s=150,zorder=2, edgecolors=['black'])
plt.title('PG1302 Binned (L+C 180 and A 100)')
plt.show()

plt.errorbar(bmjd,bmag,yerr=bmagErr,fmt='o')
plt.show()