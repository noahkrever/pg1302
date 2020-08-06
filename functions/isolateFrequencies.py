# -*- coding: utf-8 -*-
"""
Created on Sat Jun 13 19:48:38 2020

@author: noahk
"""
def isolateFrequencies(curve_file):
    import readcurves
    import h5py
    import numpy as np
    o = h5py.File("final_set.h5", "a")
    f=h5py.File(curve_file,"r")
    c = 0
    for i in f:
        mag,magErr,MJD,max_peak,y,magErry,MJDy,PLS,logf,omegaSlope,frequency = readcurves.readcurve(i,curve_file)
        
        m = np.argmax(PLS)
        max_f = frequency[m]
        
        if 5.02995849e-9 < max_f and max_f < 7.08920325e-9:
            print(max_f)
            print(c)
            o[str(c)+'/mag'] = mag
            o[str(c)+'/magErr'] = magErr
            o[str(c)+'/MJD'] = MJD
            o[str(c)+'/maxpeak'] = max_peak
            o[str(c)+'/y'] = y
            o[str(c)+'/magErry'] = magErry
            o[str(c)+'/MJDy'] = MJDy
            o[str(c)+'/PLS'] = PLS
            o[str(c)+'/logf'] = logf
            o[str(c)+'/omegaSlope'] = omegaSlope
            o[str(c)+'/frequency'] = frequency
            o[str(c)+'/sigma'] = .11
            o[str(c)+'/tau'] = 100
            
            c += 1
    return None

import h5py
import readcurves
o=h5py.File('all_significant.h5',"a")

c1 = 0
c2 = 0
c3 = 0

for i in [1,2,3,4,5,7,8,9,10]:
    filename=r'C:\Users\noahk\pg1302-research\cluster-data\significant_1000000_run_' + str(i) + '.h5'
    f=h5py.File(filename,"r")
    
    for j in f:
        c1+= 1
        mag,magErr,MJD,max_peak,y,magErry,MJDy,PLS,logf,omegaSlope,frequency = readcurves.readcurve(j,filename)
        
        o[str(j)+'/mag'] = mag
        o[str(j)+'/magErr'] = magErr
        o[str(j)+'/MJD'] = MJD
        o[str(j)+'/maxpeak'] = max_peak
        o[str(j)+'/y'] = y
        o[str(j)+'/magErry'] = magErry
        o[str(j)+'/MJDy'] = MJDy
        o[str(j)+'/PLS'] = PLS
        o[str(j)+'/logf'] = logf
        o[str(j)+'/omegaSlope'] = omegaSlope
        o[str(j)+'/frequency'] = frequency
for i in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,28,29,30]:
    filename=r'C:\Users\noahk\pg1302-research\cluster-data\significant_1000000_run_' + str(i) + '_2.h5'
    f=h5py.File(filename,"r")
    for j in f:
        c2+= 1
        mag,magErr,MJD,max_peak,y,magErry,MJDy,PLS,logf,omegaSlope,frequency = readcurves.readcurve(j,filename)
        
        o[str(int(j) * 47)+'/mag'] = mag
        o[str(int(j) * 47)+'/magErr'] = magErr
        o[str(int(j) * 47)+'/MJD'] = MJD
        o[str(int(j) * 47)+'/maxpeak'] = max_peak
        o[str(int(j) * 47)+'/y'] = y
        o[str(int(j) * 47)+'/magErry'] = magErry
        o[str(int(j) * 47)+'/MJDy'] = MJDy
        o[str(int(j) * 47)+'/PLS'] = PLS
        o[str(int(j) * 47)+'/logf'] = logf
        o[str(int(j) * 47)+'/omegaSlope'] = omegaSlope
        o[str(int(j) * 47)+'/frequency'] = frequency
        
for i in [1,2,3,4,5,6,7,8,9,10]:
    filename=r'C:\Users\noahk\pg1302-research\cluster-data\significant_1000000_run_' + str(i) + '_3.h5'
    f=h5py.File(filename,"r")
    for j in f:
        c3+= 1
        mag,magErr,MJD,max_peak,y,magErry,MJDy,PLS,logf,omegaSlope,frequency = readcurves.readcurve(j,filename)
        
        o[str(int(j) * 67)+'/mag'] = mag
        o[str(int(j) * 67)+'/magErr'] = magErr
        o[str(int(j) * 67)+'/MJD'] = MJD
        o[str(int(j) * 67)+'/maxpeak'] = max_peak
        o[str(int(j) * 67)+'/y'] = y
        o[str(int(j) * 67)+'/magErry'] = magErry
        o[str(int(j) * 67)+'/MJDy'] = MJDy
        o[str(int(j) * 67)+'/PLS'] = PLS
        o[str(int(j) * 67)+'/logf'] = logf
        o[str(int(j) * 67)+'/omegaSlope'] = omegaSlope
        o[str(int(j) * 67)+'/frequency'] = frequency
        
isolateFrequencies("all_significant.h5")
# print(c1,c2,c3)