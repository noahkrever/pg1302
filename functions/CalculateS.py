# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 17:55:47 2019

@author: noahk
"""

#CalculateS

from readcurves import readcurve
import numpy as np
import h5py
from matplotlib import pyplot as plt

File1='500_WS_PSD2_t200_c20_m14.h5'
f1 = h5py.File(File1, "r")
File2='500_WS_PSD2_t200_m0_submean.h5'
f2 = h5py.File(File2, "r")
File3='500_WS_PSD2_t200_m0_3_submean.h5'
f3 = h5py.File(File2, "r")
s0 = np.array([])
s1 = np.array([])
s2 = np.array([])
s3 = np.array([])
s4 = np.array([])
s5 = np.array([])
s6 = np.array([])
s7 = np.array([])
s8 = np.array([])
s9 = np.array([])
_2s0 = np.array([])
_2s1 = np.array([])
_2s2 = np.array([])
_2s3 = np.array([])
_2s4 = np.array([])
_2s5 = np.array([])
_2s6 = np.array([])
_2s7 = np.array([])
_2s8 = np.array([])
_2s9 = np.array([])
_3s0 = np.array([])
_3s1 = np.array([])
_3s2 = np.array([])
_3s3 = np.array([])
_3s4 = np.array([])
_3s5 = np.array([])
_3s6 = np.array([])
_3s7 = np.array([])
_3s8 = np.array([])
_3s9 = np.array([])
for i in f1:
    info = readcurve(i,File1)
    mag = info[0]
    magErr = info[1]
    MJD = info[2]
    for b in range(10):
        products = np.array([])
        for a in range(1,2001):
            if a+b < 2001:
#                prod = abs(mag[a] - np.mean(mag)) * abs(mag[a+b] - np.mean(mag))
                prod = (mag[a] - np.mean(mag)) * (mag[a+b] - np.mean(mag))
                products = np.append(products,[prod])
        S = np.mean(products)
#        S = np.median(products)

        if b == 0:
            s0 = np.append(s0,[S])
        elif b == 1:
            s1 = np.append(s1,[S])
        elif b == 2:
            s2 = np.append(s2,[S])
        elif b == 3:
            s3 = np.append(s3,[S])
        elif b == 4:
            s4 = np.append(s4,[S])
        elif b == 5:
            s5 = np.append(s5,[S])
        elif b == 6:
            s6 = np.append(s6,[S])
        elif b == 7:
            s7 = np.append(s7,[S])
        elif b == 8:
            s8 = np.append(s8,[S])
        elif b == 9:
            s9 = np.append(s9,[S])
#for j in f2:
#    info2 = readcurve(j,File2)
#    mag = info2[0]
#    magErr = info2[1]
#    MJD = info2[2]
#    for b in range(10):
#        products2 = np.array([])
#        for a in range(1,4001):
#            if a+b < 4001:
#                prod2 = abs(mag[a] - 14) * abs(mag[a+b] - 14)
#                products2 = np.append(products2,[prod2])
##        S = np.average(products2)
#        S = np.median(products2)
#
#        if b == 0:
#            _2s0 = np.append(_2s0,[S])
#        elif b == 1:
#            _2s1 = np.append(_2s1,[S])
#        elif b == 2:
#            _2s2 = np.append(_2s2,[S])
#        elif b == 3:
#            _2s3 = np.append(_2s3,[S])
#        elif b == 4:
#            _2s4 = np.append(_2s4,[S])
#        elif b == 5:
#            _2s5 = np.append(_2s5,[S])
#        elif b == 6:
#            _2s6 = np.append(_2s6,[S])
#        elif b == 7:
#            _2s7 = np.append(_2s7,[S])
#        elif b == 8:
#            _2s8 = np.append(_2s8,[S])
#        elif b == 9:
#            _2s9 = np.append(_2s9,[S])
#for k in f3:
#    info3 = readcurve(k,File3)
#    mag = info3[0]
#    magErr = info3[1]
#    MJD = info3[2]
#    for b in range(10):
#        products3 = np.array([])
#        for a in range(1,2001):
#            if a+b < 2001:
#                prod3 = abs(mag[a] - 14) * abs(mag[a+b] - 14)
#                products3 = np.append(products3,[prod3])
##        S = np.average(products3)
#        S = np.median(products3)
#
#        if b == 0:
#            _3s0 = np.append(_3s0,[S])
#        elif b == 1:
#            _3s1 = np.append(_3s1,[S])
#        elif b == 2:
#            _3s2 = np.append(_3s2,[S])
#        elif b == 3:
#            _3s3 = np.append(_3s3,[S])
#        elif b == 4:
#            _3s4 = np.append(_3s4,[S])
#        elif b == 5:
#            _3s5 = np.append(_3s5,[S])
#        elif b == 6:
#            _3s6 = np.append(_3s6,[S])
#        elif b == 7:
#            _3s7 = np.append(_3s7,[S])
#        elif b == 8:
#            _3s8 = np.append(_3s8,[S])
#        elif b == 9:
#            _3s9 = np.append(_3s9,[S])
            
                
true_20 = np.array([.0121, .0109, .0099, .00896, .00811, .00734, .00664, .006, .00543, .00492])
true_10 = [.0121, .0115, .0109, .0104, .0099, .00942, .00896, .00853, .00811, .0077]
true_5 = [.0121, .0118, .0115, .0112, .0109, .01067, .0104, .01016, .0099, .00966]
average1 = np.array([np.average(s0), np.average(s1), np.average(s2), np.average(s3), np.average(s4), np.average(s5), np.average(s6), np.average(s7), np.average(s8), np.average(s9)])
average2 = [np.average(_2s0), np.average(_2s1), np.average(_2s2), np.average(_2s3), np.average(_2s4), np.average(_2s5), np.average(_2s6), np.average(_2s7), np.average(_2s8), np.average(_2s9)]
average3 = [np.average(_3s0), np.average(_3s1), np.average(_3s2), np.average(_3s3), np.average(_3s4), np.average(_3s5), np.average(_3s6), np.average(_3s7), np.average(_3s8), np.average(_3s9)]

#print products
#print mag
fig, ax = plt.subplots()
plt.scatter(true_20, average1, color= 'green', label="20 Day Cadence, Tau = 200d")
#plt.scatter(true_10, average3, color= 'blue', label="10 Day Cadence, Tau = 200d")
#plt.scatter(true_5, average2, color= 'red', label="5 Day Cadence, Tau = 200d")
plt.legend()
lims = [
    np.min([ax.get_xlim(), ax.get_ylim()]),
    np.max([ax.get_xlim(), ax.get_ylim()]),
]

ax.plot(lims, lims, '--', alpha=0.75, zorder=0, color='black')
ax.set_aspect('equal')
ax.set_xlim(lims)
ax.set_ylim(lims)
plt.xlabel('sigma^2 * exp(-deltat/tau)')
plt.ylabel('Average S_ij')
plt.xlim(0.002, .015)
plt.ylim(0.002, .015)
plt.title('Expected S_IJ Value vs. Means for 500 Well Sampled Curves With Subtracted Means')
plt.show()

#plt.errorbar(np.arange(0,len(s0)), s0, fmt = '-o')
#plt.title('S_11')
#plt.show()
#
#plt.errorbar(np.arange(0,len(s1)), s1, fmt = '-o')
#plt.title('S_12')
#plt.show()
#
#plt.errorbar(np.arange(0,len(s2)), s2, fmt = '-o')
#plt.title('S_24')
#plt.show()
#
#plt.errorbar(np.arange(0,len(s3)), s3, fmt = '-o')
#plt.title('S_36')
#plt.show()
#
#plt.errorbar(np.arange(0,len(s4)), s4, fmt = '-o')
#plt.title('S_48')
#plt.show()
#
#plt.errorbar(np.arange(0,len(s5)), s5, fmt = '-o')
#plt.title('S_510')
#plt.show()
#
#plt.errorbar(np.arange(0,len(s6)), s6, fmt = '-o')
#plt.title('S_612')
#plt.show()
#
#plt.errorbar(np.arange(0,len(s7)), s7, fmt = '-o')
#plt.title('S_714')
#plt.show()
#
#plt.errorbar(np.arange(0,len(s8)), s8, fmt = '-o')
#plt.title('S_816')
#plt.show()
#
#plt.errorbar(np.arange(0,len(s9)), s9, fmt = '-o')
#plt.title('S_918')
#plt.show()

s_division = average1 / true_20
dt = np.array([0, 20, 40, 60, 80, 100, 120, 140, 160, 180])

plt.scatter(dt, s_division)
plt.xlabel('dt')
plt.ylabel('<S_IJ> / Expected S_IJ')
plt.title('dt Versus <S_IJ> / Expected S_IJ for 20 Day Cadence, Tau = 200')



