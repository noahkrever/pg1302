# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 16:14:48 2019

@author: noahk
"""

import numpy as np
from matplotlib import pyplot as plt
import h5py

filename1 = 'pg1302_L+C+A_10runs_randsig_randtau_info.h5'
f1=h5py.File(filename1,"r")

filename2 = 'Kozlowski_OGLE_info.h5'
f2=h5py.File(filename2,"r")
sdss_i_rhos = []
sdss_o_rhos = []
ogle_i_rhos = []
ogle_o_rhos = []
input_sigmas = []
output_sigmas = []
t_exp = 5780
for i in f1:
    i_sigma = f1['/'+str(i)+'/inputsigma']
    i_sigma = i_sigma.value
    i_tau = f1['/'+str(i)+'/inputtau']
    i_tau = i_tau.value
    o_sigma = f1['/'+str(i)+'/outputsigma']
    o_sigma = o_sigma.value
    o_tau = f1['/'+str(i)+'/outputtau']
    o_tau = o_tau.value
    rho_i = i_tau / t_exp
    rho_o = o_tau / t_exp
    rho_i_new = np.log10(rho_i)
    rho_o_new = np.log10(rho_o)
    ogle_i_rhos.append(rho_i_new)
    ogle_o_rhos.append(rho_o_new)
    input_sigmas.append(i_sigma)
    output_sigmas.append(o_sigma)
    print(i_sigma, i_tau)
    print(o_sigma,o_tau)
fig, ax = plt.subplots()
ax.scatter(ogle_i_rhos, ogle_o_rhos)
lims = [
    np.min([ax.get_xlim(), ax.get_ylim()]),
    np.max([ax.get_xlim(), ax.get_ylim()]),
]

ax.plot(lims, lims, '--', alpha=0.75, zorder=0)

ax.set_xlim(lims)
ax.set_ylim(lims)

plt.axvline(-1)
plt.axhline(-1)
plt.xlabel('Log Input Rho')
plt.ylabel('Log Output Rho')
plt.title('PG1302 LINEAR+CRTS+ASASSN With Random Tau, Random Sigma')
plt.show()

fig, ax = plt.subplots()
ax.scatter(input_sigmas, output_sigmas)
lims = [
    np.min([ax.get_xlim(), ax.get_ylim()]),
    np.max([ax.get_xlim(), ax.get_ylim()]),
]

ax.plot(lims, lims, '--', alpha=0.75, zorder=0)

ax.set_xlim(lims)
ax.set_ylim(lims)

plt.axvline(-1)
plt.axhline(-1)
plt.xlabel('Input Sigma')
plt.ylabel('Output Sigma')
plt.title('PG1302 LINEAR+CRTS+ASASSN Input vs. Output Sigmas')
plt.show()
