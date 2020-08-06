import numpy as np
from LightCurves import Sampling_Ideal_Survey, DRW_LightCurve,Beaming_LightCurve

QSOId=1

time, magErr=Sampling_Ideal_Survey(QSOId)

sigma=0.15
tau=150
mag_DRW=DRW_LightCurve(sigma,tau,time)


DatFile='LightCurvesDat_QSO_' + str(QSOId) + '.dat'

f = open(DatFile, 'w')
for kk in range(len(time)):
    f.write("%6.6f %6.6f %6.6f\n"%(time[kk] - np.min(time),
                                   mag_DRW[kk], 1.0*magErr[kk]))
f.close()



from ChiSquareMinimization import drw_likelihood

def lnprob(params,mag_DRW,magErr,time):
    return -drw_likelihood(params[0],params[1],mag_DRW,magErr,time)

nll = lambda *args: -drw_likelihood(*args)

import scipy.optimize
#
res = scipy.optimize.minimize(nll, (0.1,100),args=(mag_DRW,magErr,time))
print res.x


import emcee
ndim, nwalkers = 2, 100
pos = [res.x + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(mag_DRW,magErr,time))

sampler.run_mcmc(pos, 500)
