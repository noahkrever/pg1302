# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 12:41:36 2019

@author: noahk
"""


# from matplotlib import pyplot as plt
# from matplotlib import cm
# import matplotlib.colors as mplc
# from mpl_toolkits.mplot3d import Axes3D
# from scipy.optimize import minimize
# import numpy as np
# from matplotlib import pyplot as plt
# from ChiSquareMinimization import drw_likelihood
# from readcurves import readcurve

# # from readcurves import readcurve
# import h5py
# import zipfile
# import sys

# filename = 'final_curve_pg1302_LC.h5'
# f=h5py.File(filename,"r")
# #import scipy.optimize
# # def lnprob(params,mag_DRW,magErr,time):
# #     return -Likelihood_drw_sinusoid(params[0],params[1],params[2], params[3], params[4], mag_DRW,magErr,time)

# def lnprob(sigma,tau,mag,magErr,MJD):
#     return -drw_likelihood(sigma, tau, mag, magErr, MJD)
# ##Manual Search with Suface Plots
# #
# #rho_i_list = []
# #rho_o_list = []
# #zf = zipfile.ZipFile('pg1302_L+C+A_10runs_randsig_randtau.zip', mode='w')
# # def findL(sigma,tau,A,period,phase,i,filename):
# #     info = readcurve(i,filename)
# #     mag = info[0]
# #     magErr = info[1]
# #     MJD = info[2]
# #     return Likelihood_drw_sinusoid(sigma,tau,A,period,phase,mag,magErr,MJD)

# # def findL(sigma,tau,A,period,phase,i,filename):
# #     info = readcurve(i,filename)
# #     mag = info[0]
# #     magErr = info[1]
# #     MJD = info[2]
# #     return Likelihood_drw_sinusoid(sigma,tau,A,period,phase,mag,magErr,MJD)
    
# def findL(sigma,tau,i,filename):
#     info = readcurve(i,filename)
#     mag = info[0]
#     magErr = info[1]
#     MJD = info[2]
#     return drw_likelihood(sigma,tau,mag,magErr,MJD)


# for i in f:
#     i_sigma = f['/'+str(i)+'/sigma']
#     i_sigma = i_sigma.value
#     i_tau = f['/'+str(i)+'/tau']
#     i_tau = i_tau.value
#     # i_P = f['/'+str(i)+'/P']
#     # i_P = i_P.value
#     # i_A = f['/'+str(i)+'/A']
#     # i_A = i_A.value
#     # i_phi = f['/'+str(i)+'/phi']
#     # i_phi = i_phi.value
    
#     print(i_sigma,i_tau)
# #    fig = plt.figure()
#     sigmas = np.linspace(0,1,num=50)
# #    taus = np.logspace(0,4.76,num=100) #ranges from 0 to log(maxtau), maxtau = t_exp * 10^1
#     taus = np.linspace(0,1000,num=50)
# #    ax = fig.add_subplot(111, projection='3d')
#     X, Y = np.meshgrid(sigmas, taus)
#     likelihoods = np.array([findL(sigmas, taus, i, filename) for sigmas,taus in zip(np.ravel(X), np.ravel(Y))])
#     Z = likelihoods.reshape(X.shape)
#     #print X
#     #print Y
#     #print Z
#     p1 = np.percentile(likelihoods, 3)
#     p2 = np.percentile(likelihoods, 16)
#     p3 = np.percentile(likelihoods, 32)
#     p4 = np.percentile(likelihoods, 50)
#     p5 = np.percentile(likelihoods, 68)
#     p6 = np.percentile(likelihoods, 84)
#     p7 = np.percentile(likelihoods, 90)
#     p8 = np.percentile(likelihoods, 95)
#     p9 = np.percentile(likelihoods, 97)
#     p10 = np.percentile(likelihoods, 99)
#     p11 = np.percentile(likelihoods, 99.9)
#     percentile = np.array([p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11])
# #    surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
# #                           linewidth=0, antialiased=False)
# #    
# #    ax.set_xlabel('Sigma')
# #    ax.set_ylabel('Tau')
# #    ax.set_zlabel('Likelihood')
# #    fig.colorbar(surf)
# #    
# #    plt.show()
    
#     fig, ax0 = plt.subplots()
#     Colors0 = ['white', 'seashell', 'mistyrose', 'peachpuff', 'sandybrown', 'lightcoral', 'indianred', 'brown', 'darkred', 'skyblue']
#     cmap0= mplc.ListedColormap(Colors0)
#     levels0 = percentile
#     cmap0.set_over('dodgerblue')
#     p0 = ax0.contourf(X, Y, Z, cmap=None,levels=percentile[2:],colors=Colors0,extend='both',vmin=(Z).min(), vmax=(Z).max())
# #    CS0 = ax0.contour(X, Y, Z,levels=levels0,colors='black',linewidths=2)
#     fig.colorbar(p0, ax=ax0, extend='both')
#     plt.scatter(i_sigma,i_tau, marker = 0, color = 'yellow')
#     # plt.savefig('pg1302_L+C+A_10runs_randsig_randtau' + str(i) + '.png')
#     # zf.write('pg1302_L+C+A_10runs_randsig_randtau' + str(i) + '.png')
#     plt.show()
#     maxL = np.amax(Z)
#     result = np.where(Z==maxL)
#     msigma = float(X[result])
#     mtau = float(Y[result])
#     print(msigma, mtau)
#     print(str(i))
#     print("L max of " + str(maxL) + " is at the pair " + msigma + "," + mtau)
# #
#
##nll = lambda *args: -drw_likelihood(*args)
#
##Scipy Optimization
##    msigmas = []
##    mtaus = []
##    info = readcurve(i,filename)
##    mag = info[0]
##    magErr = info[1]
##    MJD = info[2]
##    res = scipy.optimize.minimize(lnprob, (i_sigma,i_tau),args=(mag,magErr,MJD))
##    msigma = res.x[0]
##    mtau = res.x[1]
##    print res.x
##    msigmas.append(msigma)
##    mtaus.append(mtau)
#    F2='pg1302_L+C+A_10runs_randsig_randtau_info.h5'
#    fi = h5py.File(F2, "a")
#    fi[str(i)+'/inputsigma']=i_sigma
#    fi[str(i)+'/inputtau']=i_tau
#    fi[str(i)+'/outputsigma']=msigma
#    fi[str(i)+'/outputtau']=mtau
#    fi[str(i)+'/likelihood']=maxL #note, this may not be exact because it is from the grid
#
#    fi.close()
    
#    t_exp = 2920
#    rho_i = i_tau / t_exp
#    rho_o = mtau / t_exp
#    rho_i_list.append(rho_i)
#    rho_o_list.append(rho_o)
#    
#plt.scatter(rho_i_list, rho_o_list)
#plt.xlabel('rho_in')
#plt.ylabel('rho_out')
#plt.xscale('log')
#plt.yscale('log')
#plt.title('Input vs. Output Rho for Kozlowski OGLE')
#plt.show()

    
def drwsin_MCMC(filename,Filename2):
    from ChiSquareMinimization import Likelihood_drw_sinusoid
    from matplotlib import pyplot as plt
    import h5py
    import numpy as np
    #MCMC
    # filename='100_DRW+Sinusoid.h5'
    f=h5py.File(filename,"r")
    #i_sigmas = []
    #i_taus = []
    #o_sigmas = []
    #o_taus = []
    #s_upper = []
    #s_lower = []
    #t_upper = []
    #t_lower = []
    def lnprior(theta):
        # sigma, tau, A, P, phase = theta
        # if 0.00001 < sigma < 1 and 1 < tau < 1000 and 0 < A < 1 and 0 < phase < 2*np.pi and 30 < P < 5500:
        sigma, log_tau, A, P, t0 = theta
        if 0.00001 < sigma < 1 and 0 < log_tau < 3 and 0 < A < 1 and 1 < t0 < 1000 and 30 < P < 5500:
            return 0.0
        return -np.inf
    
    def lnprob(theta, mag, magErr, MJD): 
        sigma, log_tau, A, P, t0 = theta
        lp = lnprior(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + Likelihood_drw_sinusoid(sigma, log_tau, A, P, t0, mag, magErr, MJD)
    
    for i in f:
        MJD=f['/'+str(i)+'/MJD']
        MJD=MJD.value
        mag= f['/'+str(i)+'/mag']
        mag=mag.value
        magErr= f['/'+str(i)+'/magErr']
        magErr=magErr.value
        
        i_sigma = f['/'+str(i)+'/sigma']
        i_sigma = i_sigma.value
        i_tau = f['/'+str(i)+'/tau']
        i_tau = i_tau.value
        # i_P = f['/'+str(i)+'/P']
        # i_P = i_P.value
        # i_A = f['/'+str(i)+'/A']
        # i_A = i_A.value
        # # i_phi = f['/'+str(i)+'/phi']
        # # i_phi = i_phi.value
        # i_t0 = f['/'+str(i)+'/t0']
        # i_t0 = i_t0.value
        print(i)
    #    print i_sigma
    #    print i_tau
    #    i_sigmas.append(i_sigma)
    #    i_taus.append(i_tau)
        
    #    data=np.loadtxt('C:/Python27/Data/pg1302_data_LINEAR+CRTS+ASASSN.txt')
    #    days = data[:,0]
    #    flag = data[:,3]
    #    
    #    SortedInd = np.argsort(days)
    #    days = days[SortedInd]
    #    flag = flag[SortedInd]
    #    
    #    ind=np.where(flag==0)
    #    y=mag[ind]
    #    MJDy=MJD[ind]
    #    magErry=magErr[ind]
        
        # import scipy.optimize
        #
    #    res = scipy.optimize.minimize(lnprob, (i_sigma,i_tau,i_A,i_P,i_phi),args=(mag,0,MJD))
    #    print(res.x)
    #    
    #    res.x=[i_sigma,i_tau,i_A,i_P,i_phi]
        optimal = [.11,np.log10(100),.1,1825,100]
        import emcee
        ndim, nwalkers = 5, 100
        
        pos = [optimal + 1e-4*np.random.randn(ndim) for a in range(nwalkers)]
        
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(mag,magErr,MJD))
        
        sampler.run_mcmc(pos, 1000)
        # sampler.run_mcmc(pos, 500)
        
        samples = sampler.chain[:, 100:, :].reshape((-1, ndim))
        
        chains = sampler.chain
        
        all_lnprobs = sampler.lnprobability
        
        # import corner
        # # fig = corner.corner(samples, labels=["$sigma$", "$tau$","$Amplitude$", "$Period$", "$Phase$"],
        #                     # truths=[i_sigma, i_tau, i_A, i_P, i_phi])
        # fig = corner.corner(samples, labels=["$sigma$", "$tau$","$Amplitude$", "$Period$", "$t_0$"],
        #                     truths=[i_sigma, i_tau, i_A, i_P, i_t0])
        # # fig.savefig("100runs_triangle" + str(i) + ".png")
        # fig.savefig("drwsin_triangle" + str(i) + ".png")
        
        # plt.show()
    #    
        # s_mcmc, t_mcmc, A_mcmc, P_mcmc, phase_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
        # print (s_mcmc[0], t_mcmc[0], A_mcmc[0], P_mcmc[0], phase_mcmc[0])
        s_mcmc, t_mcmc, A_mcmc, P_mcmc, t0_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
        print (s_mcmc[0], t_mcmc[0], A_mcmc[0], P_mcmc[0], t0_mcmc[0])
    # #    
    # #    s_upper.append(s_mcmc[1])
    # #    s_lower.append(s_mcmc[2])
    # #    t_upper.append(t_mcmc[1])
    # #    t_lower.append(t_mcmc[2])
    # #    
    # #    i_sigmas.append(i_sigma)
    # #    o_sigmas.append(s_mcmc[0])
    # #    
    # #    i_taus.append(i_tau)
    # #    o_taus.append(t_mcmc[0])
    # #    
    # #    for b in np.arange(0,100):
    # #        plt.plot(chains[b,:,0])
    # #        plt.title("Sigma Chains")
    # #
    # #    plt.savefig("schains" + str(i) + ".png")
    # #    plt.show()
    # #    
    # #    for c in np.arange(0,100):
    # #        plt.plot(chains[c,:,1])
    # #        plt.title("Tau Chains")
    # #    
    # #    plt.savefig("tchains" + str(i) + ".png")
    # #    plt.show()
    # #    
        # Filename2='100_DRW+Sinusoid_info_with_lnprobs.h5'
        fi = h5py.File(Filename2, "a")
    #    fi[str(i)+'/inputsigma']=i_sigma
    #    fi[str(i)+'/inputtau']=i_tau
    #    fi[str(i)+'/sigma_mean']=s_mcmc[0]
    #    fi[str(i)+'/tau_mean']=t_mcmc[0]
        fi[str(i)+'/lnprobs']=all_lnprobs
        fi[str(i)+'/chains']=chains
    
    return None

def drw_MCMC(filename,Filename2):
    from ChiSquareMinimization import drw_likelihood
    from matplotlib import pyplot as plt
    import numpy as np
    import h5py
    
    f=h5py.File(filename,"r")
    
    def lnprior(theta):
        # sigma, tau, A, P, phase = theta
        # if 0.00001 < sigma < 1 and 1 < tau < 1000 and 0 < A < 1 and 0 < phase < 2*np.pi and 30 < P < 5500:
        sigma, log_tau = theta
        if 0.00001 < sigma < 1 and 0 < log_tau < 3:
            return 0.0
        return -np.inf
    
    def lnprob(theta, mag, magErr, MJD): 
        sigma, log_tau = theta
        lp = lnprior(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + drw_likelihood(sigma, log_tau, mag, magErr, MJD)
    
    for i in f:
        MJD=f['/'+str(i)+'/MJD']
        MJD=MJD.value
        mag= f['/'+str(i)+'/mag']
        mag=mag.value
        magErr= f['/'+str(i)+'/magErr']
        magErr=magErr.value
        
        i_sigma = f['/'+str(i)+'/sigma']
        i_sigma = i_sigma.value
        i_tau = f['/'+str(i)+'/tau']
        i_tau = i_tau.value
        print(i)
        optimal = [i_sigma,np.log10(i_tau)]
        import emcee
        ndim, nwalkers = 2, 100
        
        pos = [optimal + 1e-4*np.random.randn(ndim) for a in range(nwalkers)]
        
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(mag,magErr,MJD))
        
        sampler.run_mcmc(pos, 1000)

        samples = sampler.chain[:, 100:, :].reshape((-1, ndim))
        
        chains = sampler.chain
        
        all_lnprobs = sampler.lnprobability
        
        s_mcmc, t_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
        print (s_mcmc[0], t_mcmc[0])
        fi = h5py.File(Filename2, "a")
        fi[str(i)+'/lnprobs']=all_lnprobs
        fi[str(i)+'/chains']=chains
    
    return None

    #     fi.close()
    
    #s_errors = np.array([s_lower, s_upper])
    #t_errors = np.array([t_lower, t_upper])
    #
    #fig, ax = plt.subplots()
    #ax.errorbar(i_sigmas, o_sigmas, yerr=s_errors, fmt = 'o', c='black')
    #lims = [
    #    np.min([ax.get_xlim(), ax.get_ylim()]),
    #    np.max([ax.get_xlim(), ax.get_ylim()]),
    #]
    #
    #ax.plot(lims, lims, '--', alpha=0.75, zorder=0)
    #ax.set_aspect('equal')
    #ax.set_xlim(lims)
    #ax.set_ylim(lims)
    #plt.xlabel('Input Sigma')
    #plt.ylabel('Output Sigma Mean')
    #plt.title('25 Random Sigma MCMC Results')
    #
    #fig, ax = plt.subplots()
    #ax.errorbar(i_taus, o_taus, yerr = t_errors, fmt = 'o', c='black')
    #lims = [
    #    np.min([ax.get_xlim(), ax.get_ylim()]),
    #    np.max([ax.get_xlim(), ax.get_ylim()]),
    #]
    #
    #ax.plot(lims, lims, '--', alpha=0.75, zorder=0)
    #ax.set_aspect('equal')
    #ax.set_xlim(lims)
    #ax.set_ylim(lims)
    #plt.xlabel('Input Tau')
    #plt.ylabel('Output Tau Mean')
    #plt.title('25 Random Tau MCMC Results')



    
    