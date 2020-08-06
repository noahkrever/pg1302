# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 15:49:19 2019

@author: noahk
"""

# def readMCMC(filename):
def readMCMC(filename,k):
    import h5py
    from matplotlib import pyplot as plt
    import numpy as np
    from scipy.stats import scoreatpercentile as scoretpercentile
    #filename='25_Runs_MCMC_Info_2.h5'
    f=h5py.File(filename,"r") #MCMC INFO
    # f2=h5py.File(filename2,"r")
    # i_sigmas = []
    # o_sigmas = []
    # i_taus = []
    # o_taus = []
    # s_upper = []
    # s_lower = []
    # t_upper = []
    # t_lower = []
    mcmc_sigmas = []
    mcmc_taus = []
    mcmc_amps = []
    mcmc_ps = []
    # mcmc_phis = []
    mcmc_t0s = []
    
    # bad = []
    
    best_di = {}
    probs_di = {}
    
    for i in f:
    #     # i_sigma=f2['/'+str(i)+'/inputsigma']
    #     # i_sigma=i_sigma.value
        
    #     # i_tau=f2['/'+str(i)+'/inputtau']
    #     # i_tau=i_tau.value
        
    #     i_sigma = .11
        
    #     i_tau = 100
        
    #     i_amp = .1
        
    #     i_p = 1825
        
    #     i_phi = 0
        
    #     # o_sigma=f2['/'+str(i)+'/outputsigma']
    #     # o_sigma = o_sigma.value
    #     # o_tau = f2['/'+str(i)+'/outputtau']
    #     # o_tau = o_tau.value
    # #    sigma_mean= f2['/'+str(i)+'/sigma_mean']
    # #    sigma_mean=sigma_mean.value
    # #    tau_mean = f2['/'+str(i)+'/tau_mean']
    # #    tau_mean=tau_mean.value
        lnprobs = f['/'+str(i)+'/lnprobs']
        lnprobs = lnprobs.value
        chains = f['/'+str(i)+'/chains']
        chains=chains.value
        
    #     # print(np.shape(chains))
    #     # print(np.shape(lnprobs))
        
        samples = chains[:, 300:, :].reshape((-1, k))   #-1, 2 if just sig tau
        probs = lnprobs[:, 300:].reshape(-1)
        

        s_chains = samples[:,0]
        t_chains = samples[:,1]
            
        if k == 5:
            amp_chains = samples[:,2]
            p_chains = samples[:,3]
            t0_chains = samples[:,4]
            
        
        # import corner
        # fig = corner.corner(samples, labels=["$sigma$", "$tau$"])
        
        # plt.show()

        # phi_chains = samples[:,4]
        
    
        # print(np.shape(lnprobs)) # this is (100, 500)
        # print(probs)
        # print(probs.shape)
        # print(s_chains.shape)
        b = np.argmax(probs)
    #     # coords = np.unravel_index(lnprobs.argmax(), lnprobs.shape) # returns best indices
    #     # print(coords)
        # x = coords[0]
        # y = coords[1]
        
        # best = [s_chains[b], t_chains[b], amp_chains[b], p_chains[b], t0_chains[b]]
        if k == 2:
            best = [s_chains[b], t_chains[b]]
            s_mcmc, t_mcmc = map(lambda v: (v[0], v[1], v[2]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
            print(s_mcmc, t_mcmc)
        else:
            best = [s_chains[b], t_chains[b], amp_chains[b], p_chains[b], t0_chains[b]]
            s_mcmc, t_mcmc, a_mcmc, p_mcmc, phi_mcmc = map(lambda v: (v[0], v[1], v[2]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
        # print(best)
        # mcmc_sigmas.append(best[0])
        # mcmc_taus.append(best[1])
        # mcmc_amps.append(best[2])
        # mcmc_ps.append(best[3])
        # # mcmc_phis.append(best[4])
        # mcmc_t0s.append(best[4])
        # s_mcmc, t_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
        
        # print("S percentiles : " + str(s_mcmc))
        # print("T percentiles : " + str(t_mcmc))
        probs_di[int(i)] = probs[b]
        best_di[int(i)] = best
        
        # import corner
        
        # # fig = corner.corner(samples, labels=["$sigma$", "$tau$","$Amplitude$", "$Period$", "$Phase$"],
        #                     # truths=[i_sigma, i_tau, i_A, i_P, i_phi])
        # fig = corner.corner(samples, labels=["$sigma$", "$tau$","$Amplitude$", "$Period$", "$t_0$"],
        #                     truths=[best[0], best[1], best[2], best[3], best[4]], quantiles=[.25,.50,.75])
        # # fig.savefig("100runs_triangle" + str(i) + ".png")
        # # fig.savefig("drwsin_triangle" + str(i) + ".png")
        # # 
        # fig.suptitle("PG1302 L+C DRWSIN")
        # plt.show()
        
    #     if best[0] < .08 or best[0] > .14:
    #         print("sigma, curve " + str(i) + " " + str(best))
    #         bad.append(i)
    #     if best[1] < 70 or best[1] > 130:
    #         print("tau, curve " + str(i) + " " + str(best))
    #         bad.append(i)
    #     if best[2] < .07 or best[2] > .13:
    #         print("amp, curve " + str(i) + " " + str(best))
    #         bad.append(i)
    #     if best[3] < 1500 or best[3] > 2100:
    #         print("period, curve " + str(i) + " " + str(best))
    #         bad.append(i)
    #     if best[4] > 2 and best[3] < 4:
    #         print("phi, curve " + str(i) + " " + str(best))
    #         bad.append(i)
            
    # print(bad)  #1 2, 19 3, 24 3, 31 3, 37 2, 50 2, 51 2, 54 2, 59 2, 65 2, 75 2, 83 2, 94, 2
            
        # i_sigmas.append(i_sigma)
        # i_taus.append(i_tau)
        # o_sigmas.append(o_sigma)
        # o_taus.append(o_tau)
    #    s_mcmc, t_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
        
    #    print [s_mcmc, t_mcmc]
        
    #    s_upper.append(s_mcmc[1])
    #    s_lower.append(s_mcmc[2])
    #    t_upper.append(t_mcmc[1])
    #    t_lower.append(t_mcmc[2])
    #    
    #    i_sigmas.append(i_sigma)
    #    o_sigmas.append(s_mcmc[0])
    #    
    #    i_taus.append(i_tau)
    #    o_taus.append(t_mcmc[0])
        
        
    #     n, sbins, patches = plt.hist(sigma_chains, 50)
    #     plt.show()
        
    # #    elem = np.argmax(n)
    # #    mcmc_sigma = sbins[elem]
    # #    mcmc_sigmas.append(mcmc_sigma)
        
    #     m, tbins, patches2 = plt.hist(tau_chains, 50)
    #     plt.show()
        
        
    # #    elem2 = np.argmax(m)
    # #    mcmc_tau = tbins[elem2]
    # #    mcmc_taus.append(mcmc_tau)
        
    #     o, abins, patches3 = plt.hist(amp_chains, 50)
    #     plt.show()
        
    #     p, pbins, patches4 = plt.hist(p_chains, 50)
    #     plt.show()
        
    #     q, phibins, patches5 = plt.hist(phi_chains, 50)
    #     plt.show()
        
    #     # H = plt.hist2d(sigma_chains, tau_chains, 50)
    #     # plt.show()
        
    #     H = np.histogramdd((sigma_chains, tau_chains, amp_chains, p_chains, phi_chains),50)
    #     # print(H)
    #     bincounts = H[0] # 5D array with the 50x50x50x50x50 grid
        
    #     maxb = np.amax(bincounts) # max value in the grid
        
        # print(maxb)
        
        # result = np.where(bincounts==np.amax(maxb))
        
        # print("The max bin is bin number " + str(result))
        # print(sbins[result[0]],tbins[result[1]],abins[result[2]],pbins[result[3]],phibins[result[4]])
    #    s = np.digitize(sigma_chains, sbins)
    #    t = np.digitize(tau_chains, tbins)
    #    scount = np.bincount(s)
    #    tcount = np.bincount(t) 
    #    totals = np.add(scount,tcount)
    #    maxindex = np.argmax(totals)
    #    print (maxindex)
    #     mcmc_sigmas.append(sbins[result[0]])
    #     mcmc_taus.append(tbins[result[1]])
    

    
    # plt.hist(mcmc_sigmas,bins=50)
    # plt.show()
    
    # plt.hist(mcmc_taus,bins=50)
    # plt.show()
    
    # plt.hist(mcmc_amps,bins=50)
    # plt.show()
    
    # plt.hist(mcmc_ps,bins=50)
    # plt.show()
    
    # # plt.hist(mcmc_phis,bins=50)
    # plt.hist(mcmc_t0s,bins=50)
    # plt.show()
    
        
    # fig, axs = plt.subplots(1,2)
    # axs[0].errorbar(i_sigmas, o_sigmas, fmt = 'o', c='g', label='S_in vs. Likelihood S_out')
    # axs[0].errorbar(i_sigmas, mcmc_sigmas, fmt = 'o', c='r', label='S_in vs. MCMC S_out')
    # lims = [
    #     np.min([axs[0].get_xlim(), axs[0].get_ylim()]),
    #     np.max([axs[0].get_xlim(), axs[0].get_ylim()]),
    # ]
    
    # axs[0].plot(lims, lims, '--', alpha=0.75, zorder=0)
    
    # axs[0].set_xlim(lims)
    # axs[0].set_ylim(lims)
    
    # axs[0].set_xlabel('Input Sigma')
    # axs[0].set_ylabel('Output Sigma')
    # axs[0].legend(loc=4,prop={'size': 7})
    # axs[0].set_title('S_in vs. S_out (L & MCMC)')
    
    # axs[1].errorbar(i_taus, o_taus, fmt = 'o', c='g', label='T_in vs. Likelihood T_out')
    # axs[1].errorbar(i_taus, mcmc_taus, fmt = 'o', c='r', label='T_in vs. MCMC T_out')
    # lims = [
    #     np.min([axs[1].get_xlim(), axs[1].get_ylim()]),
    #     np.max([axs[1].get_xlim(), axs[1].get_ylim()]),
    # ]
    
    # axs[1].plot(lims, lims, '--', alpha=0.75, zorder=0)
    
    # axs[1].set_xlim(lims)
    # axs[1].set_ylim(lims)
    
    # axs[1].set_xlabel('Input Tau')
    # axs[1].set_ylabel('Output Tau')
    # axs[1].legend(loc=4,prop={'size': 7})
    # axs[1].set_title('T_in vs. T_out (L and MCMC)')
    # plt.show()
    
    # plt.errorbar(i_sigmas, i_taus, fmt = 'o', c='black', label='Input Sigmas vs. Input Taus')
    # plt.errorbar(o_sigmas, o_taus, fmt = 'o', c='g', label='Likelihood Output Sigmas vs. Taus')
    # plt.errorbar(mcmc_sigmas, mcmc_taus, fmt = 'o', c='r', label='MCMC Output Sigmas vs. Taus')
    # plt.xlabel('Sigmas')
    # plt.ylabel('Taus')
    # plt.legend()
    # plt.title('All Sigmas vs. All Taus (Input, Likelihood and MCMC)')
    # plt.show()
    
    # chains = f['/'+str(31)+'/chains']
    # chains=chains.value
    
    # for b in np.arange(0,100):
    #     plt.plot(chains[b,:,0])
    #     plt.title("DRW L+C+A Sigma Chains")
    
    # # plt.savefig("schains" + str(i) + ".png")
    # plt.show()
    
    # for c in np.arange(0,100):
    #     plt.plot(chains[c,:,1])
    #     plt.title("DRW L+C+A Tau Chains")
    
    # # plt.savefig("tchains" + str(i) + ".png")
    # plt.show()
    
    # for d in np.arange(0,100):
    #     plt.plot(chains[d,:,2])
    #     plt.title("PG1302 DRWSIN L+C+A Amplitude Chains")
    
    # # plt.savefig("schains" + str(i) + ".png")
    # plt.show()
    
    # for e in np.arange(0,100):
    #     plt.plot(chains[e,:,3])
    #     plt.title("PG1302 DRWSIN L+C+A Period Chains")
    
    # # plt.savefig("tchains" + str(i) + ".png")
    # plt.show()
    
    # for f in np.arange(0,100):
    #     plt.plot(chains[f,:,4])
    #     plt.title("PG1302 DRWSIN L+C+A t_0 Chains")
    
    # # # plt.savefig("schains" + str(i) + ".png")
    # plt.show()
    
    
    return probs_di, best_di

def main():
    # from ChiSquareMinimization import drw_likelihood,Likelihood_drw_sinusoid
    
    # liu_L = Likelihood_drw_sinusoid()
    # zhu_L = 
    
    # for i in range(1,10):
    #     print(readMCMC('C:/Users/noahk/pg1302-research/cluster-data/L+C+A_mcmc_rand_sig_tau_no_outlier_info_' + str(i) + '.h5',2))
    
    # print(readMCMC('C:/Users/noahk/pg1302-research/cluster-data/final_curve_pg1302_LC_drw_info_4.h5',2))
    print(readMCMC('C:/Users/noahk/pg1302-research/cluster-data/final_curve_pg1302LCA_drw_no_outlier_info_2.h5',2))
    # print(readMCMC('C:/Users/noahk/pg1302-research/cluster-data/final_curve_pg1302LCA_drw_no_outlier_info_2.h5',2))
    # print(readMCMC('C:/Users/noahk/pg1302-research/cluster-data/final_curve10_LC_drwsin_no_outlier_info.h5',5))
    # from matplotlib import pyplot as plt
    # s = []
    # t = []
    # a = []
    # p = []
    # t0 = []
    
    # for i in range(1,99):
    #     filename =r'C:\Users\noahk\pg1302-research\cluster-data\drwsin_mcmc' + str(i) + 'info.h5'
    #     probs_di, best_di = readMCMC(filename)
    #     s.append(best_di[i][0])
    #     t.append(best_di[i][1])
    #     a.append(best_di[i][2])
    #     p.append(best_di[i][3])
    #     t0.append(best_di[i][4])
    #     print(probs_di, best_di)
    
    # # for i in range(4, 76):
    # #     filename =r'C:\Users\noahk\pg1302-research\cluster-data\drwsin_mcmc' + str(i) + 'info.h5'
    # #     probs_di, best_di = readMCMC(filename,5)
    # #     s.append(best_di[i][0])
    # #     t.append(best_di[i][1])
    # #     a.append(best_di[i][2])
    # #     p.append(best_di[i][3])
    # #     t0.append(best_di[i][4])
    # #     print(probs_di, best_di)
        
    # # for i in range(77, 99):
    # #     filename =r'C:\Users\noahk\pg1302-research\cluster-data\drwsin_mcmc' + str(i) + 'info.h5'
    # #     probs_di, best_di = readMCMC(filename,5)
    # #     s.append(best_di[i][0])
    # #     t.append(best_di[i][1])
    # #     a.append(best_di[i][2])
    # #     p.append(best_di[i][3])
    # #     t0.append(best_di[i][4])
    # #     print(probs_di, best_di)
        
    # plt.hist(s,50)
    # plt.title("Best Sigma Values, True = .11")
    # plt.show()
    
    # plt.hist(t,50)
    # plt.title("Best Tau Values, True = 100")
    # plt.show()
    
    # plt.hist(a,50)
    # plt.title("Best Amplitude Values, True = .1")
    # plt.show()
    
    # plt.hist(p,50)
    # plt.title("Best Period Values, True = 1825")
    # plt.show()
    
    # plt.hist(t0,50)
    # plt.title("Best t0 Values, True = 100")
    # plt.show()
    
    return None

if __name__ == "__main__":
    main()
    
    
