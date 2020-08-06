# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 14:36:13 2019

@author: noahk
"""
import sys

def Simulate(count,filename):
    import matplotlib.pyplot as plt
    #import h5py
    import numpy as np
#    #########################################################################################################
#    #import light curve
#    from QSOLightCurveH5FinalCatalog import QSOLightCurvesFinalCatalog
#    MJD, mag, magErr = QSOLightCurvesFinalCatalog(QSOId)
#
#
#    #sort the observations with time
#
#    SortedInd = np.argsort(MJD)
#    MJD = MJD[SortedInd]
#    mag = mag[SortedInd]
#    magErr = magErr[SortedInd]
#
#    ##########################################################################################################
#    #weighted by the photometric average observations in one night bins
#    from BinLightCurve import BinLightCurve
#    MJD, mag, magErr = BinLightCurve(MJD, mag, magErr)
#
#    MJD = MJD-np.min(MJD)
#
#####################################################
#Replace with PG1302 light curve
######################################################
#    import os
#    if os.path.exists("Realizations_Significant_Peaks.h5"):
#        os.remove("Realizations_Significant_Peaks.h5")
#    else:
#        print("The file does not exist")
    
#    Describing the desired sinusoid
     
    # P = 1833
    # A = .131
    # t0 = 21.2

    num = 0 #counter for number of significant curves
    
    for i in range(count): # how many curves do you want to generate
#        BestFitTau=np.random.uniform(50, 500)
#        BestFitSigma=np.random.uniform(0.05, 0.3)
#        temp = np.random.uniform(-3,1)
#        rho_i = 10**temp
#        t_exp = 5780 #length of linear+crts+asassn
#        BestFitTau= rho_i * t_exp
#        BestFitTau = 100
#        BestFitTau = np.random.uniform(20, 500)
#       BestFitSigma = .11
#        BestFitSigma= np.random.uniform(0.02, 0.5)
        BestFitMean=0
        data=np.loadtxt(r"C:\Users\noahk\pg1302-research\observation-data\pg1302_data_LINEAR+CRTS+ASASSN.txt")
    #    print data
        MJD=data[:,0]
#        MJD=np.linspace(0, t_exp, num=445)
        
       
        # BestFitTau = 100
        # BestFitSigma = .11
        
        
        
        # mag=data[:,1]
        magErr=data[:,2]
#        sdss = .000172
#        ogle = .000617
#        magErr=np.random.normal(0,ogle,445)
        flag = data[:,3]
        SortedInd = np.argsort(MJD)
        MJD = MJD[SortedInd]
    #    mag = mag[SortedInd]
        magErr = magErr[SortedInd]
        flag = flag[SortedInd]
        
        #DELETE OUTLIER
        
        MJD = np.delete(MJD,286)
        magErr = np.delete(magErr,286)
        flag = np.delete(flag,286)
        
        # find frequency grid
        from OptimalFrequencyGrid import OptimalFrequencyGrid
        omega, omegaSlope = OptimalFrequencyGrid(MJD)
    
        
        #################################################################################################
        #indices to downsample the uniform grid
        deltat=1
        from MinimumTimeResolution import MinimumTimeResolution
        ind, NumberOfPoints = MinimumTimeResolution(MJD,deltat)
        #
        
    #    #################################################################################################
    #    #identify best fit parameters for DRW model
    #    import h5py
    #    filename='../QSO1000Iterations'+str(k)+'/QSO1'+str(QSOId)+'_1000Iterations.h5'
    #    f=h5py.File(filename,"r")
    #    BestFitTau=f['/'+str(QSOId)+'/BestFitTau']
    #    BestFitTau=BestFitTau.value
    #    BestFitSigma=f['/'+str(QSOId)+'/BestFitSigma']
    #    BestFitSigma=BestFitSigma.value
    #
    #    BestFitMean=f['/'+str(QSOId)+'/BestFitMean']
    #    BestFitMean=BestFitMean.value
    #    f.close()
    #
    
    #############################################
    #Fix sigma and tau from Charisi+2015
    #############################################
        
        BestFitTau=np.random.uniform(1,1000)
        BestFitSigma=np.random.uniform(.001,.5)

        
        
        #    from ChiSquareMinimization import ChiSquareMinimization
        #    BestFitSigma, BestFitTau, BestFitMean, ChiSq=ChiSquareMinimization(mag, magErr, MJD)
        #    print BestFitSigma
        #    print BestFitTau
        #    print BestFitMean
        #    print ChiSq
        
        ###############################################################################################
        # DRW Simulations
        N_bootstraps=10
        from SimulateTimeSeries import SimulateTimeSeriesDRW
        x = SimulateTimeSeriesDRW(BestFitTau,BestFitSigma,NumberOfPoints,ind,MJD,magErr,BestFitMean,N_bootstraps,deltat)
        ind=np.where(flag==0)
        y=x[ind]
        MJDy=MJD[ind]
        magErry=magErr[ind]
        
        #NO OUTLIER 
        

        
#        Multiplying Lightcurves
        
#        temp = []
#        i = 0
#        for m in x:
#            t = MJD[i]  #get the day associated with that magnitude
#            sinFactor = A * np.sin((2 * np.pi * t / P) + phi)  # create factor
#            v = m - 14  #find difference between magnitude and hard-coded mean
#            vsin = v * sinFactor  # multiply the two
#            temp.append(14 + vsin) # recreate array with new variance
#            i+=1
#            
#        drwsin = np.asarray(temp)
        
        # sinusoid = A * np.sin((2 * np.pi * MJD / P) + phi)
        
        # sinusoid = A * np.sin(((2 * np.pi / P)*(MJD - t0)))
        
        # sinusoid = A * np.sin(((2 * np.pi / P)*(MJDy - t0)))
        
        # drwsin = x + sinusoid
        
        # drwsin = y + sinusoid
        
#        print("Sigma: " + str(BestFitSigma))
#        print("Tau: " + str(BestFitTau))
#        plt.errorbar(MJD, drwsin, yerr=magErr, fmt='o', ecolor='black')
        plt.errorbar(MJD, x, fmt='o', ecolor='black')
        plt.xlabel('MJD (Days)')
        plt.ylabel('Magnitude')
        plt.title("Pure DRW Curve")
#       plt.legend()
        plt.show()
        
#         plt.errorbar(MJD, sinusoid, fmt='o', c='red', ecolor='black')
#         plt.xlabel('MJD (Days)')
#         plt.ylabel('Magnitude')
#         plt.title("Pure Sinusoidal Curve")
# #       plt.legend()
#         plt.show()
        
#         plt.errorbar(MJD, drwsin, fmt='o', c='purple', ecolor='black')
#         plt.xlabel('MJD (Days)')
#         plt.ylabel('Magnitude')
#         plt.title("DRW + Sinusoid Curve")
# #       plt.legend()
#         plt.show()
#        
        # plt.errorbar(MJDy, y, yerr=magErry, fmt='o', ecolor='black')
        # plt.xlabel('MJD (Days)')
        # plt.ylabel('Magnitude')
        # plt.title("Simulated Light Curve w/ Only LINEAR+CRTS Contributions")
        # # plt.legend()
        # plt.show()
    
    #    MaxPsim=np.max(Psim,axis=0)
    #    AvgPsim=np.mean(Psim,axis=0)
    #
    #    ################################################################################################
    #    #periodogram
#         from astroML.time_series import lomb_scargle
#         P_LS = lomb_scargle(MJDy, y, magErry, omegaSlope, generalized=True)
# #        P_LS = lomb_scargle(MJD, x, magErr, omegaSlope, generalized=True)
#         frequency = omegaSlope / 6.28 / 86400
#         logf = np.log10(frequency)
# ##        finding the highest peak in the DRW curve's periodogram
#         # plt.plot(logf, P_LS, '-', c='blue', lw=1, zorder=1)
#         # plt.xlabel('Log10 of Frequency (Hz)')
#         # plt.ylabel('P_LS')
#         # plt.xlim(-8.8,-6.9)
#         # plt.title("LS Periodogram of LINEAR+CRTS DRW Curve")
#         # # plt.legend([plot1],["P_LS of LINEAR+CRTS DRW Curve"])
#         # plt.show()
# #        
#         max_peak = np.amax(P_LS)
        # print(max_peak)
        
    #####################################################################################################    
        # if max_peak > 0.8300:
        #     num+=1
        #     #we need to save everything about the curves that are significant
        #     import h5py
        #     f = h5py.File(filename, "a")
        #     f[str(i)+'/MJD']=MJD
        #     f[str(i)+'/mag']=x
        #     f[str(i)+'/magErr']=magErr
        #     f[str(i)+'/omegaSlope']=omegaSlope
        #     f[str(i)+'/maxpeak']=max_peak
        #     f[str(i)+'/PLS']=P_LS
        #     f[str(i)+'/frequency']=frequency
        #     f[str(i)+'/logf']=logf
        #     f[str(i)+'/MJDy']=MJDy
        #     f[str(i)+'/y']=y
        #     f[str(i)+'/magErry']=magErry
        #     f[str(i)+'/sigma']=BestFitSigma
        #     f[str(i)+'/tau']=BestFitTau
        #     f.close()

        #     print("Success number " + str(num) + " occured at curve " + str(i))
        import h5py
        f = h5py.File(filename, "a")
        f[str(i)+'/MJD']=MJD
        f[str(i)+'/mag']=x
        # f[str(i)+'/mag']=drwsin
        # f[str(i)+'/P']=P
        # f[str(i)+'/A']=A
        # f[str(i)+'/phi']=phi
        # f[str(i)+'/t0']=t0
        f[str(i)+'/magErr']=magErr
        
#        for saving only the linear+crts contributions
        
        f[str(i)+'/MJDy']=MJDy
        f[str(i)+'/y']=y
        f[str(i)+'/magErry']=magErry
        
        f[str(i)+'/omegaSlope']=omegaSlope
        f[str(i)+'/sigma']=BestFitSigma
        f[str(i)+'/tau']=BestFitTau
        f.close()
        
    # print("Simulated " + str(i))
    # print("Found " + str(num) + " successful curves.")
    
    return x

if __name__ == "__main__":
    Simulate(sys.argv)
    
