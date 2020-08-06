def PeriodogramRbandDRW(QSOId,k):
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
#    just insert a file containing the real data you want to run PLS ond DRW on
    data=np.loadtxt(r'C:\Users\noahk\pg1302-research\observation-data\pg1302_data_LINEAR+CRTS+ASASSN.txt')
    data2=np.loadtxt(r'C:\Users\noahk\pg1302-research\observation-data\specific_pg1302_data.txt')
#    print data
    MJD=data[:,0]
    mag=data[:,1]
    magErr=data[:,2]

    SortedInd = np.argsort(MJD)
    MJD = MJD[SortedInd]
    mag = mag[SortedInd]
    magErr = magErr[SortedInd]
    
    MJD2=data2[:,0]
    mag2=data2[:,1]
    magErr2=data2[:,2]

    SortedInd = np.argsort(MJD2)
    MJD2 = MJD2[SortedInd]
    mag2 = mag2[SortedInd]
    magErr2 = magErr2[SortedInd]
    
    # print(np.where(MJD==57207.806))
    # print(mag[286])
    
    print(magErr[286])
    print(magErr)
    MJD3 = np.delete(MJD,286)
    mag3 = np.delete(mag,286)
   
    magErr3 = np.delete(magErr,286)
    
 
#    #find frequency grid
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
    BestFitTau=100
    BestFitSigma=.11
    BestFitMean=14
    
    
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
    x=SimulateTimeSeriesDRW(BestFitTau,BestFitSigma,NumberOfPoints,ind,MJD,magErr,BestFitMean,N_bootstraps,deltat)
#    print x
    x = mag
    plt.errorbar(MJD, x, yerr=magErr, fmt='o', ecolor='black')
    plt.show()

#    MaxPsim=np.max(Psim,axis=0)
#    AvgPsim=np.mean(Psim,axis=0)
#
#    ################################################################################################
#    #periodogram
    from astroML.time_series import lomb_scargle
    P_LS = lomb_scargle(MJD, x, magErr, omegaSlope, generalized=True)
    P_LS2 = lomb_scargle(MJD2, mag2, magErr2, omegaSlope, generalized=True)
    P_LS3 = lomb_scargle(MJD3, mag3, magErr3, omegaSlope, generalized=True)
    P_LS = np.log(P_LS)
    P_LS2 = np.log(P_LS2)
    P_LS3 = np.log(P_LS3)
    ftemp = omegaSlope / (2*np.pi*86400)
    logf = np.log10(ftemp)
#     max_peak = np.amax(P_LS)
#     print(max_peak)
#     a = np.argmax(P_LS)
#     max_f = logf[a]
# #    print(max_peak)
# #    print max_peak
#     max_result = np.where(P_LS == max_peak)
#     fmax = ftemp[max_result]
#     print(fmax)
#     print(max_f)
# #    print fmax
#     half_max = max_peak/2
#     nearest = (np.abs(P_LS[0:800] - half_max)).argmin()
#     h1 = ftemp[nearest]
#     h2 = ftemp[639]
# #    print h1
# #    print h2
#     fwhm = h2-h1
#    print fwhm
    plt.plot(logf, P_LS, '--', lw=1, zorder=1, color="black", label='L+C+A w/ outlier')
    plt.plot(logf, P_LS2, '--', lw=1, zorder=1,color = 'blue', label='L+C')
    plt.plot(logf, P_LS3, '--', lw=1, zorder=1, color='red', label='L+C+A w/o outlier')
    plt.xlabel('Log Frequency (Hz)')
    plt.ylabel('P_LS')
#    plt.xlim(-8.8,-6.9)
    plt.xlim(-9,-7.2)
    plt.legend()
    # plt.xlim(0,.000000015)
    # plt.axvline(h1)
    # plt.axvline(h2)
    plt.title("LS Periodogram of PG1302 L+C+A Data")
    plt.show()
    
#    Psim2,x=SimulateTimeSeriesDRW(BestFitTau,BestFitSigma,NumberOfPoints,ind,MJD,magErr,omegaSlope,BestFitMean,N_bootstraps,deltat)
#
#    from PvalueOfPeak import BinnedFrequencies
#    NumberOfBins=50
#    PvalueBin=BinnedFrequencies(P_LS,Psim,NumberOfBins)
#
#    MinPval=np.min(PvalueBin)
#    IndMinPval=np.where(PvalueBin==MinPval)

    
    ###################################################################################################
    #Save data
    
    
#    import h5py
#    Filename='../DRW1000Iterations1/QSO'+str(QSOId)+'_1000Iterations.h5'
#    f = h5py.File(Filename, "a")
#    #    f[str(QSOId)+'/BestFitSigma']=BestFitSigma
#    #    f[str(QSOId)+'/BestFitTau']=BestFitTau
#    #    f[str(QSOId)+'/BestFitMean']=BestFitMean
#    #    f[str(QSOId)+'/ChiSquare']=ChiSq
#
#    f[str(QSOId)+'/MinPvalBin']=MinPval
#    f[str(QSOId)+'/IndMinPval']=IndMinPval
#
#    #    f[str(QSOId)+'/NumberOfPeaks']=NumberOfPeaks
#    f[str(QSOId)+'/PvalueBin']=PvalueBin
#
#    f[str(QSOId)+'/PLS']=P_LS
#    #    f[str(QSOId)+'/omegaSlope']=omegaSlope
#    f[str(QSOId)+'/MaxPsim']=MaxPsim
#    f[str(QSOId)+'/AvgPsim']=AvgPsim
#
#    #    f[str(QSOId)+'/MJD']=MJD
#    f[str(QSOId)+'/mag']=x
#    #    f[str(QSOId)+'/magErr']=magErr
#    f.close()
##k=1

PeriodogramRbandDRW(1,1)