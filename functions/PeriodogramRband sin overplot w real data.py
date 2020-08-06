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
    data=np.loadtxt('C:/Python27/Data/pg1302_data_LINEAR+CRTS+ASASSN.txt')
#    print data
    MJD=data[:,0]
#    mag=data[:,1]
    P = 2000
    mag =np.sin(2 * np.pi * MJD / P)
    magErr=data[:,2]/100


    SortedInd = np.argsort(MJD)
    MJD = MJD[SortedInd]
    mag = mag[SortedInd]
    magErr = magErr[SortedInd]
    
    data2 = np.loadtxt('C:/Python27/Data/pg1302_data_LINEAR+CRTS+ASASSN.txt')
    
    MJD2=data2[:,0]
    mag2=data2[:,1]
    magErr2=data2[:,2]

    SortedInd = np.argsort(MJD2)
    MJD2 = MJD2[SortedInd]
    mag2 = mag2[SortedInd]
    magErr2 = magErr2[SortedInd]
 
#    #find frequency grid
    from OptimalFrequencyGrid import OptimalFrequencyGrid
    omega, omegaSlope = OptimalFrequencyGrid(MJD)
    omega2, omegaSlope2 = OptimalFrequencyGrid(MJD2)

    
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
    x=SimulateTimeSeriesDRW(BestFitTau,BestFitSigma,NumberOfPoints,ind,MJD,mag,magErr,BestFitMean,N_bootstraps,deltat)
#    print x
    x = mag
    plt.errorbar(MJD, x, yerr=magErr, fmt='o', ecolor='black')
    plt.show()
    
    y = mag2
    plt.errorbar(MJD2, y, yerr=magErr2, fmt='o', ecolor='black')
    plt.show()

#    MaxPsim=np.max(Psim,axis=0)
#    AvgPsim=np.mean(Psim,axis=0)
#
#    ################################################################################################
#    #periodogram
    from astroML.time_series import lomb_scargle
    LIN_CRTS_P_LS = lomb_scargle(MJD, x, magErr, omegaSlope, generalized=True)
    entire_curve_P_LS = lomb_scargle(MJD2, y, magErr2, omegaSlope2, generalized=True)
    
    plot1, = plt.plot(omegaSlope, LIN_CRTS_P_LS, '-', c='blue', lw=1, zorder=1)
    plot2, = plt.plot(omegaSlope2, entire_curve_P_LS, '-', c='red', lw=1, zorder=1)
    plt.xlabel('Frequency (2pi/days)')
    plt.ylabel('P_LS')
    plt.title("LS Periodograms of LINEAR+CRTS+ASASSN Data Overplotted with Corresponding Sinusoid")
    plt.legend([plot1,plot2],["Sinusoid w/ P=2000", "LINEAR+CRTS+ASASSN"])

    
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