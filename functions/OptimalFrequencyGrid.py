def OptimalFrequencyGrid(MJD):
    import numpy as np
    
#    from the baseline timescale
    Tmax=np.max(MJD)
    omegaMin=2*np.pi/(Tmax)

    MeanDeltaT=np.median(np.diff(MJD))
    Hour=100
    if MeanDeltaT<Hour:
        MeanDeltaT=Hour
        OneHourFlag=1
    
#    pseudo-Nyquist frequency
    omegaMax=np.pi/MeanDeltaT

    N2=1000

    omega=np.linspace(omegaMin,omegaMax,N2)
    omegaSlope=np.logspace(np.log10(omegaMin),np.log10(omegaMax),N2)
    return omega, omegaSlope
    
    
    
