def SimulateTimeSeriesDRW(tau,sigma,NumberOfPoints,ind,MJD,magErr,Mean,N_bootstraps,deltat):
    import numpy as np
    from astroML.time_series import lomb_scargle
    #    from LS_periodogram import lomb_scargle
    from generate_power import generate_power_DRW
    from matplotlib import pyplot as plt
    from scipy.fftpack import fft
    from astroML.fourier import PSD_continuous, FT_continuous
    
    dt=deltat
    Nlength=1
    Npoints=Nlength*NumberOfPoints
    #    print 'sigma**2*tau/2'
    #    print sigma**2*tau/2
    sigma=sigma/np.sqrt(dt)
    #    tau=tau/np.sqrt(dt)
#    Psim=np.zeros(len(omega))
    #    data=np.genfromtxt("../ChainsDat/chain_215131.dat.myrun", unpack=True, skiprows = 1)
    ##    sigma=data[0,:]
    #    print sigma
    #    tau=data[1,:]
#    print tau
    
    #    from DRW_Javelin import generateTrueLC
    v=np.zeros(N_bootstraps)
    for i in range(1,N_bootstraps):
#        print i
        
        #        time,xF=generateTrueLC(sigma,tau,covfunc="drw")
        x = generate_power_DRW(Npoints, dt, sigma, tau, generate_complex=False, random_state=None)
#        print x
        x=x*np.sqrt(Npoints)
        
        #        v[i]=np.var(x)
        #        time=np.arange(Npoints)*deltat
        
        #        plt.scatter(time,x)
        #        plt.show()
        
        #        f, PSD = PSD_continuous(time, x)
        
        #        P=4*sigma**2*tau/(1+4*np.pi**2*f**2*tau**2)
        #        plt.loglog(f,PSD)
        #        plt.loglog(f,P,'r')
        #        plt.show()
        
        xtemp=x[ind-1]
        #        xtemp=xF
#        a = np.random.normal(0,1,len(magErr))
        a = np.random.normal(0,1,None)
        PhotErr=a*np.mean(magErr)
        xF=xtemp
        xF=xtemp+PhotErr
        xF=xF+Mean-np.mean(xF)

#p=lomb_scargle(MJD, xF, magErr, omega, generalized='True')
#
#Psim=np.vstack((Psim,p))
    #    print 'average'
    #    print np.mean(v)
    return xF

