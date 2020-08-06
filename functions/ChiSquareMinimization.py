
def ChiSquareMinimization(mag,magErr,MJD):
    import numpy as np
    from scipy.optimize import fmin,leastsq
    import matplotlib.pyplot as plt
    
    
#    period=np.linspace(1,np.max(MJD),10,endpoint=True)
    period=np.logspace(np.log10(1),np.log10(np.max(MJD)),20,endpoint=True)

    print (period)
#    AmplitudeMin=0.01
#    AmplitudeMax=10
#    AmplitudeMin=0.01/np.sqrt(360)
#    AmplitudeMax=10/np.sqrt(360)
    
    AmplitudeMin=0.01
    AmplitudeMax=2
    
#    AmplitudeMin=0.02
#    AmplitudeMax=0.08
    Amplitude=np.logspace(np.log10(AmplitudeMin),np.log10(AmplitudeMax),20)
#    Amplitude=np.linspace(AmplitudeMin,AmplitudeMax,10)

    print (Amplitude)
    LogL=[]
    chi=[]
    meanLC=[]
    Tau=[]
    sigma=[]
    for j in np.arange(len(period)):
        for k in np.arange(len(Amplitude)):
#  
#            p0=[Amplitude[k],period[j]]

            LnL,chiS,qhat=CovarianceMatrices(Amplitude[k],period[j],mag,magErr,MJD)

            chi.append(chiS[0][0])
            LogL.append(LnL[0][0])
            meanLC.append(qhat[0][0])
            Tau.append(period[j])
            sigma.append(Amplitude[k])

#    print LogL
    LogLL=np.reshape(LogL,(len(period),len(Amplitude)))

    LogL=np.array(LogL)
    chi=np.array(chi)
    meanLC=np.array(meanLC)
    Tau=np.array(Tau)
    sigma=np.array(sigma)
                
    SortedInd = np.argsort(LogL)
    LogL=LogL[SortedInd]            
    chi=chi[SortedInd]
    
    meanLC=meanLC[SortedInd]
    Tau=Tau[SortedInd]
    sigma=sigma[SortedInd]

    chi[:] = chi[::-1]
    LogL[:]=LogL[::-1]
    meanLC[:]=meanLC[::-1]
    Tau[:]=Tau[::-1]
    sigma[:]=sigma[::-1]
#                
    print (chi[0:10])
    print (LogL[0:10])
#    print meanLC[0:50]
    print (Tau[0:10])
    print (sigma[0:10])
#    print np.mean(Tau[0:10])
#    print np.mean(sigma[0:10])
#    Tau2=[]
#    Sigma2=[]
#    Mu2=[]
#    Ln2=[]
#    chi2=[]
#    pp=np.append(0.1,period)
#    pp=np.append(pp,10*np.max(MJD))
#                
#    vv=np.append(0.1*np.min(Amplitude),Amplitude)
#    vv=np.append(vv,10*np.max(Amplitude))
##    print pp
##
#    for l in np.arange(len(Tau[0:9])):
#        ind1=np.where(pp==Tau[l])[0]
#        pp2=np.logspace(np.log10(pp[ind1-1]),np.log10(pp[ind1+1]),10)
##        print pp2
#        ind2=np.where(vv==sigma[l])[0]
#        vv2=np.logspace(np.log10(vv[ind2-1]),np.log10(vv[ind2+1]),10)
##        print vv2
##        klsdkl
#                
#        for jj in np.arange(len(pp2)):
#            for kk in np.arange(len(vv2)):
#                LnL,chiS,qhat=CovarianceMatrices(vv2[kk],pp2[jj],mag,magErr,MJD)
##
#                chi2.append(chiS[0][0])
#                Ln2.append(LnL[0][0])
#                Mu2.append(qhat[0][0])
#                Tau2.append(pp2[jj])
#                Sigma2.append(vv2[kk])
##
##                #    print LogL
###                LogLL=np.reshape(LogL,(len(period),len(Amplitude)))
##                
#    Ln2=np.array(Ln2)
#    chi2=np.array(chi2)
#    Mu2=np.array(Mu2)
#    Tau2=np.array(Tau2)
#    Sigma2=np.array(Sigma2)
##
#    SortedInd = np.argsort(Ln2)
#    Ln2=Ln2[SortedInd]
#    chi2=chi2[SortedInd]
#                
#    Mu2=Mu2[SortedInd]
#    Tau2=Tau2[SortedInd]
#    Sigma2=Sigma2[SortedInd]
##
#    chi2[:] = chi2[::-1]
#    Ln2[:]=Ln2[::-1]
#    Mu2[:]=Mu2[::-1]
#    Tau2[:]=Tau2[::-1]
#    Sigma2[:]=Sigma2[::-1]
##
#                
#    print chi2[0:10]
#    print Ln2[0:10]
#    print Tau2[0:10]
#    print Sigma2[0:10]
#

#
#                
#                
#                
#    kljaskld
                
    
                
#        print chi
#        print LogL
#    chii=np.reshape(chi,(len(period),len(Amplitude)))
##    print chi
#    print LogL
#    print chi
#    meanLC=np.reshape(meanLC,(len(period),len(Amplitude)))
##
#    for k in np.arange(len(period)):
#                
#        plt.plot(Amplitude,LogLL[k,:])
##    plt.show()
#    from matplotlib import cm
#    import mpl_toolkits.mplot3d.axes3d as axes3d
#    
#    index=np.where(LogL>0)
#                
#    fig = plt.figure()
#    ax = fig.add_subplot(111, projection='3d')
#    ax.plot(Tau, sigma, LogL, linestyle="none", marker="o", mfc="none", markeredgecolor="red")
#    ax.plot_surface(Tau[index], sigma[index], LogL[index])
##    plt.contour(Tau,sigma,LogL)
#    plt.show()

#    ind=np.where(LogL==np.max(LogL))
#    BestFitTau = period[ind[0][0]]
#    BestFitSigma = Amplitude[ind[1][0]]
#    BestFitMean = meanLC[ind]
#    ChiSqDOF = chi[ind]
#    return BestFitSigma, BestFitTau, BestFitMean, ChiSqDOF

    return sigma[0], Tau[0], meanLC[0],chi[0]



import scipy.linalg as sl
import numpy as np

def drw_likelihood(sigma,log_tau,mag,magErr,MJD):
    if sigma < 0.00001 or sigma > 100:
        return -np.inf
    # elif tau < 1 or tau > 1000:
    elif log_tau < 0 or log_tau > 3:
        return -np.inf
    else:
#from kozlowski 2010
        tx, ty = np.meshgrid(MJD, MJD)
    
        S = np.exp(- np.abs(tx-ty) / 10 ** log_tau)
    
        S *= sigma**2 #eq.(1)
        
        N=np.diag(magErr**2,k=0)
        
        C=S+N #defined below (A4)
        
    
        cf = sl.cho_factor(C)
        CInv = sl.cho_solve(cf, np.eye(C.shape[0]))
        LogDetC = np.sum(2*np.log(np.diag(cf[0]))) #first determinant in A8
    
    
        L = np.ones(len(MJD))
        CL = np.dot(CInv, L)
        Cq =  1.0 / np.dot(L.T, CL) #eq. (A7)
    
        detCqInv = np.abs(1.0 / Cq) #second determinant in A8
    
        Cf = CInv - np.outer(CL, CL.T) * Cq #eq. (A7)
    
        chiSq=np.transpose(mag).dot(Cf).dot(mag) #exp eq. (A8)
    
        LogLikelihood=-0.5*LogDetC-0.5*np.log(detCqInv)-chiSq/2 #eq. (A8)
        
    return np.real(LogLikelihood)




# def Likelihood_drw_sinusoid(sigma,tau,A,Period,phase,mag,magErr,MJD):
def Likelihood_drw_sinusoid(sigma,log_tau,A,Period,t0,mag,magErr,MJD):
#def Likelihood_drw_sinusoid(sigma,tau,Period,mag,magErr,MJD):
#drw likelihood from kozlowski 2010

    # Sine_Wave=A*np.sin(2*np.pi/Period*MJD+phase)
    Sine_Wave=A*np.sin((2*np.pi/Period)*(MJD - t0))
#    Sine_Wave=np.sin(2*np.pi/Period*MJD)

    y=mag-Sine_Wave         #this changes if we change generation


    tx, ty = np.meshgrid(MJD, MJD)
    S = np.exp(- np.abs(tx-ty) / 10 ** log_tau)
    
    S *= sigma**2 #eq.(1)
    
    N=np.diag(magErr**2,k=0)
    
    C=S+N #defined below (A4)
    
    
    cf = sl.cho_factor(C)
    CInv = sl.cho_solve(cf, np.eye(C.shape[0]))
    LogDetC = np.sum(2*np.log(np.diag(cf[0]))) #first determinant in A8
    
    
    L = np.ones(len(MJD))
    CL = np.dot(CInv, L)
    Cq =  1.0 / np.dot(L.T, CL) #eq. (A7)
    
    detCqInv = np.abs(1.0 / Cq) #second determinant in A8
    
    Cf = CInv - np.outer(CL, CL.T) * Cq #eq. (A7)
    
    chiSq=np.transpose(y).dot(Cf).dot(y) #exp eq. (A8)
    
    LogLikelihood=-0.5*LogDetC-0.5*np.log(detCqInv)-chiSq/2 #eq. (A8)  #check second minus sign bc already in 251
    
    return np.real(LogLikelihood)
