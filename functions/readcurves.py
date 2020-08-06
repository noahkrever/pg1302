# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 16:58:37 2019

@author: noahk
"""
import h5py
import matplotlib.pyplot as plt
import numpy as np

peaks = []
frequencies = []

A_peaks = []
A_frequencies = []

test_peaks = []
test_frequencies = []

def readcurve(i,filename):
    # filename=r'C:\Users\noahk\pg1302-research\cluster-data'
#    filename='500_WS_PSD2_t200_m0_2.h5'
    f=h5py.File(filename,"r")
    MJD=f[str(i)+'/MJD']
    MJD=MJD.value
    mag= f[str(i)+'/mag']
    mag=mag.value
    magErr= f[str(i)+'/magErr']
    magErr=magErr.value
    omegaSlope = f[str(i)+'/omegaSlope']
    omegaSlope=omegaSlope.value
    # PLS = f[str(i)+'/PLS']
    # PLS=PLS.value
    # max_peak = f[str(i)+'/maxpeak']
    # max_peak=max_peak.value
    # MJDy=f[str(i)+'/MJDy']
    # MJDy=MJDy.value
    # y= f[str(i)+'/y']
    # y=y.value
    # magErry= f[str(i)+'/magErry']
    # magErry=magErry.value
    # logf= f[str(i)+'/logf']
    # logf=logf.value
    # frequency=f[str(i)+'/frequency']
    # frequency=frequency.value

    
    # data=np.loadtxt(filename)
    # days = data[:,0]
    # flag = data[:,3]
    
    # SortedInd = np.argsort(days)
    # days = days[SortedInd]
    # flag = flag[SortedInd]
    
    # ind=np.where(flag==0)
    # y=mag[ind]
    # MJDy=MJD[ind]
    # magErry=magErr[ind]
    
#    print mag
#    print y
#    print MJD
#    print mag
#    print magErr
#    print omegaSlope
    
    # from OptimalFrequencyGrid import OptimalFrequencyGrid
    
    # omega, omegaSlope = OptimalFrequencyGrid(MJD)

    # from astroML.time_series import lomb_scargle
    
#    P_LS = lomb_scargle(MJDy, y, magErry, omegaSlope, generalized=True)
#    ftemp = omegaSlope / (2*np.pi*86400)
#    logf = np.log10(ftemp) 
    
    # plt.errorbar(MJD, mag, yerr=magErr, fmt='o', ecolor='black')
    # plt.xlabel('MJD (Days)')
    # plt.ylabel('Magnitude')
    # plt.title("Simualated DRW+Sinusoid curve " + str(i))
    # plt.legend()
    # plt.show()

#    plt.errorbar(MJDy, y, yerr=magErry, fmt='o', ecolor='black')
#    plt.xlabel('MJD (Days)')
#    plt.ylabel('Magnitude')
#    plt.title("LINEAR+CRTS+ASASSN DRW Simulated Light Curve")
#    plt.legend()
#    plt.show()

#    plot1, = plt.plot(ftemp, P_LS, '-', c='blue', lw=1, zorder=1)
#    plt.xlabel('Frequency (Hz)')
#    plt.ylabel('P_LS')
#    plt.xlim(.00000001,.00001)
#    plt.xlim(-8.8,-6.9) - fits an appropriate x axis for log of f
#    plt.title("LS Periodogram of Successful DRW Curve")
#    plt.legend([plot1],["P_LS of Successful DRW Curve"])
#    plt.show()

#    max_peak = np.amax(P_LS)
#    
#    result = np.where(P_LS == max_peak)
#    frequency = ftemp[result]
#    log = np.log10(frequency)
#    
#    frequencies.append(log)
#    peaks.append(max_peak)
    
#    return (y,magErry,MJDy)
    # return y,magErry,MJDy,omegaSlope
    return mag,magErr,MJD,omegaSlope
#    print "The max peak of this curve is " + str(max_peak) + " and it occurs at a frequency of " + str(frequency) + " Hz."
#    f.close()

#filename='Realizations_Significant_Peaks.h5'
#filename='Rerun_Realizations_Significant_Peaks.h5'
def main():
    import h5py
    from readMCMC import readMCMC
    s_in_LC = []
    s_out_LC = []
    t_in_LC = []
    t_out_LC = []
    
    s_in_LCA = []
    s_out_LCA = []
    t_in_LCA = []
    t_out_LCA = []
    
    s_lower = []
    s_upper = []
    t_lower = []
    t_upper = []
    for i in range(1,10):
        # filename1 = r'C:\Users\noahk\pg1302-research\cluster-data\L+C_mcmc_rand_sig_tau_info_' + str(i) + '.h5'
        filename1 = r'C:\Users\noahk\pg1302-research\cluster-data\L+C+A_mcmc_rand_sig_tau_no_outlier_info_' + str(i) + '.h5'
        f=h5py.File(filename1,"r")
        
        filename2 = 'L+C+A_mcmc_rand_sig_tau_no_outlier_' + str(i) + '.h5'
        f2 = h5py.File(filename2,"r")
        
        # filename3 = r'C:\Users\noahk\pg1302-research\cluster-data\L+C+A_mcmc_rand_sig_tau_info_' + str(i) + '.h5'
        # f3 = h5py.File(filename3,"r")

        for j in f:
            LCprobs_di, LCbest_di, s_perc, t_perc = readMCMC(filename1,2)
            sigma=f2[str(j)+'/sigma'].value
            tau=f2[str(j)+'/tau'].value
            
            # LCAprobs_di, LCAbest_di = readMCMC(filename3,2)
            
            s_in_LC.append(sigma)
            t_in_LC.append(tau)
            
            s_in_LCA.append(sigma)
            t_in_LCA.append(tau)
            

            # print(best_di)
            # print(best_di[int(i)])
            s_out_LC.append(LCbest_di[int(i)][0])
            t_out_LC.append(10 ** LCbest_di[int(i)][1])
            
            # s_out_LCA.append(LCAbest_di[int(i)][0])
            # t_out_LCA.append(LCAbest_di[int(i)][1])
            # print(best_di[j])
            # s_out.append(best_di[i][0])
            # t_out.append(best_di[i][1])
            print("original: " + str(t_perc))
            t_perc = [10**x for x in t_perc]
            print(t_perc)
            s_lower.append(s_perc[2])
            s_upper.append(s_perc[1])
            
            t_lower.append(t_perc[2])
            t_upper.append(t_perc[1])
            
    
    s_lower = np.array(s_lower)
    s_upper = np.array(s_upper)
    
    s_err = np.vstack((s_lower,s_upper))
    
    t_lower = np.array(t_lower)
    t_upper = np.array(t_upper)

    t_err = np.vstack((t_lower,t_upper))
    
    
    plt.errorbar(s_in_LC,s_out_LC, yerr=s_err, color='blue', fmt='o', label='L+C+A S in vs. MCMC S out')
    # plt.scatter(s_in_LCA,s_out_LCA, color='red',label='L+C+A S in vs. MCMC S out')
    plt.xlabel("True Sigma")
    plt.ylabel("MCMC Output Sigma")
    plt.plot([0,.55],[0,.55],linestyle='--')
    plt.title("PG1302 L+C+A Simulated DRW Curves w/o Outlier Sigma in vs. out")
    plt.legend()
    plt.show()

    plt.errorbar(t_in_LC,t_out_LC, yerr=t_err, color='blue', fmt='o', label='L+C+A T in vs. MCMC T out')
    # plt.scatter(t_in_LCA,t_out_LCA, color='red', label='L+C+A T in vs. MCMC T out')
    plt.xlabel("True Tau")
    plt.ylabel("MCMC Output Tau")
    plt.plot([0,1100],[0,1100],linestyle='--')
    plt.title("PG1302 L+C+A Simulated DRW Curves w/o Outlier Tau in vs. out")
    plt.legend()
    plt.show()
    # from astroML.time_series import lomb_scargle
    # from statistics import mean
    # # from OptimalFrequencyGrid import OptimalFrequencyGrid
    # count = 0
    # for i in [1,2,3,4,5,6,7,8,9,10]:
    #     filename=r'C:\Users\noahk\pg1302-research\cluster-data\significant_1000000_run_' + str(i) + '_2.h5'
    #     f=h5py.File(filename,"r")
        
    #     for i in f:
    #         count+=1
    #         print(count)
    #         mag,magErr,MJD,max_peak,y,magErry,MJDy,PLS,logf,omegaSlope = readcurve(i,filename)
    #         peaks.append(max_peak)
    
    #         print(max_peak)
    #         a = np.argmax(PLS)
    #         max_f = logf[a]
    #         print(max_f)
    #         frequencies.append(max_f)
            
    
    #         # plt.plot(logf, PLS, '-', c='blue', lw=1, zorder=1)
    #         # plt.xlabel('Log10 of Frequency (Hz)')
    #         # plt.ylabel('P_LS')
    #         # plt.xlim(-8.8,-6.9)
    #         # plt.title("LS Periodogram of LINEAR+CRTS DRW Curve")
    #         # # plt.legend([plot1],["P_LS of LINEAR+CRTS DRW Curve"])
    #         # plt.show()
    #         # plt.errorbar(MJDy, y, yerr=magErry, ecolor='black', fmt='o', color='blue', zorder=2.0, label='L+C Component')
    #         # plt.errorbar(MJD, mag, yerr=magErr, ecolor='black', fmt='o', color='red', zorder=1.0, label='ASASSN Component')
    #         # plt.title("Significant Simulated DRW Curve " + str(count))
    #         # plt.xlabel('MJD (Days)')
    #         # plt.ylabel('Magnitude')
    #         # plt.legend()
    #         # plt.savefig("significant_simulated_DRW_curve" + str(count) + ".png")
    #         # plt.show()
            
    
    #         A_PLS = lomb_scargle(MJD, mag, magErr, omegaSlope, generalized=True)
    # #        P_LS = lomb_scargle(MJD, x, magErr, omegaSlope, generalized=True)
    #         A_frequency = omegaSlope / 6.28 / 86400
    #         A_logf = np.log10(A_frequency)
    # ##        finding the highest peak in the DRW curve's periodogram
    #         # plt.plot(logf, A_PLS, '-', c='blue', lw=1, zorder=1)
    #         # plt.xlabel('Log10 of Frequency (Hz)')
    #         # plt.ylabel('P_LS')
    #         # plt.xlim(-8.8,-6.9)
    #         # plt.title("LS Periodogram of L+C+A DRW Curve " + str(count))
    #         # # plt.legend([plot1],["P_LS of LINEAR+CRTS DRW Curve"])
    #         # plt.savefig("significant_curve_LCA_periodogram" + str(count) + ".png")
    #         # plt.show()
    # #        
    #         A_max_peak = np.amax(A_PLS)
    #         b = np.argmax(A_PLS)
    #         A_max_f = A_logf[b]
    #         A_peaks.append(A_max_peak)
    #         A_frequencies.append(A_max_f)
            
            
    
    # # orig_mag,orig_magErr,orig_MJD = readcurve(0,'1_t0_test.h5')
    # # plt.errorbar(orig_MJD, orig_mag, yerr=orig_magErr, fmt='o', color='red', ecolor='black')
    # # plt.show()
    # #flat_f = []
    # #
    # #for arr in frequencies:
    # #    arr.tolist()
    # #    for item in arr:
    # #        flat_f.append(item)
    # #
    # # plt.hist(peaks, 20)
    # # plt.xlabel('Peak Magnitude')
    # # plt.ylabel('Count')
    # # plt.show()
    
    # # plt.hist(frequencies, 20)
    # # plt.xlabel('Highest Frequency (Log of Hz)')
    # # plt.ylabel('Count')
    # # plt.show()
    
    # # plt.scatter(test_frequencies, test_peaks)
    # # plt.show()
            
    # # from statistics import mean
    # print(mean(peaks))
    # print(mean(A_peaks))
    
    # plt.scatter(-8.20894230751607,.83, color='lawngreen', edgecolors='black', s=100, zorder=113.0, label="PG1302 L+C Only")
    # plt.scatter(-8.222359528836039,.754771676, color='yellow', edgecolors='black', s=100, zorder=113.0, label="PG1302 L+C+A")
    # plt.scatter(frequencies, peaks, color='blue', label="LINEAR+CRTS Curves")
    # plt.scatter(A_frequencies, A_peaks, color='red', label="LINEAR+CRTS+ASASSN Curves")
    # plt.xlabel('Log of Frequency (Hz)')
    # plt.ylabel('Peak Magnitude')
    # plt.legend(prop={'size': 10})
    # plt.show()
    
    # plt.scatter(list(range(1,113)),peaks, color='blue', label="L+C Curves")
    # plt.scatter(list(range(1,113)),A_peaks, color='red', label="L+C+A Curves")
    # plt.scatter(28, .83, color='lawngreen', edgecolors='black', s=100, zorder=113.0, label="PG1302 L+C=.83")
    # plt.scatter(28, .754771676, color='yellow', edgecolors='black', s=100, zorder=113.0, label="PG1302 L+C+A=.755")
    # plt.hlines(mean(peaks), xmin=0, xmax=56, linestyles='dashed', colors=['deepskyblue'], label="L+C Peak Mean=.841")
    # plt.hlines(mean(A_peaks), xmin=0, xmax=56, linestyles='dashed', colors=['lightcoral'], label="L+C+A Peak Mean=.804")
    # plt.xlabel('Curve Number (1 - 56)')
    # plt.ylabel('Peak Magnitude')
    # plt.xlim(0,56)
    # plt.legend(prop={'size': 10})
    # plt.show()
    
    # plt.scatter(list(range(1,113)),peaks, color='blue', label="L+C Curves")
    # plt.scatter(list(range(1,113)),A_peaks, color='red', label="L+C+A Curves")
    # plt.scatter(85, .83, color='lawngreen', edgecolors='black', s=100, zorder=133.0, label="PG1302 L+C=.83")
    # plt.scatter(85, .754771676, color='yellow', edgecolors='black', s=100, zorder=113.0, label="PG1302 L+C+A=.755")
    # plt.hlines(mean(peaks), xmin=57, xmax=112,linestyles='dashed', colors=['deepskyblue'], label="L+C Peak Mean=.841")
    # plt.hlines(mean(A_peaks), xmin=57, xmax=112,linestyles='dashed', colors=['lightcoral'], label="L+C+A Peak Mean=.804")
    # plt.xlabel('Curve Number (57 - 112)')
    # plt.ylabel('Peak Magnitude')
    # plt.xlim(57,112)
    # plt.legend(prop={'size': 10})
    # plt.show()
    
    # plt.scatter(peaks, A_peaks, color='blue', label='DRW Significant Curves')
    # plt.xlabel('L+C Peaks')
    # plt.ylabel('L+C+A Peaks')
    # plt.plot([.82, .9],[.82, .9],linestyle='--', color='black')
    # plt.scatter(.83, .754771676, color="lawngreen", label='PG1302', edgecolors='black', s=100, zorder=113.0)
    # plt.xlim(.82,.9)
    # plt.ylim(.72,.9)
    # plt.legend(prop={'size': 10})
    # plt.show()
    
    # plt.scatter(-8.20894230751607,.83, color='lawngreen', edgecolors='black', s=100, zorder=113.0, label="PG1302 L+C Only")
    # plt.scatter(-8.222359528836039,.754771676, color='yellow', edgecolors='black', s=100, zorder=113.0, label="PG1302 L+C+A")
    # plt.scatter(frequencies, peaks, color='blue', label="LINEAR+CRTS Curves")
    # plt.scatter(A_frequencies, A_peaks, color='red', label="LINEAR+CRTS+ASASSN Curves")
    # plt.xlabel('Log of Frequency (Hz)')
    # plt.ylabel('Peak Magnitude')
    # plt.legend(prop={'size': 10})
    # plt.show()
    
    # plt.scatter(list(range(1,134)),peaks, color='blue', label="L+C Curves")
    # plt.scatter(list(range(1,134)),A_peaks, color='red', label="L+C+A Curves")
    # plt.scatter(28, .83, color='lawngreen', edgecolors='black', s=100, zorder=134.0, label="PG1302 L+C=.83")
    # plt.scatter(28, .754771676, color='yellow', edgecolors='black', s=100, zorder=134.0, label="PG1302 L+C+A=.755")
    # plt.hlines(mean(peaks), xmin=0, xmax=66, linestyles='dashed', colors=['deepskyblue'], label="L+C Peak Mean=.844")
    # plt.hlines(mean(A_peaks), xmin=0, xmax=66, linestyles='dashed', colors=['lightcoral'], label="L+C+A Peak Mean=.805")
    # plt.xlabel('Curve Number (1 - 66)')
    # plt.ylabel('Peak Magnitude')
    # plt.xlim(0,66)
    # plt.legend(prop={'size': 10})
    # plt.show()
    
    # plt.scatter(list(range(1,134)),peaks, color='blue', label="L+C Curves")
    # plt.scatter(list(range(1,134)),A_peaks, color='red', label="L+C+A Curves")
    # plt.scatter(85, .83, color='lawngreen', edgecolors='black', s=100, zorder=134.0, label="PG1302 L+C=.83")
    # plt.scatter(85, .754771676, color='yellow', edgecolors='black', s=100, zorder=134.0, label="PG1302 L+C+A=.755")
    # plt.hlines(mean(peaks), xmin=67, xmax=133,linestyles='dashed', colors=['deepskyblue'], label="L+C Peak Mean=.844")
    # plt.hlines(mean(A_peaks), xmin=67, xmax=133,linestyles='dashed', colors=['lightcoral'], label="L+C+A Peak Mean=.805")
    # plt.xlabel('Curve Number (67 - 133)')
    # plt.ylabel('Peak Magnitude')
    # plt.xlim(67,133)
    # plt.legend(prop={'size': 10})
    # plt.show()
    
    # plt.scatter(peaks, A_peaks, color='blue', label='DRW Significant Curves')
    # plt.xlabel('L+C Peaks')
    # plt.ylabel('L+C+A Peaks')
    # plt.plot([.82, .9],[.82, .9],linestyle='--', color='black')
    # plt.scatter(.83, .754771676, color="lawngreen", label='PG1302', edgecolors='black', s=100, zorder=113.0)
    # plt.xlim(.82,.9)
    # plt.ylim(.72,.9)
    # plt.legend(prop={'size': 10})
    # plt.show()
    #
    # def sinusoid(MJD):
    #     return .137 * np.sin(((2 * np.pi / 2118)*(MJD - 72)))
    
    # def sinusoid2(MJDy):
    #     return .131 * np.sin(((2 * np.pi / 1833)*(MJDy - 21.2)))
    
    # mag,magErr,MJD,omegaSlope = readcurve(0,'final_curve_pg1302_LCA.h5')
    # y,magErry,MJDy,omegaSlopey = readcurve(0,'final_curve_pg1302_LC.h5')
    
    # plt.errorbar(MJD,mag-14.82,yerr=magErr, fmt='o', color='lightcoral',zorder=1)
    # plt.errorbar(MJDy,y-14.82,yerr=magErry, fmt='o', color='deepskyblue',zorder=2)
    # plt.plot(MJD,sinusoid(MJD),color='red',zorder=4)
    # plt.plot(MJDy,sinusoid2(MJDy),color='blue',zorder=3)
    
    # plt.show()
    
    # peaks = []
    # frequencies = []
    # for i in f:
    #     mag,magErr,MJD,max_peak,y,magErry,MJDy,PLS,logf,omegaSlope,frequency = readcurve(i,filename)
    #     peaks.append(max_peak)
    #     m = np.argmax(PLS)
    #     frequencies.append(logf[m])
    #     print(i)
    
    # filename='final_curve34_LC.h5'
    
    # y,magErry,MJDy,omegaSlopey = readcurve(34,filename)
    
    
    # filename='final_curve34_LCA.h5'
    
    # mag,magErr,MJD,omegaSlope = readcurve(34,filename)
    
    # plt.title("Curve 34 DDBIC = 4.22")
    # plt.errorbar(MJD,mag,yerr=magErr,fmt='o',color='red',ecolor='black')
    # plt.errorbar(MJDy,y,yerr=magErry,fmt='o',color='blue',ecolor='black')
    
    # plt.show()
    
    # filename='final_curve44_LC.h5'
    
    # y,magErry,MJDy,omegaSlopey = readcurve(44,filename)
    
    
    # filename='final_curve44_LCA.h5'
    
    # mag,magErr,MJD,omegaSlope = readcurve(44,filename)
    
    # plt.title("Curve 44 DDBIC = 4.64")
    # plt.errorbar(MJD,mag,yerr=magErr,fmt='o',color='red',ecolor='black')
    # plt.errorbar(MJDy,y,yerr=magErry,fmt='o',color='blue',ecolor='black')
    
    # plt.show()
    # data=np.loadtxt(r'C:\Users\noahk\pg1302-research\observation-data\pg1302_data_LINEAR+CRTS+ASASSN.txt')
#     data=np.loadtxt(r'C:\Users\noahk\pg1302-research\observation-data\specific_pg1302_data.txt')
# #    print data
#     MJD=data[:,0]
#     mag=data[:,1]
#     magErr=data[:,2]

#     SortedInd = np.argsort(MJD)
#     MJD = MJD[SortedInd]
#     mag = mag[SortedInd]
#     magErr = magErr[SortedInd]
 
# #    #find frequency grid
#     from OptimalFrequencyGrid import OptimalFrequencyGrid
#     omega, omegaSlope = OptimalFrequencyGrid(MJD)

    
#     #################################################################################################
#     #indices to downsample the uniform grid
#     deltat=1
#     from MinimumTimeResolution import MinimumTimeResolution
#     ind, NumberOfPoints = MinimumTimeResolution(MJD,deltat)
#     #
    
# #    #################################################################################################
# #    #identify best fit parameters for DRW model
# #    import h5py
# #    filename='../QSO1000Iterations'+str(k)+'/QSO1'+str(QSOId)+'_1000Iterations.h5'
# #    f=h5py.File(filename,"r")
# #    BestFitTau=f['/'+str(QSOId)+'/BestFitTau']
# #    BestFitTau=BestFitTau.value
# #    BestFitSigma=f['/'+str(QSOId)+'/BestFitSigma']
# #    BestFitSigma=BestFitSigma.value
# #
# #    BestFitMean=f['/'+str(QSOId)+'/BestFitMean']
# #    BestFitMean=BestFitMean.value
# #    f.close()
# #

# #############################################
# #Fix sigma and tau from Charisi+2015
# #############################################
#     BestFitTau=100
#     BestFitSigma=.11
#     BestFitMean=14
    
    
#     #    from ChiSquareMinimization import ChiSquareMinimization
#     #    BestFitSigma, BestFitTau, BestFitMean, ChiSq=ChiSquareMinimization(mag, magErr, MJD)
#     #    print BestFitSigma
#     #    print BestFitTau
#     #    print BestFitMean
#     #    print ChiSq
    
#     ###############################################################################################
#     # DRW Simulations
#     N_bootstraps=10
#     from SimulateTimeSeries import SimulateTimeSeriesDRW
#     # x=SimulateTimeSeriesDRW(BestFitTau,BestFitSigma,NumberOfPoints,ind,MJD,mag,magErr,BestFitMean,N_bootstraps,deltat)
# #    print x
#     x = mag
#     plt.errorbar(MJD, x, yerr=magErr, fmt='o', ecolor='black')
#     plt.show()

# #    MaxPsim=np.max(Psim,axis=0)
# #    AvgPsim=np.mean(Psim,axis=0)
# #
# #    ################################################################################################
# #    #periodogram
#     from astroML.time_series import lomb_scargle
#     P_LS = lomb_scargle(MJD, x, magErr, omegaSlope, generalized=True)
#     ftemp = omegaSlope / (2*np.pi*86400)
#     logf = np.log10(ftemp)
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
#     h1 = np.log10(ftemp[nearest])
#     h2 = np.log10(ftemp[639])
#     print(h1)
#     print(h2)
# #    print h1
# #    print h2
#     fwhm = h2-h1
# #    print fwhm
#     plt.plot(logf, P_LS, '-', c='black', lw=1, zorder=1)
#     plt.xlabel('Frequency (Hz)')
#     plt.ylabel('P_LS')
# #    plt.xlim(-8.8,-6.9)
#     plt.xlim(-8.8,-6.9)
#     # plt.xlim(0,.000000015)
#     plt.axvline(h1)
#     plt.axvline(h2)
    # plt.title("LS Periodogram of PG1302 L+C Data")
    # plt.scatter(-8.20894230751607,.83, color='lawngreen', edgecolors='black', s=100, zorder=113.0, label="PG1302 L+C Only")
    # plt.scatter(frequencies,peaks,color='blue',label="Simulated DRW's in Frequency Range")
    # plt.xlabel('Log of Frequency (Hz)')
    # plt.ylabel('Max Peak')
    # # plt.legend()
    # plt.show()
    return None


if __name__ == '__main__':
    main()
    
    
    
    
    
    
    ##plt.hist(flat_f, bins=np.arange(min(flat_f), max(flat_f) + .012, .012))
    #plt.hist(flat_f, 5)
    #plt.xlabel('Log10 of Frequency (Hz)')
    #plt.ylabel('Count')
    #plt.show()
    #
    #plt.plot(flat_f, peaks, 'o')
    #plt.xlabel('Log of Frequency (Hz)')
    #plt.ylabel('Peak Magnitude')
    

    
    
