# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 18:04:45 2020

@author: noahk
"""

def BIC(curve_file,LC_drw_MCMC_output, LCA_drw_MCMC_output, LC_drwsin_MCMC_output, LCA_drwsin_MCMC_output):
    
    from readMCMC import readMCMC
    from readcurves import readcurve
    import numpy as np
    import h5py
    
    BIC_di = {}
    
    f=h5py.File(curve_file,"r")
    
    LC_drw_logL_di, LC_drw_best_di = readMCMC(LC_drw_MCMC_output,2)
    LC_drwsin_logL_di, LC_drwsin_best_di = readMCMC(LC_drwsin_MCMC_output,5)
    
    LCA_drw_logL_di, LCA_drw_best_di = readMCMC(LCA_drw_MCMC_output,2)
    LCA_drwsin_logL_di, LCA_drwsin_best_di = readMCMC(LCA_drwsin_MCMC_output,5)
    
    for i in f:
        # mag,magErr,MJD,omegaSlope = readcurve(i,curve_file) #adjusted for different lengths already
        
        LC_drw_BIC = 2 * np.log(201) - (2 * LC_drw_logL_di[int(i)])
        # print("LC DRW = 2 * " + str(np.log(201)) + " - (2 * " + str(LC_drw_logL_di[int(i)]) + ") = " + str(LC_drw_BIC))
        LC_drwsin_BIC = 5 * np.log(201) - (2 * LC_drwsin_logL_di[int(i)])
        # print("LC DRWSIN = 5 * " + str(np.log(201)) + " - (2 * " + str(LC_drwsin_logL_di[int(i)]) + ") = " + str(LC_drwsin_BIC))
        
        LCA_drw_BIC = 2 * np.log(406) - (2 * LCA_drw_logL_di[int(i)])
        # print("LCA DRW = 2 * " + str(np.log(407)) + " - (2 * " + str(LCA_drw_logL_di[int(i)]) + ") = " + str(LCA_drw_BIC))
        LCA_drwsin_BIC = 5 * np.log(406) - (2 * LCA_drwsin_logL_di[int(i)])
        # print("LCA DRWSIN = 5 * " + str(np.log(407)) + " - (2 * " + str(LCA_drwsin_logL_di[int(i)]) + ") = " + str(LCA_drwsin_BIC))
        
        LC_drw_del_BIC = LC_drwsin_BIC - LC_drw_BIC
        LCA_drw_del_BIC = LCA_drwsin_BIC - LCA_drw_BIC
        
        LC_drwsin_del_BIC = LC_drw_BIC - LC_drwsin_BIC
        LCA_drwsin_del_BIC = LCA_drw_BIC - LCA_drwsin_BIC
        
        # print("LC DBIC = -1042.1 - -1076.85 = 34.75")
        # print("LCA DBIC = -1754.6 - -1792.69 = 38.09")
        
        # print("DDBIC = 38.09 - 34.75 = 3.34")
        
        BIC_di[int(i)] = LC_drw_del_BIC, LCA_drw_del_BIC,LC_drwsin_del_BIC, LCA_drwsin_del_BIC
        
        
    return BIC_di

def BICpg(curve_file,LC_drw_MCMC_output, LCA_drw_MCMC_output, LC_drwsin_MCMC_output, LCA_drwsin_MCMC_output):
    
    from readMCMC import readMCMC
    from readcurves import readcurve
    import numpy as np
    import h5py
    
    BIC_di = {}
    
    f=h5py.File(curve_file,"r")
    
    LC_drw_logL_di, LC_drw_best_di = readMCMC(LC_drw_MCMC_output,2)
    LC_drwsin_logL_di, LC_drwsin_best_di  = readMCMC(LC_drwsin_MCMC_output,5)
    
    LCA_drw_logL_di, LCA_drw_best_di = readMCMC(LCA_drw_MCMC_output,2)
    LCA_drwsin_logL_di, LCA_drwsin_best_di = readMCMC(LCA_drwsin_MCMC_output,5)
    
    for i in f:
        # mag,magErr,MJD,omegaSlope = readcurve(i,curve_file) #adjusted for different lengths already
        
        LC_drw_BIC = 2 * np.log(201) - (2 * LC_drw_logL_di[int(i)])
        print("LC DRW = 2 * " + str(np.log(201)) + " - (2 * " + str(LC_drw_logL_di[int(i)]) + ") = " + str(LC_drw_BIC))
        LC_drwsin_BIC = 5 * np.log(201) - (2 * LC_drwsin_logL_di[int(i)])
        print("LC DRWSIN = 5 * " + str(np.log(201)) + " - (2 * " + str(LC_drwsin_logL_di[int(i)]) + ") = " + str(LC_drwsin_BIC))
        
        LCA_drw_BIC = 2 * np.log(406) - (2 * LCA_drw_logL_di[int(i)])
        print("LCA DRW = 2 * " + str(np.log(406)) + " - (2 * " + str(LCA_drw_logL_di[int(i)]) + ") = " + str(LCA_drw_BIC))
        LCA_drwsin_BIC = 5 * np.log(406) - (2 * LCA_drwsin_logL_di[int(i)])
        print("LCA DRWSIN = 5 * " + str(np.log(406)) + " - (2 * " + str(LCA_drwsin_logL_di[int(i)]) + ") = " + str(LCA_drwsin_BIC))
        
        LC_drw_del_BIC = LC_drwsin_BIC - LC_drw_BIC
        LCA_drw_del_BIC = LCA_drwsin_BIC - LCA_drw_BIC
        
        LC_drwsin_del_BIC = LC_drw_BIC - LC_drwsin_BIC
        LCA_drwsin_del_BIC = LCA_drw_BIC - LCA_drwsin_BIC
        
        # print("LC DBIC = -1042.1 - -1076.85 = 34.75")
        # print("LCA DBIC = -1754.6 - -1792.69 = 38.09")
        
        # print("DDBIC = 38.09 - 34.75 = 3.34")
        
        BIC_di[int(i)] = LC_drw_del_BIC, LCA_drw_del_BIC,LC_drwsin_del_BIC, LCA_drwsin_del_BIC
        
        
    return BIC_di

def main():
    from matplotlib import pyplot as plt
    from readMCMC import readMCMC
    from readcurves import readcurve
    import numpy as np
    lcdrw = []
    lcadrw = []
    lcdrwsin = []
    lcadrwsin = []
    diff = []

    for i in range(111):
        a,b,c,d,e = ('final_curve' + str(i)+'_LCA_no_outlier.h5',r'C:\Users\noahk\pg1302-research\cluster-data\final_curve' + str(i) + '_LC_drw_info.h5',r'C:\Users\noahk\pg1302-research\cluster-data\final_curve' + str(i) + 'LCA_drw_no_outlier_info.h5', r'C:\Users\noahk\pg1302-research\cluster-data\final_curve' + str(i) + '_LC_drwsin_info.h5', r'C:\Users\noahk\pg1302-research\cluster-data\final_curve' + str(i) + '_LCA_drwsin_no_outlier_info.h5')
        di = BIC(a,b,c,d,e)
        lcdrw.append(di[i][0])
        lcadrw.append(di[i][1])
        lcdrwsin.append(di[i][2])
        lcadrwsin.append(di[i][3])
        diff.append(di[i][3] - di[i][2]) #L+C+A DBIC - L+C DBIC
        # print(di)
        # print(di)
        
    for i in range(112,123):
        a,b,c,d,e = ('final_curve' + str(i)+'_LCA_no_outlier.h5',r'C:\Users\noahk\pg1302-research\cluster-data\final_curve' + str(i) + '_LC_drw_info.h5',r'C:\Users\noahk\pg1302-research\cluster-data\final_curve' + str(i) + 'LCA_drw_no_outlier_info.h5', r'C:\Users\noahk\pg1302-research\cluster-data\final_curve' + str(i) + '_LC_drwsin_info.h5', r'C:\Users\noahk\pg1302-research\cluster-data\final_curve' + str(i) + '_LCA_drwsin_no_outlier_info.h5')
        di = BIC(a,b,c,d,e)
        lcdrw.append(di[i][0])
        lcadrw.append(di[i][1])
        lcdrwsin.append(di[i][2])
        lcadrwsin.append(di[i][3])
        diff.append(di[i][3] - di[i][2]) #L+C+A DBIC - L+C DBIC
        # print(di)
        # print(di)
    
    pgdi = BICpg('pg1302_L+C+A_no_outlier.h5',r'C:\Users\noahk\pg1302-research\cluster-data\final_curve_pg1302_LC_drw_no_outlier_logtau_info_4.h5',r'C:\Users\noahk\pg1302-research\cluster-data\final_curve_pg1302LCA_drw_no_outlier_logtau_info.h5',r'C:\Users\noahk\pg1302-research\cluster-data\final_curve_pg1302_LC_drwsin_no_outlier_logtau_info_4.h5',r'C:\Users\noahk\pg1302-research\cluster-data\final_curve_pg1302_LCA_drwsin_no_outlier_logtau_info.h5')
    # print(pgdi)
    # print("LC drw " + str(readMCMC(r'C:\Users\noahk\pg1302-research\cluster-data\final_curve_pg1302_LC_drw_info.h5',2)))
    # print("LC drwsin " + str(readMCMC(r'C:\Users\noahk\pg1302-research\cluster-data\final_curve_pg1302_LC_drwsin_info.h5',5)))
    # print("LCA drw " + str(readMCMC(r'C:\Users\noahk\pg1302-research\cluster-data\final_curve_pg1302LCA_drw_info.h5',2)))
    # print("LCA drwsin " + str(readMCMC(r'C:\Users\noahk\pg1302-research\cluster-data\final_curve_pg1302_LCA_drwsin_info.h5',5)))
    
    # mag,magErr,MJD,omegaSlope = readcurve(0,'pg1302_LINEAR+CRTS+ASASSN.h5')
    mag,magErr,MJD,omegaSlope = readcurve(0,'pg1302_L+C+A_no_outlier.h5')
    SortedInd = np.argsort(MJD)
    MJD = MJD[SortedInd]
    mag = mag[SortedInd]
    magErr = magErr[SortedInd]
    
    y,magErry,MJDy,omegaSlopey = readcurve(0,'final_curve_pg1302_LC.h5')
    
    # print(len(mag))
    # print(len(y))
    
    
    plt.errorbar(MJD,mag,yerr=magErr, fmt='o')
    plt.show()
    
    plt.errorbar(MJDy,y,yerr=magErry, fmt='o')
    plt.show()
    
    pgdiff = pgdi[0][3] - pgdi[0][2]
    print("DDBIC for pg1302 is " + str(pgdiff))
    
    count = len([i for i in diff if i < pgdiff]) 
    
    print("number of curves less than pg1302 = " + str(count))
    # print([i for i in diff if i < pgdiff])
    print(diff)
    # print("indeces")
    # print([j for j in range(123) if diff[j] > pgdiff])
    
    plt.hist(lcdrw)
    plt.title("LINEAR+CRTS DRW Simulations Delta BIC (DRWSIN - DRW)")
    plt.show()
    
    plt.hist(lcadrw)
    plt.title("LINEAR+CRTS+ASASSN DRW Simulations Delta BIC (DRWSIN - DRW)")
    plt.show()
    
    plt.hist(lcdrwsin)
    plt.title("LINEAR+CRTS DRW Simulations Delta BIC (DRW - DRWSIN)")
    plt.show()
    
    plt.hist(lcadrwsin)
    plt.title("LINEAR+CRTS+ASASSN DRW Simulations Delta BIC (DRW - DRWSIN)")
    plt.show()
    
    plt.scatter(lcdrwsin, lcadrwsin, c='blue')
    plt.plot([-12.5,15],[-12.5,15],'--')
    plt.title("DRW - DRWSIN (Evidence Against DRW) DBICS for Significant DRW Curves")
    plt.xlabel("L+C DRW - DRWSIN")
    plt.ylabel("L+C+A DRW - DRWSIN")
    plt.show()
    
    plt.hist(diff,30)
    plt.axvline(pgdiff)
    plt.title("DDBIC = ( L+C+A [BIC{DRW} - BIC{DRW+SIN}] - L+C [BIC{DRW} - BIC{DRW+SIN}] )")
    plt.show()
    
    
main()