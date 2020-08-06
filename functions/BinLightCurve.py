def BinLightCurve(MJD,mag,magErr):
    import numpy as np
    #MJD changes at midnight in Greenwich. We want to make sure that observations taken at the same night in
    #Palomar are binned together (UT-8)! E.g. observations can be from 5-ish in the evening, which corresponds to 1+-1 in UT
    #to the next morning, say 6-ish, which is in the afternoon in UT. So we need to bin observations with the same MJD
    #same day observations
    #    print np.shape(MJD)
    MJD2=np.floor(MJD);
    MJD_Unique, MJD_UnIndices = np.unique(MJD2, return_inverse=True)
    
    #    print np.shape(MJD_Unique)
    MJD_temp=np.zeros(len(MJD_Unique))
    mag_temp=np.zeros(len(MJD_Unique))
    magErr_temp=np.zeros(len(MJD_Unique))

#    weighted average based on photometric error____MAKE SURE THIS IS CORRECT
    for i in np.arange(len(MJD_Unique)):
        iii=np.where(MJD2==MJD_Unique[i])
        W=np.sum(1.0/magErr[iii]**2)
        wei=1.0/(W*magErr[iii]**2)
    
        MJD_temp[i]=np.average(MJD[iii], axis=None, weights=wei, returned=False)
        mag_temp[i]=np.average(mag[iii], axis=None, weights=wei, returned=False)
        #        magErr_temp[i]=np.average(magErr[iii], axis=None, weights=wei, returned=False)
        magErr_temp[i]=np.sqrt(np.sum(wei**2*magErr[iii]**2))

    MJD=MJD_temp
    mag=mag_temp
    magErr=magErr_temp

    return MJD,mag,magErr




def BinLightCurveFiveDayInterval(MJD,mag,magErr, dt):
    import numpy as np
    #MJD changes at midnight in Greenwich. We want to make sure that observations taken at the same night in
    #Palomar are binned together (UT-8)! E.g. observations can be from 5-ish in the evening, which corresponds to 1+-1 in UT
    #to the next morning, say 6-ish, which is in the afternoon in UT. So we need to bin observations with the same MJD
    #same day observations
    #    print np.shape(MJD)
    MJD2=np.floor(MJD);
    # print np.max(MJD2)+5
    # print np.ceil((np.max(MJD2)+5-np.min(MJD2))/5)
    LinearGrid=np.arange(np.ceil((np.max(MJD2)+dt+1-np.min(MJD2))/dt))*dt+np.min(MJD2)
    
    # print LinearGrid
#    LinearGrid=np.linspace(np.min(MJD2),np.max(MJD2),np.ceil((np.max(MJD2)-np.min(MJD2))/5))

    MJD_Unique, MJD_UnIndices = np.unique(MJD2, return_inverse=True)
    
    #    print np.shape(MJD_Unique)
    MJD_temp=np.zeros(len(LinearGrid)-1)
    mag_temp=np.zeros(len(LinearGrid)-1)
    magErr_temp=np.zeros(len(LinearGrid)-1)
    #    print MJD2
    #    weighted average based on photometric error____MAKE SURE THIS IS CORRECT
    for i in np.arange(len(LinearGrid)-1):
        #        print i
        iii=np.where((MJD2>=LinearGrid[i])&(MJD2<LinearGrid[i+1]))
        #        print np.shape(iii)
        if np.shape(iii)[1]>=1:
            W=np.sum(1.0/magErr[iii]**2)
            wei=1.0/(W*magErr[iii]**2)



            MJD_temp[i]=np.average(MJD[iii], axis=None, weights=wei, returned=False)
            mag_temp[i]=np.average(mag[iii], axis=None, weights=wei, returned=False)
    #        magErr_temp[i]=np.average(magErr[iii], axis=None, weights=wei, returned=False)
            magErr_temp[i]=np.sqrt(np.sum(wei**2*magErr[iii]**2))
        else:
            MJD_temp[i]=-5
            mag_temp[i]=-5
            magErr_temp[i]=-5

#    print MJD_temp
#    print mag_temp
#    print magErr_temp
    kk=np.where(MJD_temp==-5)
    MJD=np.delete(MJD_temp,kk)
    kk=np.where(mag_temp==-5)
    mag=np.delete(mag_temp,kk)
    kk=np.where(magErr_temp==-5)
    magErr=np.delete(magErr_temp,kk)
    
    #    print MJD
    
    return MJD,mag,magErr



def BinLightCurve_Epoch(MJD,mag,magErr,dt_min):
    import numpy as np

    SortedInd = np.argsort(MJD)
    MJD = MJD[SortedInd]
    mag = mag[SortedInd]
    magErr = magErr[SortedInd]

    dt=np.diff(MJD)
    ii=np.where(dt>dt_min)

    epoch_left=np.append(0,ii[0]+1)
    epoch_right=np.append(ii[0],len(MJD)-1)

    
    #    print np.shape(MJD_Unique)
    mjd_epoch=np.zeros(len(epoch_left))
    mag_epoch=np.zeros(len(epoch_left))
    std_epoch=np.zeros(len(epoch_left))
    
    mjd_bin=np.zeros(len(epoch_left))
    mag_bin=np.zeros(len(epoch_left))
    mag_err_bin=np.zeros(len(epoch_left))
    
    #    weighted average based on photometric error____MAKE SURE THIS IS CORRECT
    for kk in np.arange(len(epoch_left)):
     
        
        mjd_epoch[kk]=np.mean(MJD[epoch_left[kk]:epoch_right[kk]]);
        mag_epoch[kk]=np.mean(mag[epoch_left[kk]:epoch_right[kk]]);
        std_epoch[kk]=np.std(mag[epoch_left[kk]:epoch_right[kk]]);
        mjd_temp=MJD[epoch_left[kk]:epoch_right[kk]];
        mag_temp=mag[epoch_left[kk]:epoch_right[kk]];
        mag_err_temp=magErr[epoch_left[kk]:epoch_right[kk]];
        
        ii1=np.where(mag_temp>mag_epoch[kk]+2*std_epoch[kk]);
        ii2=np.where(mag_temp<mag_epoch[kk]-2*std_epoch[kk]);
        ii=np.append(ii1[0],ii2[0])
        mjd_temp=np.delete(mjd_temp,ii)
        mag_temp=np.delete(mag_temp,ii)
        mag_err_temp=np.delete(mag_err_temp,ii)



        if len(mjd_temp)<2:
            mjd_bin[kk]=-9;
            mag_bin[kk]=-9;
            mag_err_bin[kk]=-9;
        else:
            W=sum(1.0/mag_err_temp**2);
            wei=1.0/(W*mag_err_temp**2);
    
            mjd_bin[kk]=np.sum(wei*mjd_temp);
            mag_bin[kk]=np.sum(wei*mag_temp);
            mag_err_bin[kk]=np.sqrt(np.sum(wei**2*mag_err_temp**2));

    mjd_bin=np.delete(mjd_bin,np.where(mjd_bin==-9))
    mag_bin=np.delete(mag_bin,np.where(mag_bin==-9))
    mag_err_bin=np.delete(mag_err_bin,np.where(mag_err_bin==-9))


    return mjd_bin,mag_bin,mag_err_bin



