def MinimumTimeResolution(MJD,deltat):
    import numpy as np
    
    N=round((np.max(MJD)-np.min(MJD))/deltat);
    N=N+1
    #    linear time grid with resolution deltat
    Grid=np.linspace(int(np.min(MJD)),int(np.max(MJD)),int(N))
    
    #    find the indices for the Grid array that will give the times in MJD array \
    #    within deltat
    ind=np.searchsorted(Grid,MJD,"left")
    
    NumberOfPoints=len(Grid)
    return ind, NumberOfPoints


