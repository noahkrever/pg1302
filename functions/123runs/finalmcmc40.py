
import time
import calculateL
    
start = time.time()
calculateL.drwsin_MCMC(r"/rigel/astro/users/ndk2115/pg1302/data/final_curve40_LC.h5",r"/rigel/astro/users/ndk2115/pg1302/outputs/final_curve40_LC_drwsin_no_outlier_info.h5")
calculateL.drw_MCMC(r"/rigel/astro/users/ndk2115/pg1302/data/final_curve40_LC.h5",r"/rigel/astro/users/ndk2115/pg1302/outputs/final_curve40_LC_drw_no_outlier_info.h5")
calculateL.drwsin_MCMC(r"/rigel/astro/users/ndk2115/pg1302/data/final_curve40_LCA.h5",r"/rigel/astro/users/ndk2115/pg1302/outputs/final_curve40_LCA_drwsin_no_outlier_info.h5")
calculateL.drw_MCMC(r"/rigel/astro/users/ndk2115/pg1302/data/final_curve40_LCA.h5",r"/rigel/astro/users/ndk2115/pg1302/outputs/final_curve40LCA_drw_no_outlier_info.h5")
end = time.time()
print("The run took " + str(end - start) + " seconds.")
    