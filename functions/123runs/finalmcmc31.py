
import time
import calculateL
    
start = time.time()
calculateL.drwsin_MCMC(r"/rigel/astro/users/ndk2115/pg1302/data/final_curve31_LC.h5",r"/rigel/astro/users/ndk2115/pg1302/outputs/final_curve31_LC_drwsin_no_outlier_info.h5")
calculateL.drw_MCMC(r"/rigel/astro/users/ndk2115/pg1302/data/final_curve31_LC.h5",r"/rigel/astro/users/ndk2115/pg1302/outputs/final_curve31_LC_drw_no_outlier_info.h5")
calculateL.drwsin_MCMC(r"/rigel/astro/users/ndk2115/pg1302/data/final_curve31_LCA.h5",r"/rigel/astro/users/ndk2115/pg1302/outputs/final_curve31_LCA_drwsin_no_outlier_info.h5")
calculateL.drw_MCMC(r"/rigel/astro/users/ndk2115/pg1302/data/final_curve31_LCA.h5",r"/rigel/astro/users/ndk2115/pg1302/outputs/final_curve31LCA_drw_no_outlier_info.h5")
end = time.time()
print("The run took " + str(end - start) + " seconds.")
    