
import time
import calculateL
    
start = time.time()
calculateL.drwsin_MCMC(r"/rigel/astro/users/ndk2115/pg1302/data/final_curve71_LC.h5",r"/rigel/astro/users/ndk2115/pg1302/outputs/final_curve71_LC_drwsin_no_outlier_info.h5")
calculateL.drw_MCMC(r"/rigel/astro/users/ndk2115/pg1302/data/final_curve71_LC.h5",r"/rigel/astro/users/ndk2115/pg1302/outputs/final_curve71_LC_drw_no_outlier_info.h5")
calculateL.drwsin_MCMC(r"/rigel/astro/users/ndk2115/pg1302/data/final_curve71_LCA.h5",r"/rigel/astro/users/ndk2115/pg1302/outputs/final_curve71_LCA_drwsin_no_outlier_info.h5")
calculateL.drw_MCMC(r"/rigel/astro/users/ndk2115/pg1302/data/final_curve71_LCA.h5",r"/rigel/astro/users/ndk2115/pg1302/outputs/final_curve71LCA_drw_no_outlier_info.h5")
end = time.time()
print("The run took " + str(end - start) + " seconds.")
    