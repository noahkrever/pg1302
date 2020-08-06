
import time
import calculateL
    
start = time.time()
calculateL.drwsin_MCMC(r"/rigel/astro/users/ndk2115/pg1302/data/final_curve_pg1302_LC.h5",r"/rigel/astro/users/ndk2115/pg1302/outputs/final_curve_pg1302_LC_drwsin_info_2.h5")
calculateL.drw_MCMC(r"/rigel/astro/users/ndk2115/pg1302/data/final_curve_pg1302_LC.h5",r"/rigel/astro/users/ndk2115/pg1302/outputs/final_curve_pg1302_LC_drw_info_2.h5")
calculateL.drwsin_MCMC(r"/rigel/astro/users/ndk2115/pg1302/data/pg1302_L+C+A_no_outlier.h5",r"/rigel/astro/users/ndk2115/pg1302/outputs/final_curve_pg1302_LCA_drwsin_no_outlier_info.h5")
calculateL.drw_MCMC(r"/rigel/astro/users/ndk2115/pg1302/data/pg1302_L+C+A_no_outlier.h5",r"/rigel/astro/users/ndk2115/pg1302/outputs/final_curve_pg1302LCA_drw_no_outlier_info.h5")
end = time.time()
print("The run took " + str(end - start) + " seconds.")
    