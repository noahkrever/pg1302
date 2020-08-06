
import time
import calculateL
        
start = time.time()
calculateL.drw_MCMC(r"/rigel/astro/users/ndk2115/pg1302/data/L+C+A_mcmc_rand_sig_tau_no_outlier_2.h5",r"/rigel/astro/users/ndk2115/pg1302/outputs/L+C+A_mcmc_rand_sig_tau_no_outlier_info_2.h5")
end = time.time()
print("The run took " + str(end - start) + " seconds.")
        