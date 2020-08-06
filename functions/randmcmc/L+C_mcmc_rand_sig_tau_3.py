
import time
import calculateL
        
start = time.time()
calculateL.drw_MCMC(r"/rigel/astro/users/ndk2115/pg1302/data/L+C_mcmc_rand_sig_tau_3.h5",r"/rigel/astro/users/ndk2115/pg1302/outputs/L+C_mcmc_rand_sig_tau_info_3.h5")
end = time.time()
print("The run took " + str(end - start) + " seconds.")
        