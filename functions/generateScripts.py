# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 19:12:20 2020

@author: noahk
"""

# for i in range(123):
#     output1 = open("finalmcmc" + str(i) + ".sh", "w")
    
#     output1.write(
# """#!/bin/sh
# #SBATCH --account=astro
# #SBATCH --job-name= finalmcmc""" + str(i) + """
# #SBATCH --mail-type ALL
# #SBATCH --mail-user ndk2115@columbia.edu
# #SBATCH -c 1
# #SBATCH --time=0-04:00
# #SBATCH --mem-per-cpu=1gb
    
# module load anaconda
# source activate /rigel/astro/users/ndk2115/python-envs/pg
# python3 finalmcmc""" + str(i) + """.py""")
    
# output2 = open("finalmcmcpg1302_no_outlier.py", "w")

# output2.write(
#     """
# import time
# import calculateL
    
# start = time.time()
# calculateL.drwsin_MCMC(r"/rigel/astro/users/ndk2115/pg1302/data/final_curve_pg1302_LC.h5",r"/rigel/astro/users/ndk2115/pg1302/outputs/final_curve_pg1302_LC_drwsin_info_2.h5")
# calculateL.drw_MCMC(r"/rigel/astro/users/ndk2115/pg1302/data/final_curve_pg1302_LC.h5",r"/rigel/astro/users/ndk2115/pg1302/outputs/final_curve_pg1302_LC_drw_info_2.h5")
# calculateL.drwsin_MCMC(r"/rigel/astro/users/ndk2115/pg1302/data/pg1302_L+C+A_no_outlier.h5",r"/rigel/astro/users/ndk2115/pg1302/outputs/final_curve_pg1302_LCA_drwsin_no_outlier_info.h5")
# calculateL.drw_MCMC(r"/rigel/astro/users/ndk2115/pg1302/data/pg1302_L+C+A_no_outlier.h5",r"/rigel/astro/users/ndk2115/pg1302/outputs/final_curve_pg1302LCA_drw_no_outlier_info.h5")
# end = time.time()
# print("The run took " + str(end - start) + " seconds.")
#     """)
# output2.close()

# output1 = open("finalmcmcpg1302_no_outlier.sh", "w")
    
# output1.write(
# """#!/bin/sh
# #SBATCH --account=astro
# #SBATCH --job-name= finalmcmcpg1302_no_outlier
# #SBATCH --mail-type ALL
# #SBATCH --mail-user ndk2115@columbia.edu
# #SBATCH -c 1
# #SBATCH --time=0-02:00
# #SBATCH --mem-per-cpu=1gb
    
# module load anaconda
# source activate /rigel/astro/users/ndk2115/python-envs/pg
# python3 finalmcmcpg1302_no_outlier.py""")

# output1.close()

for i in range(10):
    output1 = open("L+C+A_mcmc_rand_sig_tau_" + str(i) + ".sh", "w")
        
    output1.write(
    """#!/bin/sh
#SBATCH --account=astro
#SBATCH --job-name= L+C+A_mcmc_rand_sig_tau_""" + str(i) + """
#SBATCH --mail-type ALL
#SBATCH --mail-user ndk2115@columbia.edu
#SBATCH -c 1
#SBATCH --time=0-00:45
#SBATCH --mem-per-cpu=1gb
    
module load anaconda
source activate /rigel/astro/users/ndk2115/python-envs/pg
python3 L+C+A_mcmc_rand_sig_tau_""" + str(i) + """.py""")
    
    output1.close()

#     output2 = open("L+C+A_mcmc_rand_sig_tau_" + str(i) + ".py", "w")

#     output2.write(
#         """
# import time
# import calculateL
        
# start = time.time()
# calculateL.drw_MCMC(r"/rigel/astro/users/ndk2115/pg1302/data/L+C+A_mcmc_rand_sig_tau_""" + str(i) + """.h5",r"/rigel/astro/users/ndk2115/pg1302/outputs/L+C+A_mcmc_rand_sig_tau_info_""" + str(i) + """.h5")
# end = time.time()
# print("The run took " + str(end - start) + " seconds.")
#         """)
#     output2.close()