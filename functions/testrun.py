# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 21:38:38 2020

@author: noahk
"""

import time
import SimulatePeriodograms

start = time.time()
SimulatePeriodograms.Simulate(1,r"/rigel/astro/users/ndk2115/pg1302/outputs/significant_test_run.h5")
end = time.time()
print("The run took " + str(end - start) + " seconds.")
