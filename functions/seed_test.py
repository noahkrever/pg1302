# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 21:05:19 2020

@author: noahk
"""

import time
import SimulatePeriodograms

start = time.time()
SimulatePeriodograms.Simulate(5,"seed_test_2.h5")
end = time.time()
print("The run took " + str(end - start) + " seconds.")
