# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 20:36:16 2020

@author: noahk
"""

import time
import SimulatePeriodograms

start = time.time()
SimulatePeriodograms.Simulate(100,"100_drwsin_habanero.h5")
end = time.time()
print("The run took " + str(end - start) + " seconds.")
