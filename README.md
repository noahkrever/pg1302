In the functions folder, here are the major functions that I believe are worth chechking out along with their descriptions:

- BIC.py
  * This contains the functions I use to calculate BIC. It uses the results from another function, readMCMC, which returns the max L parameters and likelihood value and uses these values to execute the calculation: BIC = k * ln(n) - 2 * ln(L), and subsequent calculations for delta BIC, and delta delta BIC, as well as print work and plot results with pg1302.

- calculateL.py
  * This contains functions to make the gridsearch/surface plots to find maxL params and values, as well as the MCMC runs for two different (DRW and DRWSIN) models. It uses functions from the file ChiSquareMinimization to get likelihood as formulated by Kozlowski et al. The lnprobs and chains are saved in an h5 file and saved.

- ChiSquareMinimization.py
  * This contains code largely untouched by me that Maria sent. It has functions for getting loglikelihood according to Kozlowski et al, altered slightly for the two models. I have edited these functions to run on log tau, so you will see a 10 ** logtau value where tau was originally.

- readcurves.py
  * This is what I use to quickly return some of the main features of each curve from an h5 file, mag, magErr, MJD, etc.

- readMCMC.py
  * This is the previously mentioned function that unpacks MCMC chains and probabilities and returns the best parameters and max L value.
  
- SimulatePeriodograms.py
  * This is an alteration of PeriodogramRband.py that I now use as a function generate simulated DRW and DRWSIN curves, and it can also be altered to simply save pg1302's data at the different legnths. In its current state it is generating random sigma and tau, no outlier DRW curves at the full L+C+A length.
  
Some other functions of note:
- I use separateCurves.py to separate a large h5 file with many curves into individual curves to be analyzed in parallel on the cluster.
- I am using and will use binpg1302.py and binlightcurves.py to attempt to replicate the Liu et al results (there is currently an issue with the returned magErr array).
- The folders of functions are the very repetitive, iteratively written ones I needed to run things in parallel on the cluster.
  
The rest of the functions are mainly things I needed to simulate different batches of curves for various tests, and a lot of them either involve cut and dry manipulation of h5 files and data, or using the SimulatePeriodograms.py function to simulate something.

Note that all of the functions have various things commented and uncommented, various changes to certain lines and specific files being read in, etc, as I just put this all on github in the somewhat random current state that I had everything in today. Additionally, as I mentioned before, I began this project over a year ago with a knowledge of best practice/comfort in Python far inferior to where it is now with all of the projects and work and classes I've done, so there are some inconsistencies in style which may be confusing.
 
Please reach out if something else is unclear,
 
Thanks!
