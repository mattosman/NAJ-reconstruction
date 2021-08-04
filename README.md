This folder contains the collection of MATLAB scripts needed to run the North Atlantic Jet Stream (NAJ) ice core-based reconstruction described in Osman *et al.* (2021). <br>

The main “driving function” is “NAJ_rec_driver.m”.   Users can define parameters that produce *slightly* different versions of the reconstruction at the top of this script via options that define which NAJ target to use (NOAA20C or ERA20C).  Parameters include an option to detrend the observed NAJ prior to calibration; an option to detrend the Greenlandic proxies; and an option to explore different proxy imputation techniques (see Table S2 and discussions in the Supplementary Information of Osman *et al.*, 2021).  It is not recommended to change the “oldYr” and “calibYr" to values other than 2000 and 1900, respectively. Note that missing proxy values have already been imputed for convenience. <br>

Please note that the CCA-based reconstruction enables an option to manually fix the number of CCA models used in each nest of the reconstruction.  To use this option, change ‘Aut’ to ‘Man’ in Line 131 of “NAJ_rec_driver.m”. <br>

Running this script generates a .mat output file (~4-5mB) in the folder “output”.
Last successfully tested in MATLAB 2020b on July 25, 2021.  If running with the cross validation significance testing (i.e., setting “crossValSig = true” on Line 20 of “NAJ_rec_driver.m”), a 13-century reconstruction takes about 4 hours to run on my 2018 MacBook Pro (2.7 GHz Quad-Core Intel Core i7). Without cross validation significance testing (“crossValSig = false”), it takes ~20-30 minutes. <br>

*Reference*: <br>
Osman, M.B., Coats, S., Das, S.B. *et al*. A thirteen-century context for North Atlantic jet stream projections. *PNAS*, *accepted*, 2021.
