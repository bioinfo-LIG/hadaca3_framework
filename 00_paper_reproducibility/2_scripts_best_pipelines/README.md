# Find here the best method recomposed
All methods are writen with codabench format and were resubmited to the data challenge. 

## Joker method : 
This is the method submitted by the winning team of Hadaca3. 
It is called the joker method in memory of the team J that posted this method. 

## Best late integration methods: 

* benchmark_late_integration_limeain 
    RNA(LogNorm, scpseudobulk, lm) - 
    scRNA(scConcat,scpseudobulk) - 
    MET(ID, maxdiscriminant, RLRpoisson)-
    Late_integration (limean)

* benchmark_late_integration_TunedJ
    RNA(Scale, ID,RLRpoisson) - 
    MET(ID,maxdiscriminant,RLRpoisson) - 
    Late_integration(TunedJ)

* benchmark_late_integration_limeanRMSE:
    RNA(LogNorm, Toastbulknbfs, RLRpoisson) -
    MET(ID, maxdiscriminant, RLRpoisson) -
    Late_integration(limeanRMSE)

* benchmark_unimodal_onlyMet:
   Met(ID, maxdiscriminant, RLRpoisson)

* benchmark_unimodal_onlyRNA:
    RNA(LogNorm, Toastbulknbfs, RLRpoisson)


## Best Early integration methods : 

* benchmark_early_integration_concatnoscale 
    RNA(ID, scpseudobulk) -
    scRNA(sccluster,  SCcluster) -
    MET(ID, mostmethylated) -
    Early_integration(concatnoscale) -
    Decovolution(RLRpoisson)

* benchmark_early_integration_concatscale
    RNA(LogNorm, scpseudobulk) -
    scRNA(scConcat, scpseudobulk) -
    MET(ID, Toastpercent) -
    Early_integration(concatscale) -
    Decovolution(RLR)

* benchmark_early_integration_Kernel
    RNA(LogNorm, Toastbulknbfs) -
    MET(Scale, SPLSDA) -
    Early_integration(Kernel) -
    Decovolution(epic)

* benchmark_early_integration_omicade4bulk
    RNA(Scale, Toastvst) -
    MET(Scale, SPLSDA) -
    Early_integration(omicade4bulk) -
    Decovolution(lm)

* benchmark_early_integration_OT 
    RNA(LogNorm, Toastbulknbfs) -
    MET(LogNorm, ID) -
    Early_integration(OT) -
    Decovolution(RLR)
