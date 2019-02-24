# model_pcr_error

model_pcr_error simulates the creation of 'new sequence' due to PCR errors  
It requires a few arguments, but there are defaults if you just want to see how it works  

Required parameters are:  
fragment_length: the length (bp) of PCR product you wish to simulate, default 200  
error_rate: the error rate (per bp/per duplication) of the polymerase, default 0.0001  
pcr_cycles: the number of pcr cycles you wish to simulate. Note, no exhaustion of reagents is currently modeled. Default, 20  
starting_copies: number of starting copies of targeted template, default 10  
stochasticity_threshold: the number of copies below which doubling per cycle is not guaranteed due to fluid dynamics, default 20  
maxseqmodeled: once the copy number gets above this threshold for a given haplotype, error is no longer explicitly modeled, increasing computational speed, default 10000  
replicates: the number of replicate simulated PCRS you wish to run, default 1 

outputs are graphs of counts of wild type versus PCR error sequences by PCR cycle, and the simulated data  

To use this code, copy the R function from https://github.com/laninsky/model_pcr_error/blob/master/model_pcr_error.R, paste it into your R console or R studio,  and then run by:
```
model_pcr_error (fragment_length,error_rate,pcr_cycles,starting_copies,stochasticity_threshold, maxseqmodeled,replicates)

#e.g. to model 366 bp fragment, 100 starting copies, and 10 replicates (w/ defaults for everything else
model_pcr_error(366,starting_copies = 100,replicates=10)
```
