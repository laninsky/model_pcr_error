model_pcr_error <- function (fragment_length=200,error_rate=0.0001,pcr_cycles=20,starting_copies=10,stochasticity_threshold=20,maxseqmodeled=10000,replicates=1) {

# Printing info about this function
print("model_pcr_error simulates the creation of 'new sequence' due to PCR errors")
print("It requires a few arguments, but there are defaults if you just want to see how it works")
print("")
print("Required parameters are:")
print("fragment_length: the length (bp) of PCR product you wish to simulate, default 200")
print("error_rate: the error rate (per bp/per duplication) of the polymerase, default 0.0001")
print("pcr_cycles: the number of pcr cycles you wish to simulate. Note, no exhaustion of reagents is currently modeled. Default, 20")
print("starting_copies: number of starting copies of targeted template, default 10")
print("stochasticity_threshold: the number of copies below which doubling per cycle is not guaranteed due to fluid dynamics, default 20")
print("maxseqmodeled: once the copy number gets above this threshold for a given haplotype, error is no longer explicitly modeled, increasing computational speed, default 10000")
print("replicates: the number of replicate simulated PCRS you wish to run, default 1")
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
print("outputs are graphs of counts of wild type versus PCR error sequences by PCR cycle, and the simulated data")

# Loading in tidyverse
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')

# Creating a variable to capture our output
outputdf <- NULL

# For each replicate simulated PCR....
for (j in 1:replicates) {
  # Printing where we are up to
  print(paste("Up to replicate ",j," out of ",replicates," total replicates",sep=""))
  # Creating a variable to record our output in
  tempoutput <- matrix(0,ncol=1,nrow=(pcr_cycles+1))
  # Intializing it with our starting numbers of template
  tempoutput[1,1] <- starting_copies
  # For each PCR cycle...
  for (i in 2:(pcr_cycles+1)) {
    # Printing where we are up tp
    print(paste("Up to PCR cycle ",(i-1)," out of ",pcr_cycles," total cycles",sep=""))
    # Getting the current number of different haplotypes in the dataset
    currentdim <- dim(tempoutput)[2]
    # For each unique haplotype (wild-type or PCR errors)
    for (l in 1:currentdim) {
      # Heuristic for once numbers get big - no point synthesizing new haplotypes b/c
      # they won't be able to appreciably compete
      if (tempoutput[(i-1),l] > maxseqmodeled) {
        # Number from previous cycle is increased by 2, minus number expected to be lost to error
        tempoutput[i,l] <- tempoutput[(i-1),l]*(2-(1-((1-error_rate)^fragment_length)))
      } else {
        # If haplotypes have fewer copies than masequencemodeled, explicitly modeling error for each copy of that sequence
        for (m in 1:tempoutput[(i-1),l]) {
          # Comparing random number to chance a fragment has NO errors
          if (runif(1) >= ((1-error_rate)^fragment_length)) {
            # If an error is simulated, the copy number for that sequence goes up by 1
            tempoutput[i,l] <- tempoutput[i,l]+1
            # And a new haplotype is created
            newhap <- c(rep(0,(i-1)),1,rep(0,(pcr_cycles-i+1)))
            tempoutput <- cbind(tempoutput,newhap)
          } else {
            # That copy of that sequence is duplicated
            tempoutput[i,l] <- tempoutput[i,l]+2
          }
        }  
      }
    }
    # Additional stochasticity step for very low copy number things (b/c of fluid dynamics)
    for (x in 1:dim(tempoutput)[2]) {
      # For each haplotype, if it is below the stochastity threshold, there is a chance
      # proportional to its frequency relative to the threshold that it will not be duplicated this cycle
      if (runif(1) > tempoutput[i,x]/stochasticity_threshold) {
          tempoutput[i,x] <- round(tempoutput[i,x]/2)
      }
    }  
    
    # sorting based on copy number prevalence after last cycle (if more than one error seq present)
    if (dim(tempoutput)[2] > 2) {
      tosort <- tempoutput[dim(tempoutput)[1],-1]
      names(tosort) <- c(1:length(tosort))
      tosort <- sort(tosort,decreasing=TRUE)
      tosubset <- (as.numeric(names(tosort))+1)
      tempoutput <- tempoutput[,c(1,tosubset)]
    
      # Only taking 20 most prevalent haplotypes after each cycle to speed up calculations
      if(dim(tempoutput)[2] > 20) {
        tempoutput <- tempoutput[,1:20]
      }
    }
  }
  # At the end of that replicate, binding the PCR cycles to tempoutput
  tempoutput <- cbind(c(0:20),tempoutput)
  
  # If this replicate has less than 10 haplotypes, filling in additional columns with NAs
  if (dim(tempoutput)[2] < 10) {
    newcols <- matrix(NA,nrow=21,ncol=(10-dim(tempoutput)[2]))
    tempoutput <- cbind(tempoutput,newcols)
  }
  # Binding the results of this replicate - including the replicate id (j) - into the final output variable
  outputdf <- rbind(outputdf,cbind(j,tempoutput[,1:10]))
} 

# Naming the output                                                
colnames(outputdf) <- c("replicate","cycle","wt","error1","error2","error3","error4","error5","error6","error7","error8")

# Converting it into long format for ggplot2                 
longformat <- gather(as.data.frame(outputdf),haplotype,copy_count,wt:error8)

# Reodering the plotting order of the haplotypes
longformat$haplotype <- factor(longformat$haplotype, levels = c("wt","error1","error2","error3","error4","error5","error6","error7","error8"))

# Creating a folder to hold the results of this combination of parameters
outputfolder <- paste("pcr_error_params",fragment_length,error_rate,pcr_cycles,starting_copies,stochasticity_threshold,replicates,sep="_")

dir.create(outputfolder)

# Plotting by haplotype, facet wrapped by replicate  
ggplot(longformat,aes(x=cycle,y=copy_count,color=haplotype)) + 
  geom_line(size=1.5) + facet_wrap(~ replicate, ncol=1) + scale_y_log10()

ggsave(paste(outputfolder,"/by_replicate.pdf",sep=""), width = 8.5, height = 22, units = "in")

# Taking just the wt and three most prevalent error haplotypes
top3seqerrors <- longformat %>% filter(haplotype == "wt"| haplotype == "error1"| haplotype == "error2"| haplotype == "error3")

# Plotting by replicate, facet wrapping by haplotype (for wt and three most common error haplotypes)
ggplot(top3seqerrors,aes(x=cycle,y=copy_count,color=factor(replicate))) + 
  geom_line(size=1.5) + facet_wrap(~ haplotype, ncol=1) + scale_y_log10()

ggsave(paste(outputfolder,"/by_haplotype.pdf",sep=""), width = 8.5, height = 11, units = "in")

# Saving simulation output
write.table(outputdf,paste(outputfolder,"/simulated_data.txt",sep=""),quote=FALSE,row.names = FALSE)

# Giving a summary of the parameters used in this run
print("Your parameters are:")
print(paste("fragment_length:",fragment_length))
print(paste("error_rate:",error_rate))
print(paste("pcr_cycles:",pcr_cycles))
print(paste("starting_copies:",starting_copies))
print(paste("stochasticity_threshold:",stochasticity_threshold))
print(paste("replicates:",replicates))
print("Your results have been written to:")
print(paste("pcr_error_params",fragment_length,error_rate,pcr_cycles,starting_copies,stochasticity_threshold,replicates,sep="_"))

}
