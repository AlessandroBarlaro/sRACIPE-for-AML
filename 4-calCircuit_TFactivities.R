remove(list = ls()) 

library(NetAct)

# Load circuit summary
#-------------------- 
circuit_metrics <- read.csv(file = './circuits/summary.circuits.csv', row.names = 1)

# select section of the metrics specific to the simulated circuits:
circuit_metrics.sim <- circuit_metrics[!circuit_metrics$DupStatus,]

# Load circuits
circuits <- readRDS('./circuits/circuits.uniq.rds')

# select TOP n TFs from each of the THREE methods: NetAct, MARINa, and RI
#------------------------------------------------------------------------
coreTFs.list <- readRDS('./coreTFs.rds')

# Construct possible circuit set for the TFs based on each TF-target DB in the TF-target DB set
#---------------------------------------------------------------------------------------------- 
targetDB <- readRDS('./targetDB.list.rds')

# Load ExpressionSet object and DE results ---- START 
# Load brain array expression data as an ExpressionSet object
#------------------------------------------------------------
fname.eset.brain_array <- './eset.brain_array.rda'
load(fname.eset.brain_array) # loads object eset.brain_array  

# Load DE results 
#---------------- 
fname.de.results <- './de.results.rda'
load(fname.de.results) # loads object de.results 
# Load ExpressionSet object and DE results ---- END 

# Calculate activities and save them to a file
#--------------------------------------------- 
TF_activities <- list()
circuit_idx <- names(circuits)[1]
for(circuit_idx in names(circuits)){
   print(circuit_idx)
   #break()
   # obtain TF activities using NetAct:
   fr <- strsplit(circuit_idx, split = '-')[[1]][1]
   top.TFs.count <- strsplit(circuit_idx, split = '-')[[1]][2] 
   coreTFs <- coreTFs.list[[fr]][[top.TFs.count]][['COMB']]
   targetDB = targetDB.list[[fr]] 
   length(coreTFs)
   coreTFs <- intersect(coreTFs, names(targetDB)) 
   length(coreTFs) 
   a = TF_Activity(tfs = coreTFs,
                   GSDB = targetDB,  
                   eset = eset.brain_array,
                   DErslt = de.results  #DErslt=de.results$Overall
   )
   #data.REF <- a$all_activities
   TF_activities[[circuit_idx]] <- a$all_activities
   #break()
}

fname.TFactivities <- './circuits/circuit_TFactivities.rds'
saveRDS(TF_activities, fname.TFactivities) 
