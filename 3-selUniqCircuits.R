remove(list = ls())

outdir <- './circuits/'

circuits <- readRDS(file = './circuits/circuits.rds')
length(names(circuits))

circuit.summary <- read.csv(file = paste0(outdir, './summary.circuits.csv'), row.names = 1 )

circuits.uniq <- list()
# Select the non-duplicate circuits for simulation
for(idx_name in names(circuits)){ 
   if(!circuit.summary[idx_name, "DupStatus"]){  
      circuits.uniq[[idx_name]]  <- circuits[[idx_name]] 
   }
}

length(names(circuits.uniq))

fname.circuits <- paste0(outdir, 'circuits.uniq.rds')
saveRDS(circuits.uniq, file = fname.circuits) 
