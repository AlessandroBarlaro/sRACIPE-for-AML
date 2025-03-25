outdir <- './circuits/'
dir.create(outdir)

fname.eset.brain_array <- '../data.tfs/eset.brain_array.rda'
load(fname.eset.brain_array)

fname.de.results <- '../data.tfs/de.results.rda'
load(fname.de.results)

coreTFs.list <- readRDS('../tfSets/data/coreTFs.rds')
names(coreTFs.list)

targetDB.list <- readRDS('../databases/targetDB.list.rds')

circuits <- list()   
for(fr in names(coreTFs.list)){ 
   print(fr)
   coreTFs.byTOPtfCount <- coreTFs.list[[fr]] 
   for(top.TFs.count in names(coreTFs.byTOPtfCount)){
      print(top.TFs.count) 
      coreTFs.byMethod <- coreTFs.byTOPtfCount[[top.TFs.count]] 
      coreTFs <- coreTFs.byMethod[["COMB"]] # estrarre TF combinati
      
      circuits.tmp <- infer.circuits.by.targetDB(coreTFs = coreTFs,
                                                 targetDB = targetDB.list[[fr]],
                                                 eset.brain_array = eset.brain_array,
                                                 de.results = de.results, 
                                                 int.strengths <- INTERACTION.STRENGTHS, 
                                                 SUBNETWORK.SIZE.TSH=SUBNETWORK.SIZE.TSH)  
      # estarre e salvare i circuiti 
      for(mi.tsh in names(circuits.tmp)){
         idx_name <- paste(fr, top.TFs.count, mi.tsh, sep = '-')
         print(idx_name)
         circuits[[idx_name]] <- circuits.tmp[[mi.tsh]]
      }
   }
}

fname.circuits <- paste0(outdir, 'circuits.rds')
saveRDS(circuits, file = fname.circuits) 

