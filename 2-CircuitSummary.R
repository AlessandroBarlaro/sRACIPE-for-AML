remove(list = ls())

#prima funzione
# Calculates specific features of a list of circuits 
# indexed by a conjugate key "no_top_tfs"-"feature_ratio_cutoff"-"absCor"
cal.circuit_metrics.globally <- function(networkList){ 
  # create a matrix to save features of each circuit: 
  circuit.features <- c('FeatureRatio', 'TopTFs', 'AbsCor', 
                        'Nodes', 'Interactions', 'PosInt', 
                        'Connected', 'Transitivity', 'MeanDistance')
  
  networkMetrics <- as.data.frame(matrix(nrow = length(circuits), ncol = length(circuit.features)))
  colnames(networkMetrics) <-  circuit.features   
  rownames(networkMetrics) <- names(circuits) 
  
  # loop through each circuit and calculate the features of the circuit:
  for(circuit_idx in names(networkList)){
    feature_ratio_cutoff  <- strsplit(circuit_idx, '-', 2)[[1]][1] 
    no_top_tfs <- strsplit(circuit_idx, '-', 2)[[1]][2] 
    absCor <- strsplit(circuit_idx, '-', 2)[[1]][3] 
    
    circuit <- networkList[[circuit_idx]]
    
    g <- graph_from_data_frame(circuit, 
                               directed = TRUE, 
                               vertices = NULL) 
    stat.vector <- c(feature_ratio_cutoff = as.numeric(feature_ratio_cutoff),
                     no_top_tfs = as.numeric(no_top_tfs),  
                     absCor = as.numeric(absCor), 
                     Nodes = length(union(circuit$Source,circuit$Target)),
                     Interactions = length(circuit$Source),
                     PosInt = length(which(circuit$Type ==1)),
                     Connected = igraph::is_connected(g), 
                     Transitivity = igraph::transitivity(g),
                     MeanDistance = igraph::mean_distance(g, directed = TRUE, 
                                                          unconnected = FALSE) 
    )
    networkMetrics[circuit_idx, ] <- stat.vector
  } 
  return(networkMetrics)
}

library(igraph)

outdir <- './circuits/'

circuits <- readRDS(file = './circuits/circuits.rds')
length(names(circuits)) 


circuit_metrics <- cal.circuit_metrics.globally(networkList=circuits)
circuit_metrics.tmp <- circuit_metrics


circuit.metrics.sorted <- circuit_metrics[with(circuit_metrics, order(Nodes, Interactions, PosInt)), ]
circuit_metrics <- circuit.metrics.sorted

max(circuit.metrics.sorted$Nodes)
min(circuit.metrics.sorted$Nodes)
sum(circuit.metrics.sorted$Nodes==4)


# Find duplicate topologies - based on THREE network features 
# node_count, interaction_count.similarity, posInt.similarity: 
# screening by metric node_count 
node_count.similarity <- duplicated(circuit_metrics$Nodes)
sum(node_count.similarity)  # 533 
 
# screening by metric interaction_count.similarity
interaction_count.similarity <- duplicated(circuit_metrics$Interactions) 
interaction_count.similarity
sum(interaction_count.similarity) # 461

# screening by metric posInt.similarity
posInt.similarity <- duplicated(circuit_metrics$PosInt)  
sum(posInt.similarity) # 517
 

# Find duplicates
node_interaction_posInt.similarity <- node_count.similarity & interaction_count.similarity & posInt.similarity
sum(node_interaction_posInt.similarity) # 396

circuit_metrics$DupStatus <- node_interaction_posInt.similarity


# Add a column to indicate which circuit is simulated
SimIdx <- vector(mode = "character", length = dim(circuit_metrics)[1])
SimIdx[1] <- rownames(circuit_metrics)[1]
for(j in 2:dim(circuit_metrics)[1]){
   if(circuit_metrics$DupStatus[j]) SimIdx[j] <- SimIdx[j-1] 
   else SimIdx[j] <- rownames(circuit_metrics)[j] 
}
circuit_metrics$SimIdx <- SimIdx

# Attach the topology id (row name) of the non-duplicated topology 
# to the duplicated topologies with that non-duplicated topology
#---------------------------------------------------------------
write.csv(circuit_metrics, file = paste0(outdir, './summary.circuits.csv'), 
          row.names = T, quote = F)
write.csv(circuit_metrics.tmp, file = paste0(outdir, './summary.circuits.asorted.csv'), 
          row.names = T, quote = F)

