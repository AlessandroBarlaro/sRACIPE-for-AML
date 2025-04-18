library(devtools)
install_github("lusystemsbio/NetAct", dependencies=T, build_vignettes = T)

library(igraph)

# prima funzione
# Infer circuits for a set of core TFs based on a given 
# TF-target DB
infer.circuits.by.targetDB <- function(coreTFs, 
                                       targetDB, 
                                       eset.brain_array, 
                                       de.results, 
                                       int.strengths, 
                                       SUBNETWORK.SIZE.TSH){ 
  ## parameters:
  # coreTFs = coreTFs
  # targetDB = targetDB.list[[fr]]  
  # eset.brain_array = eset.brain_array
  # de.results = de.results
  # int.strengths <- INTERACTION.STRENGTHS
  
  ## databases and data structures of this function: 
  # coreTFs
  # targetDB 
  # eset.brain_array 
  # de.results  
  # functions: 
  # createNetworkMiValue
  # retain_uniq_networks 
  # cal.network_metrics.pre_sim 
  # select.largest.connected.subgraph 
  
  # Construct initial network 
  # (1) Construct initial network from TF-target DB
  #------------------------------------------------
  network <- data.frame(Source = rep(names(lapply(targetDB, function(x)names(x))),
                                     times = unlist(lapply(targetDB, function(x)length(x)))),
                        Target = unlist(targetDB))
  
  
  # (2) Retain subnetwork restricted to only core TFs
  #--------------------------------------------------
  # retain TFs in TF-target DB
  #---------------------------
  
  coreTFs <- intersect(coreTFs, names(targetDB)) 
  
  
  # Retain subnetwork restricted to only core TFs
  #----------------------------------------------
  network <- network[intersect(which(network$Source %in% coreTFs),
                               which(network$Target %in% coreTFs)),] 
  
  # (3) Retain subnetwork restricted to activities 
  # Calculate TF activities  
  #------------------------
  a = TF_Activity(tfs = coreTFs,
                  GSDB = targetDB,  
                  eset = eset.brain_array,
                  DErslt = de.results  #DErslt=de.results$Overall
  ) 
  TFactivities <- a$all_activities #a.activities
  
  
  # Retain subnetwork restricted to activities #forse è inutile
  #---------------------------------------------------------
  network <- network[intersect(which(network$Source %in% rownames(TFactivities)),
                               which(network$Target %in% rownames(TFactivities))),]
  
  
  # (4) Add interactions (+ve or -ve type) to the network
  #------------------------------------------------------
  # find network genes
  networkGenes <- union(network$Source, network$Target)
  length(networkGenes)
  
  # subset expression data related to network genes
  TFactivities <- t(TFactivities[networkGenes,])
  dim(TFactivities)
  
  # Infer and assign interaction types
  #------------------------------------ 
  # calculate gene expression correlations
  actCor <- cor(TFactivities, method = "s") 
  dim(actCor)
  #head(rownames(actCor))
  #head(actCor, 2)
  
  # calculate correlation between Source and Target
  int.type <- integer(length = length(network$Source))
  for(i in seq_along(network$Source)){
    int.type[i] <- actCor[as.character(network$Source[i]), as.character(network$Target[i])]
  }
  
  # correlation
  network$Cor <- int.type
  
  # Deterimine interaction type from the Source:Target correlations
  # Assign interaction sign based on correlation
  # Convert correlations to interaction type 
  # >0: 1 (activation)
  # <0: 2 (inhibition)
  int.type[int.type>0] <- 1
  int.type[int.type<=0] <- 2
  network$Type <- int.type
  
  # assign absolute value of correlations to MI
  network$Mi <- abs(network$Cor)
  
  
  # Infer CANDIDATE networks from initial network 
  #----------------------------------------------
  print("Infer CANDIDATE networks")
  networkList  <- list() 
  for(int.strength in int.strengths){
    networkList[[as.character(int.strength)]] <- createNetworkMiValue(tmpNetwork = network, 
                                                                      coreTf = networkGenes, 
                                                                      posMi = int.strength, 
                                                                      negMi = int.strength, 
                                                                      name="net")
  } 
  
  #names(networkList)
  #length(names(networkList))
  
  
  # Remove duplicated networks based on number of nodes and number of interactions 
  # (apply only once)
  #--------------------------------------------------------------------------------
  #length(networkList)
  networkList <- retain_uniq_networks(networkList)
  #length(networkList)
  
  
  # Refine sampled networks: 
  # Retain largest connected component having
  # nodes more than SUBNETWORK.SIZE.TSH  
  ##----------------------------------------- 
  networkList.refined <- list() 
  for(int.strength in names(networkList)){
    print(int.strength) 
    network.cur <- networkList[[int.strength]]  
    network.comp <- select.largest.connected.subgraph(network.cur=network.cur, 
                                                      network = network, 
                                                      subnetwork.size.TSH = SUBNETWORK.SIZE.TSH)
    if(!is.null(network.comp))
      networkList.refined[[int.strength]] <- network.comp
  }
  
  names(networkList.refined)
  length(names(networkList.refined))
  
  return(networkList.refined)
}

#seconda funzione
createNetworkMiValue <- function(tmpNetwork=network, 
                                 coreTf = coreTf, 
                                 posMi = posMi,
                                 negMi = negMi, name = "A549", 
                                 simulate = FALSE, allowAllSameSign = FALSE,
                                 minCoreTFs = 0, minPosInt=0,minNegInt=0){
  # tmpNetwork <- tmpNetwork[union(which(tmpNetwork$Source %in% coreTf),
  #                                which(tmpNetwork$Source %in% coreTf)),]
  # sort the network based on mutual information
  tmpNetwork <- tmpNetwork[order(tmpNetwork$Mi, decreasing = TRUE),]
  # remove any duplicate interactions
  tmpNetwork <- tmpNetwork[!duplicated(tmpNetwork[c("Source","Target","Type")]),]
  # Select the positive network
  posNetwork <- tmpNetwork[tmpNetwork$Cor >0,]
  posNetwork <- posNetwork[which(posNetwork$Mi > posMi),]
  # Select the negative network
  negNetwork <- tmpNetwork[tmpNetwork$Cor < 0,]
  negNetwork <- negNetwork[which(negNetwork$Mi > negMi),]
  
  # Exit if there are either no positive or negative interaction.
  # Can be removed if all positive or all negative networks are ok
  if(!allowAllSameSign){
    if(dim(posNetwork)[1] == 0 | dim(negNetwork)[1] == 0) return()
  }
  # Combine the network
  tmpNetwork <- rbind(posNetwork,negNetwork)
  print("check if there are any duplicate interactions with same or conflicting sign")
  print(sum(duplicated(tmpNetwork[,c("Source","Target")])))
  print(" check if there are any duplicate interactions with same sign")
  print(sum(duplicated(tmpNetwork[,c("Source","Target","Type")])))
  # Remove duplicate interactions
  tmpNetwork <- tmpNetwork[!duplicated(tmpNetwork[,c("Source","Target","Type")]),]
  print("Check if any conflicting interaction is there")
  print(tmpNetwork[duplicated(tmpNetwork[,c("Source","Target")]),])
  print(tmpNetwork[duplicated(tmpNetwork[,c("Source","Target")], fromLast = TRUE),])
  # Remove conflicting interactions
  delRows <- anyDuplicated(tmpNetwork[,c("Source","Target")])
  delRows <- c(delRows, anyDuplicated(tmpNetwork[,c("Source","Target")], fromLast = TRUE))
  delRows <- delRows[delRows>0]
  if(length(delRows>0)){ tmpNetwork <- tmpNetwork[-delRows,]}
  
  # Remove signaling interactions
  # Uncomment for iterative removal
  # interactionRemoved = length(tmpNetwork$Source)
  # while(interactionRemoved>0)
  {
    tmpVar = length(tmpNetwork$Source)
    tmpNetwork <- tmpNetwork[(which(tmpNetwork$Source %in% tmpNetwork$Target )),]
    # if targets only nodes are also to be removed.
    # tmpNetwork <- tmpNetwork[which(tmpNetwork$Target %in% tmpNetwork$Source),]
    # interactionRemoved = tmpVar - length(tmpNetwork$Source)
  }
  
  tmpNetwork <- tmpNetwork[!duplicated(tmpNetwork),]
  #tmpNetwork <- removeNodePairs(tmpNetwork) 
  print("Number of interactions")
  print(length(tmpNetwork$Source))
  networkTfs <- union(tmpNetwork$Source,tmpNetwork$Target)
  #print("Core Tfs included in the network")
  #print(networkTfs[which(networkTfs %in% coreTfGenes)])
  
  #if(length(networkTfs) < minCoreTFs) simulate = FALSE
  
  require(igraph)
  circuit <-  tmpNetwork[,c("Source", "Target", "Type")]
  g <- graph_from_data_frame(circuit, directed = TRUE, vertices = NULL)
   
  return(tmpNetwork) 
}

#terza funzione
# Remove duplicated networks based on number of nodes and number of interactions
#------------------------------------------------------------------------------
retain_uniq_networks <- function(networkList){
  # find unique networks
  network.similarity <- cal.network_metrics.pre_sim(networkList = networkList) 
  
  # sort the networks by three network properties:
  network.similarity.sorted <- network.similarity[with(network.similarity, order(Nodes, Interactions, PosInt)), ] 
  network.similarity <- network.similarity.sorted
  
  node_count.similarity <- duplicated(network.similarity$Nodes)
  sum(node_count.similarity)   
  node_count.similarity   
  interaction_count.similarity <- duplicated(network.similarity$Interactions)
  interaction_count.similarity
  sum(interaction_count.similarity)
  
  posInt.similarity <- duplicated(network.similarity$PosInt) 
  
  # find duplicated networks:
  node_interaction.similarity <- node_count.similarity & interaction_count.similarity & posInt.similarity
  sum(node_interaction.similarity)
  
  # retain the unique networks 
  MIs.retained <- network.similarity$MI[!node_interaction.similarity]
  MIs.retained <- unlist(lapply(MIs.retained, function (x) as.character(x)))
  networkList <-  networkList[MIs.retained]   
  return(networkList)
}

#quarta funzione
cal.network_metrics.pre_sim <- function(networkList){
  networkMetrics <- NULL
  for(mi in names(networkList)){
    #print(mi) 
    circuit <- networkList[[mi]]
    
    posMi <- as.numeric(mi)
    negMi <- as.numeric(mi)
    
    #print(dim(circuit))
    g <- graph_from_data_frame(circuit, 
                               directed = TRUE, 
                               vertices = NULL)
    networkMetricsTmp <- data.frame(MI=as.numeric(mi), 
                                    Nodes= length(union(circuit$Source,circuit$Target)),
                                    Interactions=length(circuit$Source),
                                    PosInt = length(which(circuit$Type ==1)),
                                    Connected = igraph::is_connected(g), 
                                    Transitivity = igraph::transitivity(g),
                                    MeanDistance = igraph::mean_distance(g, directed = TRUE, 
                                                                         unconnected = FALSE), 
                                    MiPos = posMi, 
                                    MiNeg = negMi)  
    networkMetrics <- rbind(networkMetrics, networkMetricsTmp)
  }   
  return(networkMetrics)
}

#quinta funzione
select.largest.connected.subgraph <- function(network.cur,
                                              network,
                                              subnetwork.size.TSH = 0.80){
  #node.total <- union(network$Source, network$Target) 
  node.total <- union(network.cur$Source, network.cur$Target)
  length(node.total)
  g <- igraph::graph_from_data_frame(network.cur, 
                                     directed = TRUE, 
                                     vertices = NULL)    
  clu <- igraph::components(g)
  idx <- which(clu$csize==max(clu$csize))
  
  g.groups <- igraph::groups(clu) 
  
  # if there are more than one subnetworks of the largest size, select the first one:
  nodes.largest_component <- g.groups[[idx[1]]] 
  
  #prop.nodes <- length(nodes.largest_component)/length(node.total)
  
  if(length(nodes.largest_component)/length(node.total) <= subnetwork.size.TSH){
    return() #return(NULL)
  }
  
  status.SOURCE <- network$Source %in% nodes.largest_component
  status.TARGET <- network$Target %in% nodes.largest_component
  status.BOTH <- status.SOURCE & status.TARGET
  network.sele <- network[status.BOTH,] 
  return(network.sele)
}





INTERACTION.STRENGTHS  <- seq(0.00, 0.95, 0.05)
SUBNETWORK.SIZE.TSH <- 0.80 

outdir <- './circuits/'
dir.create(outdir)

download.file(url = "https://github.com/AlessandroBarlaro/sRACIPE-for-AML/blob/main/eset.brain_array.rda", destfile = "./eset.brain_array.rda")                                
fname.eset.brain_array <- './eset.brain_array.rda'
load(fname.eset.brain_array)

download.file(url = "https://github.com/AlessandroBarlaro/sRACIPE-for-AML/blob/main/de.results.rda", destfile = "./de.results.rda")                           
fname.de.results <- './de.results.rda'
load(fname.de.results)

download.file(url = "https://github.com/AlessandroBarlaro/sRACIPE-for-AML/blob/main/coreTFs.rds", destfile = "./coreTFs.rds")                                
coreTFs.list <- readRDS('./coreTFs.rds')

download.file(url = "https://github.com/AlessandroBarlaro/sRACIPE-for-AML/blob/main/targetDB.list.rds", destfile = "./targetDB.list.rds")                                
targetDB.list <- readRDS('./targetDB.list.rds')

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
    # estrarre e salvare i circuiti 
    for(mi.tsh in names(circuits.tmp)){
      idx_name <- paste(fr, top.TFs.count, mi.tsh, sep = '-')
      print(idx_name)
      circuits[[idx_name]] <- circuits.tmp[[mi.tsh]]
    }
  }
}


fname.circuits <- paste0(outdir, 'circuits.rds')
saveRDS(circuits, file = fname.circuits) 

