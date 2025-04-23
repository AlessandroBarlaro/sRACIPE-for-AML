remove(list = ls())

#1)simulazione dei circuiti

library(sRACIPE)
#caricare l'oggetto contenente la topologia dei circuiti
circuits <- readRDS(file = './circuits/circuits.uniq.rds') 

#creare una nuova directory per i circuiti simulati
outdir <- './circuits.sim/'
dir.create(outdir)

set.seed(100) 
#per velocizzare le simulazioni, eseguire la computazione parallela 
library(doParallel)
library(foreach)
NO_AVAILABLE_CORES <- detectCores() 
print(NO_AVAILABLE_CORES)
cl <- makeCluster(NO_AVAILABLE_CORES)
class(cl)
registerDoParallel(cl)
getDoParWorkers()

#eseguire le simulazioni per ogni circuito
rv = foreach(idx_name = names(circuits), .inorder = TRUE) %dopar% { 
  library(sRACIPE)
  circuitSimulated <- sracipeSimulate(circuit = circuits[[idx_name]][,c("Source", "Target", "Type")], 
                                      plotToFile = FALSE,
                                      numModels = 10000, 
                                      plots = FALSE,
                                      stepper = "RK4",
                                      integrateStepSize = 0.05) 
  fname.sim.circuit <- paste(outdir, './circuit_simulated_', idx_name ,'.rds', sep = '') 
  saveRDS(circuitSimulated, file = fname.sim.circuit)
}

stopCluster(cl)

#2)calcolo accuratezza

#creazione nuova directory per i risultati di heatmapsimilarity
outdir.hS <- './circuits.hS/'
dir.create(outdir.hS)

library(NetAct)

#caricare dataframe riepilogativo dei circuiti
circuit_metrics <- read.csv(file = './circuits/summary.circuits.csv', row.names = 1)

#selezionare solamente i circuiti con DupStatus = true
circuit_metrics.sim <- circuit_metrics[!circuit_metrics$DupStatus,]

#rimuovere i circuiti con nodi < 15
circuit_metrics.sim <- circuit_metrics.sim[circuit_metrics.sim$Nodes>=15, ]

#caricare oggetto con i migliori TF calcolati con i tre metodi
coreTFs.list <- readRDS('.coreTFs.rds')

#caricare oggetto possibili target TF
targetDB.list <- readRDS('./targetDB.list.rds')

#caricare le espressioni sperimentali
fname.eset.brain_array <- './eset.brain_array.rda'
load(fname.eset.brain_array)  
 
#caricare oggetto de.results 
fname.de.results <- './de.results.rda'
load(fname.de.results)

#caricare il cluster cut dei dati sperimentali
tmp <- read.csv("https://raw.githubusercontent.com/AlessandroBarlaro/sRACIPE-for-AML/refs/heads/main/clusterCut.REF.csv", row.names = 1)
clusterCut.REF <- as.integer(tmp$x)
names(clusterCut.REF) <- as.character(rownames(tmp))

#si ricavano i data REF e data SIM e poi si esegue la comparazione
for(circuit_idx in rownames(circuit_metrics.sim)){
  
  #utilizzo di NetAct per ottenere l'attività dei TF (data REF)
  fr <- strsplit(circuit_idx, split = '-')[[1]][1]
  top.TFs.count <- strsplit(circuit_idx, split = '-')[[1]][2]  
  coreTFs <- coreTFs.list[[fr]][[top.TFs.count]][['COMB']]
  targetDB = targetDB.list[[fr]] 
  coreTFs <- intersect(coreTFs, names(targetDB))
  length(coreTFs) 
  a = TF_Activity(tfs = coreTFs,
                  GSDB = targetDB,  
                  eset = eset.brain_array,
                  DErslt = de.results  
  )
  data.REF <- a$all_activities
  
  # caricare le simulazioni dei circuiti
  racipe <- readRDS(file = paste('./circuits.sim/circuit_simulated_', 
                                 circuit_idx, '.rds', sep = ''))
  # normalizzare i risultati
  racipe <- sracipeNormalize(racipe) 
  
  # estrazione espressioni simulate (data SIM)
  data.sim <- assay(racipe,1)  
  
  # utilizzo di sracipeHeatmapSimilarity
  hS <- sracipeHeatmapSimilarity(dataReference = data.REF,
                                 dataSimulation = data.sim, 
                                 returnData = T, 
                                 clusterCut = clusterCut.REF)  
  
  # aggiungere i nuovi parametri di similarità nel dataframe
  circuit_metrics.sim[circuit_idx, "ClusterA"] <- hS$simulated.cluster.freq[2] 
  circuit_metrics.sim[circuit_idx, "ClusterB"] <- hS$simulated.cluster.freq[3] 
  circuit_metrics.sim[circuit_idx, "Accuracy"] <- (hS$simulated.cluster.freq[2] + hS$simulated.cluster.freq[3])
  circuit_metrics.sim[circuit_idx, "KLdist"] <- hS$KL 
  circuit_metrics.sim[circuit_idx, "ClusterA2"] <- hS$cluster.similarity[2] 
  circuit_metrics.sim[circuit_idx, "ClusterB2"] <- hS$cluster.similarity[3] 
  #AVGdist è presente nella funzione heatmapsimilarity update, ma non sembra essere necessario
  #circuit_metrics.sim[circuit_idx, "AvgDist"] <- hS$AvgDist 
  
  # salvataggio dei risultati di heatmapsimilarity
  fname.hS <- paste(outdir.hS, 'hS_', circuit_idx, '.rds', sep = '')
  saveRDS(hS, file = fname.hS)
}

#creare excel con nuovi dati su accuratezza
write.csv(circuit_metrics.sim, file = paste0(outdir, "./summary.circuits.sim.csv"),
          row.names = T, quote = F) 

#3)calcolo flessibilità

#si crea un dataframe per i valori di flessibilità
flexibility.df <- as.data.frame(matrix(nrow = length(circuits), ncol = 2)) 
colnames(flexibility.df) <- c('circuit.idx', 'flexibility')
rownames(flexibility.df) <- names(circuits)

#si crea una funzione per calcolare la distanza euclidea tra le distribuzioni 
#delle simulazioni KD e distribuzioni WT e riporta la loro media
calDistance <- function(racipe.kd){
  prop.wt <- racipe.kd$WT
  dist.df <- as.data.frame(matrix(nrow = (length(racipe.kd)-1), ncol = 2)) 
  colnames(dist.df) <- c('tf', 'dist')
  rownames(dist.df) <- names(racipe.kd)[2:length(racipe.kd)]
  for(tf in names(racipe.kd)[2:length(racipe.kd)]){
    #print(tf)
    prop.kd <- racipe.kd[[tf]] 
    d <- sqrt((as.numeric(prop.wt[1])-as.numeric(prop.kd[1]))^2 + 
                (as.numeric(prop.wt[2])-as.numeric(prop.kd[2]))^2)
    dist.df[tf,] <- c(tf, d)
    #break()
  }
  
  avg.dist <- sum(as.numeric((dist.df$dist)))/length(dist.df$dist)   
  return(avg.dist)
}

#si esegue analisi knockdown con sracipeknockdown
for(circuit_idx in names(circuits)){
  print(circuit_idx)
  
  # si caricano i circuiti simulati e si normalizzano
  racipe <- readRDS(file = paste('./circuits.sim/circuit_simulated_', 
                                 circuit_idx, '.rds', sep = '')) 
  racipe <- sRACIPE::sracipeNormalize(racipe) 
  
  # impostazioni per knockdown
  racipe.kd <- sRACIPE::sracipeKnockDown(racipe, plotToFile = FALSE,
                                         plotBarPlot = FALSE, 
                                         plotHeatmap = FALSE, 
                                         reduceProduction = 10
  )
  avg.dist <- calDistance(racipe.kd)
  flexibility.df[circuit_idx, ] <- c(circuit_idx, avg.dist)
  #break()
}

#si crea la sottocartella results e si salva l'excel della flessibilità
outdir <- './results/'
dir.create(outdir)
fname <- paste(outdir, 'circuit.flexibility.csv', sep = '')
write.csv(flexibility.df, file = fname, quote = FALSE, row.names = FALSE)

#4)creare un nuovo dataframe sommario+flessibilità

#si carica il sommario
circuit_metrics.sim <- read.csv(file = './results/summary.circuits.sim.csv', row.names = 1)

#si caricano i risultati della flessibilità
flexibility.df <- read.csv(file = './results/circuit.flexibility.csv')
rownames(flexibility.df) <- flexibility.df$circuit.idx

#si ordina la flessibilità per i nomi dei circuiti del sommario
flexibility.ordered <- flexibility.df[rownames(circuit_metrics.sim), ]

#si controlla che le righe siano uguali
sum(rownames(circuit_metrics.sim) == rownames(flexibility.ordered))

#si crea un nuovo dataframe con tutte le colonne del sommario e la colonna 
#flexibility della flessibilità
circuit_metrics.sim.new <- cbind(circuit_metrics.sim, flexibility.ordered$flexibility)
colnames(circuit_metrics.sim.new) <- c(colnames(circuit_metrics.sim), 'flexibility')

#si salva il dataframe
outdir <- './results/'
dir.create(outdir)
write.csv(circuit_metrics.sim.new, file = paste0(outdir, "./summary.circuits.sim.flex.csv"),
          row.names = T, quote = F)

#salto lo script cal.maxSigValByMethod
#identificazione miglior circuito accuratezza/flessibilità

#caricare sommario+flessibilità
circuit_metrics.sim <- read.csv(file = './results/summary.circuits.sim.flex.csv', row.names = 1)

#funzione per ordinare dataframe per accuratezza e flessibilità
sortByTwoIndices.acc_and_flex <- function(circuit_metrics.sim){
  # accuratezza:
  data.sortedByAccuracy <- circuit_metrics.sim[order(circuit_metrics.sim$Accuracy, decreasing = T), ]
  data.sortedByAccuracy$idxAccuracy <- 1:dim(data.sortedByAccuracy)[1]
  
  # flessibilità:
  data.sortedByFlexibility <- data.sortedByAccuracy[order(data.sortedByAccuracy$flexibility, decreasing = T), ]
  data.sortedByFlexibility$idxFlexibility <- 1:dim(data.sortedByFlexibility)[1] # adding index
  
  # aggiunta degli indici:
  data.sortedByFlexibility$idxBoth <- data.sortedByFlexibility$idxAccuracy+data.sortedByFlexibility$idxFlexibility
  
  # ordine dato dagli indici congiunti:
  data.sortedByBoth <- data.sortedByFlexibility[order(data.sortedByFlexibility$idxBoth, decreasing = F),]   
  return(data.sortedByBoth)
}

#ordinare tramite funzione sortByTwoIndices.acc_and_flex
circuit_metrics.sim.sorted <- sortByTwoIndices.acc_and_flex(circuit_metrics.sim) 
write.csv(circuit_metrics.sim.sorted, 
          file = './results/summary.circuits.sim.sortedByAcc_flex.csv',  
          row.names = T, quote = F) 
