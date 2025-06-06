library(sRACIPE)
library(tidyverse)
library(glmnet)

#scaricare la cartella del github contenente tutte le simulazioni skd e dkd
install.packages("gert")
library(gert)
git_clone("https://github.com/AlessandroBarlaro/sRACIPE-for-AML", path = "github")
origine <- "github"
destinazione <- "analisi.kd"
file.rename(from = origine, to = destinazione)


NO_CLUSTERS.REF <- 2

outdir  <- './data/'

outdir.kd.clustered <- './data.sim.clustered/' 
dir.create(outdir.kd.clustered)


# Load trained model 

model <- readRDS('./data/model.glmnet.rds')

# Load training data 

NO_META_COLS = 1
fname_data <- './data/racipe.models.wt.labeled.csv'
mydata = read.csv(file=fname_data, header=TRUE)
dim(mydata)

# Load mean and sd of models (with no knockdown) 

fname_data <- './data/mean_sd.racipe.models.wt.csv'
mean_sd.wt.df = read.csv(file=fname_data, header=TRUE)

means.wt <- mean_sd.wt.df$mean
names(means.wt) <- mean_sd.wt.df$tf
sds.wt <- mean_sd.wt.df$sd
names(sds.wt) <- mean_sd.wt.df$tf

# Load kd simulations 

data.dir.sim <- './analisi.kd/singleknockdown/' 
fnames <- sort(list.files(data.dir.sim))
fnames
length(fnames)
fname <- fnames
# funzioni per estrarre i nomi dei gnei
xtract_gene_name <- function(fname){
  fname.suffix <- strsplit(fname, split = '-', fixed = TRUE)[[1]][2]
  gene_name <- strsplit(fname.suffix, split = '.', fixed = TRUE)[[1]][1]
  return(gene_name)
}

# funzione per la normalizzazione
normalize_by_wt_mean_and_sd <- function(geneExpression, means, sds){
  geneExpression <- sweep(geneExpression, 1, means, FUN = "-")
  geneExpression <- sweep(geneExpression, 1, sds, FUN = "/")  
  return(geneExpression)
}

gene_names <- sapply(fnames, function(fname) xtract_gene_name(fname)) 
gene_names
length(gene_names)

cluster_props <- as.data.frame(matrix(nrow = (length(fnames)+1), ncol = NO_CLUSTERS.REF)) 

rownames(cluster_props) <- c('Untreated', gene_names)
colnames(cluster_props) <- c('Cluster_1', 'Cluster_2')


fname <- fnames[1]
for(fname in fnames){
  print(fname)
  gene_name <- xtract_gene_name(fname)
  racipe <- readRDS(file = paste(data.dir.sim, fname, sep = ''))   
  
  geneExpression.kd <- assay(racipe,1) 
  geneExpression.kd <- log2(1+geneExpression.kd) 
  
  ndata.kd <- t(normalize_by_wt_mean_and_sd(geneExpression=geneExpression.kd, 
                                            means = means.wt, sds=sds.wt)) 
  
  probabilities <- model %>% predict(newx = ndata.kd, type='response') 
  predicted.classes <- ifelse(probabilities > 0.5, 1, 2) 
  
  # obtain cluster proportions: 
  cluster_props[gene_name,] <- c(sum(predicted.classes==1), 
                                 sum(predicted.classes==2)) 
  
  # save ndata.kd with cluster id attached
  ndata.kd.clustered <- cbind(predicted.classes, ndata.kd)  
  colnames(ndata.kd.clustered) <- c('CLUSTER_NO', colnames(ndata.kd))
  fname.ndata <- paste(outdir.kd.clustered, 'racipe-', gene_name, '.csv', sep = '')
  write.csv(ndata.kd.clustered, fname.ndata, quote = F, row.names = F)
  # break()
}

# Calculate cluster proportions for WT models using nnet preditor

probabilities <- model %>% predict(newx = as.matrix(mydata[, 2:ncol(mydata)]), type='response') 
predicted.classes <- ifelse(probabilities > 0.5, 1, 2) # WT models
cluster_props["Untreated",] <- c(sum(predicted.classes==1), 
                                 sum(predicted.classes==2))   

fname.out <- paste(outdir , 'cluster_props.csv', sep = '')
write.csv(cluster_props, file = fname.out, quote = F)



