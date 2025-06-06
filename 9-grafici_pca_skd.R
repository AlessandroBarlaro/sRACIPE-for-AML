#grafici pca

NO_CLUSTERS.REF <- 2
#PROB_DIFF_THRESHOLD <- 0.05
# XLIMIT <- c(-6, 14)
# YLIMIT <- c(-10, 8)
XLIMIT <- c(-11, 5)
YLIMIT <- c(-4, 5)

# cluster_names <- c('CL0', 'CL1', 'CL2', 'CL3')
# names(cluster_names) <- c('Uncertain', 'AML', 'Hybrid', 'Untreated') 
# color_vector = c("black", "red", "dodgerblue", "orange") # CL0, CL1, CL2, CL3

cluster_names <- c('CL1', 'CL2')
names(cluster_names) <- c('Untreated', 'AML') 
color_vector = c("dodgerblue", "red") # CL1, CL2

#names(cluster_names) <- c('AML', 'Untreated') 
#color_vector = c("red", "dodgerblue") # CL1, CL2

figdir <- './figs.kd.clusters/'
dir.create(figdir)


# Load pca data 

datadir <- './data/' 
#fname.in <- paste(datadir, 'pca.data.rds', sep = '')
#pca.data <- readRDS(file = fname.in)
#pca <- pca.data$pca
NO_META_COLS = 1
fname_data <- './data/racipe.models.wt.labeled.csv'
mydata = read.csv(file=fname_data, header=TRUE)
dim(mydata)
mydata.saved <- mydata

# remove data in CL3 (NULL models)
mydata <- mydata[mydata$CLUSTER_NO!=3, ] 
dim(mydata)
ndata <- mydata[, (NO_META_COLS+1):ncol(mydata)] 
pca <- prcomp(ndata) 

# PC contributions 

eigs <- pca$sdev^2
stat <- rbind(SD = sqrt(eigs),
              Proportion = eigs/sum(eigs),
              Cumulative = cumsum(eigs)/sum(eigs))

# Load clustered knockdown data 

datadir.sim <- './data.sim.clustered/' 
fnames <- sort(list.files(datadir.sim))
fnames
length(fnames)

fname <- fnames[1]
xtract_gene_name <- function(fname){
  fname.suffix <- strsplit(fname, split = '-', fixed = TRUE)[[1]][2]
  gene_name <- strsplit(fname.suffix, split = '.', fixed = TRUE)[[1]][1]
  return(gene_name)
}

gene_names <- sapply(fnames, function(fname) xtract_gene_name(fname)) 
gene_names
length(gene_names)

library(car) 

# WIDTH <- 6 #8
# HEIGHT <- 5 #6
WIDTH <- 6  
HEIGHT <- 5 

fname <- fnames[1] 
for(fname in fnames){
  print(fname)
  # Load normalized and clustered knockdown simulation data 
  
  ndata.kd <- read.csv(file = paste(datadir.sim, fname, sep = ''))
  
  # project the knockdown simulation data on pca from the wt simulations
  
  prData.kd <- scale(ndata.kd[, 2:ncol(ndata.kd)] , pca$center, pca$scale) %*%  pca$rotation
  
  # Plot the projected data on PC1-PC2 from wt simulation data
  
  gene_name <- xtract_gene_name(fname) # extract knockdown pair from file name
  #figname <- paste(figdir, 'clusters_on_PCs.kd-', gene_name,'.pdf', sep = '') 
  figname <- paste(figdir, 'clusters_on_PCs.kd-', gene_name,'-', WIDTH, 'x', HEIGHT,'.pdf', sep = '') 
  pdf(file = figname, width = WIDTH, height = HEIGHT, paper = 'special')
  
  mainstr <- paste(' Models: single knockdown ', '(', gene_name, ')', sep = '')
  
  plot(pca$x[,1], pca$x[,2], #xlim=XLIMIT, ylim=YLIMIT, 
       main=mainstr, col='gray', pch=19, cex=0.05, 
       xlab=paste('PC1 (', format(stat[2,1]*100, digits = 5), '%)' , sep = ''), 
       ylab=paste('PC2 (', format(stat[2,2]*100, digits = 4), '%)' , sep = '')
  ) 
  #for(cluster.no in rev(0:NO_CLUSTERS.REF)){ 
  for(cluster.no in 1:NO_CLUSTERS.REF){
    prCluster=prData.kd[ndata.kd$CLUSTER_NO==cluster.no,]  
    if(is.vector(prCluster)) { 
      points(prCluster[1], prCluster[2], pch=19,cex=0.5, col=color_vector[cluster.no]) 
    }
    else{
      points(prCluster[,1], prCluster[,2],pch=19,cex=0.5, col=color_vector[cluster.no])  
    }
    
    # dataEllipse(prCluster[,1], prCluster[,2],
    #             levels=c(0.8), add=TRUE, col=color_vector[cluster.no+1],
    #             center.pch = 19,center.cex = 1.5, pch=1, cex=0.05, plot.points=TRUE,
    #             fill=TRUE, fill.alpha=0.1)
    # break()
  }
  #legend(-1, (YLIMIT[1]+7), names(cluster_names), color_vector, bty='n') 
  dev.off()
  #break()
}


# Plot legend

WIDTH <- 3
HEIGHT <- 3
figname <- paste(figdir, 'legend-', WIDTH, 'x', HEIGHT,'.pdf', sep = '')
pdf(file = figname, width = WIDTH, height = HEIGHT, paper = 'special')

plot(pca$x[,1], pca$x[,2], xlim=XLIMIT, ylim=YLIMIT, 
     col='white', pch=19, cex=0.005, xlab='PC1', ylab='PC2') 
legend(XLIMIT[1], (YLIMIT[2]-1), names(cluster_names), color_vector, bty='n') 
dev.off()
