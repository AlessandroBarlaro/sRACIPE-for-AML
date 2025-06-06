library(sRACIPE)
library(NetAct)
library(gplots)

circuit_idx <- '0.09-32-0.85'
NO_MODELS <- 10000

# scaricare file hS del circuito migliore
download.file(url = "https://github.com/AlessandroBarlaro/sRACIPE-for-AML/raw/refs/heads/main/hS_0.09-32-0.85.rds",
              destfile="./hS_0.09-32-0.85.rds")

outdir <- './data/'
dir.create(outdir)

# scaricare la simualzione del circuito migliore
download.file(url = "https://github.com/AlessandroBarlaro/sRACIPE-for-AML/raw/refs/heads/main/circuit_simulated_0.09-32-0.85.rds", 
              destfile = "./circuit_simulated_0.09-32-0.85")
racipe <- readRDS(file ="./circuit_simulated_0.09-32-0.85")
circuit_tpo <- sracipeCircuit(racipe)

# salvare topologia del circuito

fname.out <- paste(outdir, 'circuit-', circuit_idx, '.tpo' ,sep = '')
write.table(circuit_tpo, file = fname.out, sep = '\t', quote = F, row.names = F)

# salvare i fattori di trascrizione della rete

TFs_in_circuit <- union(circuit_tpo$Source, circuit_tpo$Target)

fname.out <- paste(outdir, 'TFs-', circuit_idx, '.txt' ,sep = '')
write.table(TFs_in_circuit, file = fname.out, sep = '\t', quote = F, row.names = F)

# caricare il file hS del circuito

fname.hS <- "hS_0.09-32-0.85.rds"
hS <- readRDS(file = fname.hS)

hS$simulated.cluster.freq * NO_MODELS
sum(hS$simulated.cluster.freq * NO_MODELS)

# salvare le simulazioni con le informazioni sui cluster

dim(hS$dataSimulation)
data.sim <- t(hS$dataSimulation)
dim(data.sim) 
cluster.labels <- c(rep('1', hS$simulated.cluster.freq[2]*NO_MODELS),  
                    rep('2', hS$simulated.cluster.freq[3]*NO_MODELS), 
                    rep('3', (as.numeric(hS$simulated.cluster.freq[1]*NO_MODELS))+1)) 

length(cluster.labels) 
sum(cluster.labels=='3') # 642 - Hybrid
sum(cluster.labels=='1') # 2443 - AML
sum(cluster.labels=='2') # 6915 - Untreated

data.sim.labeled <- cbind(cluster.labels, data.sim) 
colnames(data.sim.labeled) <- c('CLUSTER_NO', colnames(data.sim))
fname.out <- paste(outdir, 'racipe.models.wt.labeled.csv', sep = '')
write.csv(data.sim.labeled, fname.out, row.names = FALSE, quote = F)







#-----------------------------

geneExpression <- assay(racipe, 1) 
geneExpression <- log2(1+geneExpression)
means <- rowMeans(geneExpression)
sds <-  apply(geneExpression, 1, sd)

# salvare i dati non perturbati
mean_sd.df <- cbind(names(means), means, sds)
colnames(mean_sd.df) <- c('tf' ,'mean', 'sd')

fname.out <- paste(outdir, 'mean_sd.racipe.models.wt.csv', sep='')
write.csv(mean_sd.df, file = fname.out, quote = F, row.names = FALSE)


#-----------------------------

ndata.df <- read.csv(file = "./data/racipe.models.wt.labeled.csv")
mean_sd.df <- read.csv(file = "./data/mean_sd.racipe.models.wt.csv", row.names = 1)
data.df <- as.data.frame(matrix(nrow = nrow(ndata.df), ncol = ncol(ndata.df)))
colnames(data.df) <- colnames(ndata.df)

col.no <- 1
for(gname in colnames(ndata.df)[2:ncol(ndata.df)]){
  print(gname)
  data.df[, gname] <- ndata.df[, gname] * mean_sd.df[gname, "sd"] + mean_sd.df[gname, "mean"] 
  #break()
}

data.df$CLUSTER_NO <- ndata.df$CLUSTER_NO


fname.out <- paste(outdir, 'racipe.models.unnormalized.csv', sep = '')
write.csv(data.df, fname.out, row.names = FALSE, quote = F)  


#------------------------------

library(glmnet)

NO_CLUSTERS.REF <- 2

# Load training data 

NO_META_COLS = 1
fname_data <- './data/racipe.models.wt.labeled.csv'
mydata = read.csv(file=fname_data, header=TRUE)
dim(mydata)
mydata <- mydata[mydata$CLUSTER_NO!=3, ]
dim(mydata)

mydata$CLUSTER_NO <- as.numeric(mydata$CLUSTER_NO==1)
mydata$CLUSTER_NO <- factor(mydata$CLUSTER_NO)
unique(mydata$CLUSTER_NO)

# randomize data 
set.seed(1)
mydata <- mydata[sample(nrow(mydata)), ]

x <- model.matrix(CLUSTER_NO ~., mydata)[,-1] 
y <- ifelse(mydata$CLUSTER_NO == 1, 1, 0)
head(y)


ALPHA <- 0 # Ridge regression

# perform 10 fold cross validation to obtain best lamba (penalty) 
cv.glmnet.obj <- cv.glmnet(x, y, alpha = ALPHA, family = "binomial") 

# create glmnet model using the optimum lambda
model.glmnet <- glmnet(x, y, alpha = ALPHA, family = "binomial",
                       lambda = cv.glmnet.obj$lambda.1se #cv.glmnet.obj$lambda.min
)

fname.out <- paste(outdir, 'cv.glmnet.obj.rds', sep = '') 
saveRDS(cv.glmnet.obj, file = fname.out)

fname.out <- paste(outdir, 'model.glmnet.rds', sep = '') 
saveRDS(model.glmnet, file = fname.out)
