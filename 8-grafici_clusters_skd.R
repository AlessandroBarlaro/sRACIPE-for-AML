#grafici

cluster_names <- c('CL1', 'CL2')
#names(cluster_names) <- c('AML','Untreated')
names(cluster_names) <- c('Untreated', 'AML')

datadir <- './data/'
figdir <- './figs.props/'
dir.create(figdir)

# Load cluster proportions
#-----------------------
fname.data <- paste(datadir, 'cluster_props.csv', sep = '') 
cluster_props <- read.csv(file = fname.data, row.names = 1)

# Take percent 
cluster_props <- cluster_props/100

# Sort the data by AML cluster
cluster_props.sorted <- cluster_props[order(cluster_props$Cluster_1, decreasing = F), ]

# Change column names to phenotypes
#----------------------------------
colnames(cluster_props) <- names(cluster_names)
colnames(cluster_props.sorted) <- names(cluster_names)

# Plot using ggplot
#================== 
# melt data
df <- as.data.frame(cbind(rownames(cluster_props.sorted),
                          cluster_props.sorted))

colnames(df) <- c('TF', colnames(cluster_props.sorted))

library(reshape2)
df2 <- melt(df, id=c('TF')) 
colnames(df2) <- c('TF', 'Cluster', 'prop')

df2$TF <- factor(df2$TF, levels = as.character(df$TF))

library(ggplot2)
p1 <- ggplot(df2, aes(x=TF, y=prop, fill=Cluster)) + 
  geom_bar(stat="identity")  +  
  coord_flip() +
  #scale_fill_manual(values=c("red", "black", "dodgerblue", "orange")) + 
  #scale_fill_manual(values=c("red", "dodgerblue")) + 
  scale_fill_manual(values=c("dodgerblue", "red")) + 
  labs(x='Transcription Factor', y='Percent Models') 

p1

WIDTH <- 4  
HEIGHT <- 8 
figname <- paste(figdir, 'cluster_props-', WIDTH, 'x', HEIGHT, '.pdf', sep = '')
ggsave(figname, device = 'pdf', width = WIDTH, height = HEIGHT )

WIDTH <-  6  
HEIGHT <- 8  

figname <- paste(figdir, 'cluster_props-', WIDTH, 'x', HEIGHT, '.pdf', sep = '')
ggsave(figname, device = 'pdf', width = WIDTH, height = HEIGHT )



