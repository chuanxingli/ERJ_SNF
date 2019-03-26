# This is the analysis code for run multi-omics analysis in our ERJ paper in a psedo dataset.

fpath = "<Your Path to Folder>"
fpath = "~/Documents/GitHub/ERJ_SNF/" 
rpath = "~/Documents/GitHub/ERJ_SNF/results/"

library(SNFtool);library(pbapply);library(parallel)
source(paste0(fpath,'r/generate_random_dataL.R'))
source(paste0(fpath,'/r/generate_input.R'))
source(paste0(fpath,'r/internal.R'))

#===============================================================================
# 1. Prepare multi-omics psedo dataset

## 1.1 generate a dataset with multi omic, unequal sample and feature size.
### The format in dataL is row for samples and columns for features

raw.equal <- generate_random_dataL(10,3) # equal samples for all omics
raw <- generate_random_dataL(nomics=10,ngroup=4,equal = FALSE,ssize = 30,gnum = seq(100,1000,by=100))
### unequal samples for all omics
sapply(raw$dataL,dim)

## 1.2 prepare formated data for SNF
input <- generate_input(raw$dataL,raw$label,omics=names(raw$dataL),minS=10)
#===============================================================================

# 2. SNF similarity matrix

### global parameters for  SNF
vk=3:10;va=seq(0.3,0.8,by=0.1);vt=c(20,30);
vk=5:10;va=seq(0.3,0.8,by=0.2);vt=c(20,30);

### weight matrix for each data block for each K and alpha combinations
Wi <- weightMatrix(input$dist,vk,va)

### Fused similarity matrix for all possible multi-omic combinations
W=SNFweight(input,Wi,vt)

#===============================================================================
# 3. supervised clustering
## get predictions
label=input$label

# Calculate the number of cores
no_cores <- detectCores() - 3

# Initiate cluster
# !not run long time
t1 = Sys.time()
cl <- makeCluster(no_cores)
clusterExport(cl, "W")
clusterExport(cl, "label")
clusterExport(cl, "fpath")
resc <- parLapply(cl,1:length(W$value),function(a){
        
        library(SNFtool);library(pbapply)
        source(paste0(fpath,'r/internal.R'))
        res=SNF_loocv(W$value[[a]],label,ssize=5,permn=100);return(res)})
        # pernm could be increased to get a stable result
time=Sys.time()
stopCluster(cl)
t2=Sys.time()

## formated results
# nmi(Normalized Mutual Information) = accuracy

unlist(resc) %>% cbind(W$par,nmi=.) %>% merge(.,input$pindex) -> res.nmi

write.csv2(res.nmi,row.names=F,file = paste0(rpath,"res_supervised.csv"))
res.nmi <- read.csv2(paste0(rpath,"res_supervised.csv"))
#===============================================================================

# 4. unsupervised clustering with fix number of groups
C = 4 # fix number of groups
label = input$label
res.nmi = vector("numeric",length=length(W$value))
no_cores <- detectCores() - 3
# Initiate cluster

cl <- makeCluster(no_cores)
clusterExport(cl, "W")
clusterExport(cl, "label")
clusterExport(cl, "C")
resc <- parLapply(cl,1:length(W$value),function(tempi){
        library(SNFtool);library(pbapply)
        group = spectralClustering(W$value[[tempi]], C);
        names(group)=rownames(W$value[[tempi]])
        return(calNMI(group, label[names(group)]))
})
stopCluster(cl)

unlist(resc) %>% cbind(W$par,nmi=.) %>% merge(.,input$pindex) -> res.nmi

write.csv2(res.nmi,row.names=F,file = paste0(rpath,"res_unsupervised_fix.csv"))
res.nmi <- read.csv2(paste0(rpath,"res_unsupervised_fix.csv"))
#===============================================================================
# 5. unsupervised clustering with flexible number of groups

label = input$label
res.nmi = vector("numeric",length=length(W$value))
no_cores <- detectCores() - 3
# Initiate cluster

cl <- makeCluster(no_cores)
clusterExport(cl, "W")
clusterExport(cl, "label")
resc <- parLapply(cl,1:length(W$value),function(tempi){
        library(SNFtool);library(pbapply)
        tempc <- estimateNumberOfClustersGivenGraph(W = W$value[[tempi]],3:15)
        group = spectralClustering(W$value[[tempi]], tempc[[1]]);
        names(group)=rownames(W$value[[tempi]])
        return(calNMI(group, label[names(group)]))
})
stopCluster(cl)

unlist(resc) %>% cbind(W$par,nmi=.) %>% merge(.,input$pindex) -> res.nmi

write.csv2(res.nmi,row.names=F,file = paste0(rpath,"res_unsupervised.csv"))
res.nmi <- read.csv2(paste0(rpath,"res_unsupervised.csv"))