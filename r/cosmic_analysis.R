# This is the analysis code for run multi-omics analysis in our ERJ paper in a psedo dataset.

fpath = "<Your Path to Folder>"
fpath = "~/Documents/GitHub/ERJ_SNF/" 
rpath = "~/Documents/GitHub/ERJ_SNF/results/"
load("/Users/cli/Documents/proj/copd/data/cosmic.RData")
library(SNFtool);library(pbapply);library(parallel)

source(paste0(fpath,'/r/generate_input.R'))
source(paste0(fpath,'r/internal.R'))

#===============================================================================
# 1. Prepare multi-omics karolinska cosmic dataset
# load processed cosmic.RData

# source(paste0(fpath,"r/cosmic_process.R"))
# if (file.exists(paste0(fpath,"data/cosmic.RData"))==FALSE) {
#         
#         cosmic<-cosmic_process(par,del_sample)
#         fdate <- Sys.Date()
#         cosmic$BAL_P_iTRAQ_ONE_impute<-iTRAQ_impute(cosmic$BAL_P_iTRAQ_ONE,10,"./BAL_P_iTRAQ_ONE") # imputation of iTRAQ data
#         setwd(fpath)
#         cosmic$BAL_P_iTRAQ_ALL_impute<-iTRAQ_impute(cosmic$BAL_P_iTRAQ_ONE,10,"./BAL_P_iTRAQ_ALL") # imputation of iTRAQ data
#         setwd(fpath)
#         save(cosmic,fdate,file=paste0(fpath,"data/cosmic.RData"))
# }
load(paste0(fpath,"data/cosmic.RData"))
### The format in cosmic is column for samples and rows for features

## prepare label
group4 <- cosmic$datasummary$cgroup; names(group4) <- rownames(cosmic$datasummary); # label for four groups
group3 <- droplevels(subset(group4,group4!="EC",drop=TRUE))
# group3 <- droplevels(subset(group3,group3!="NH",drop=TRUE))
label <- as.factor(as.numeric(group3)); names(label) <- names(group3) # SH + SC
gs <- intersect(rownames(cosmic$datasummary)[cosmic$datasummary$gender=="female"],names(label)) # females
label <- label[gs]

rm(group3,group4,gs)

omics <- c("BAL_T_mRNA","BAL_T_miR","BAL_P_DIGE","BAL_P_iTRAQ_ALL_impute","BEC_P_TMT_impute","EXO_T_miR","Serum_M_Non_targeted","Serum_M_Oxylip","BALF_M_Oxylip")
#input <- datapre(cosmic,label,minS=5,mRNArsd=5.5)
# minS - minimize sample numbers in subgroups = 5

dataL <- lapply(cosmic[omics],t)
## 1.2 prepare formated data for SNF
input <- generate_input(dataL,label,omics=omics,minS=5)

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
C = 3 # fix number of groups
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