## The main code for supervised SNF analysis of multi-omics integration in Karolinska COSMIC study
fpath = "<Your Path to Main Folder>"
fpath = "/Users/cli/Documents/GitHub/ERJ_SNF/"; 
setwd(fpath) 
dpath = "./data/cosmic.RData"
rpath ="./results/supervise/"# path of results
if (!dir.exists(rpath)){
        dir.create(rpath)
}

library(SNFtool);library(pbapply);library(parallel)
source("./r/internal.R")
source('./r/datapre2017-08-28.R')
source('./r/weightMatrix2016-11-10.R')
source('./r/SNFweight2017-04-24.R')
source('./r/groupPredictW.R')
source('./r/SNF_loocv2016-12-07.R')



### predict by SNF
vk=3:10;va=seq(0.3,0.8,by=0.1);vt=c(20,30);

Wi <- weightMatrix(input$dist,vk,va)
# weight matrix for each data block for each K and alpha combinations
W=SNFweight(input,Wi,vt)

# save(W,Wi,file=paste0(rpath,"weight.Rdata")) # around 8 mins

# load(paste0(rpath,"weight.Rdata"))
# load(paste0(rpath,"input.Rdata"))

## get predictions
label=input$label

# Calculate the number of cores
no_cores <- detectCores() - 3

# Initiate cluster
# not run longtime
t1 = Sys.time()
cl <- makeCluster(no_cores)
clusterExport(cl, "W")
clusterExport(cl, "label")
        resc <- parLapply(cl,1:length(W$value),function(a){

                library(SNFtool);library(pbapply)
                source('./r/SNF_loocv2016-12-07.R')
                source('./r/groupPredictW.R')
                res=SNF_loocv(W$value[[a]],label,ssize=5,permn=100);return(res)})
time=Sys.time()
stopCluster(cl)
t2=Sys.time()

## formated results
# nmi(Normalized Mutual Information) = accuracy

unlist(resc) %>% cbind(W$par,nmi=.) %>% merge(.,input$pindex) -> res.nmi

save(res.nmi,file=paste0(rpath,"resNmi.Rdata"))