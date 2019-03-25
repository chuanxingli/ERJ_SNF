## The main code for unsupervised SNF analysis of multi-omics integration in Karolinska COSMIC study
fpath = "<Your Path to Main Folder>"
fpath = "/Users/cli/Documents/GitHub/ERJ_SNF/"; 
setwd(fpath) 
dpath = "./data/cosmic.RData"
rpath ="./results/unsupervise/"# path of results
if (!dir.exists(rpath)){
        dir.create(rpath)
}
library(SNFtool);library(pbapply);library(parallel)
source('./r/datapre2017-08-28.R')
source('./r/weightMatrix2016-11-10.R')
source('./r/SNFweight2017-04-24.R')

## input files
load(dpath)

## prepare data for smoker and COPD smoker female
group4 <- cosmic$datasummary$cgroup; names(group4) <- rownames(cosmic$datasummary); # label for four groups
group3 <- droplevels(subset(group4,group4!="EC",drop=TRUE))
# group3 <- droplevels(subset(group3,group3!="NH",drop=TRUE))
label <- as.factor(as.numeric(group3)); names(label) <- names(group3) # SH + SC
gs <- intersect(rownames(cosmic$datasummary)[cosmic$datasummary$gender=="female"],names(label)) # females
label <- label[gs]
# two groups: female SH (1, n = 20), female SC (2, n = 12)
rm(group3,group4,gs)

# omics <- c("HLA_typing","BAL_T_mRNA","BAL_T_miR","BAL_P_DIGE","BAL_P_iTRAQ_ALL_impute","BEC_T_miR","BEC_P_TMT_impute","EXO_T_miR","Serum_M_Non_targeted","Serum_M_Oxylip","Serum_M_Biocrates","Serum_M_Kynurenine","Serum_M_Sphingolipid","BALF_M_Oxylip")
input <- datapre(cosmic,label,minS=5,mRNArsd=5.5)
# minS - minimize sample numbers in subgroups = 5
# mRNArsd - RSD threshold for mRNAs. If 0, means no filter = 5.5.


rm(cosmic, label, datapre)


vk=3:10;va=seq(0.3,0.8,by=0.1);vt=c(20,30);

Wi <- weightMatrix(input$dist,vk,va)
W=SNFweight(input,Wi,vt)

save(input,file=paste0(rpath,"input.rdata"))
save(Wi,W,file=paste0(rpath,"weight.rdata"))

C = 3
label = input$label
res.nmi = vector("numeric",length=length(W$value))
no_cores <- detectCores() - 3
# Initiate cluster
t1 = Sys.time()
cl <- makeCluster(no_cores)
clusterExport(cl, "W")
clusterExport(cl, "label")
resc <- parLapply(cl,1:length(W$value),function(tempi){
        library(SNFtool);library(pbapply)
        group = spectralClustering(W$value[[tempi]], 3);
        names(group)=rownames(W$value[[tempi]])
        return(calNMI(group, label[names(group)]))
})
stopCluster(cl)

res.nmi <- unlist(resc)
res.nmi <- cbind.data.frame(W$par,nmi=res.nmi)
res.nmi = cbind(res.nmi,input$pindex[match(res.nmi[,"omics"],input$pindex[,"omics"]),])
res.nmi = res.nmi[,-6] # delete duplicate omics column
save(res.nmi,file=paste0(rpath,"resnmi.rdata"))
