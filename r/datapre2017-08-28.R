## Description: Dreparation of input data for SNF analysis
## Usage: input <- datapre(cosmic,label,minS,mRNArsd)
# input:
# cosmic - list of multi-omics data from COSMIC project
# label - factor vector of sample groups
# minS - minimize sample numbers in subgroups
# mRNArsd - RSD threshold for mRNAs. If 0, means no filter.
# output:
# dataL - list of omics data, row for samples, column for features
# label - factor vector of sample groups
# index - list of omics index for each multi-omics combination 
# sampleL - list of sample IDs shared by multi-omics in index
# pindex - data.frame of summary of each multi-omics combinations with "omics" (IDs in index and sampleL),"Min_sample","N_sample","N_omics"  
# dist - list of dist matrix corresponding to dataL

## The data processing
# Three miR expression profiles filtered by LOQ >= 5.5
# BAL_T_mRNA expression filtered by LOQ>=5,5 and RSD > mRNArsd
# minimized sample number in subgroup is minS
# Serum metabolic data are combined to Serum_M_targeted



datapre <- function(cosmic,label,minS,mRNArsd){
        
        # omics will be used
        omics <- c("HLA_typing","BAL_T_mRNA","BAL_T_miR","BAL_P_DIGE","BAL_P_iTRAQ_ALL_impute","BEC_T_miR","BEC_P_TMT_impute","EXO_T_miR","Serum_M_Non_targeted","Serum_M_Oxylip","Serum_M_Biocrates","Serum_M_Sphingolipid","Serum_M_Kynurenine","BALF_M_Oxylip")
        
        ## prepare dataL
        dataL <- cosmic[omics]
        # miRNA profiles filter by ave >= 5.5
        teo <- c("BAL_T_miR","BEC_T_miR","EXO_T_miR")
        for (tei in 1:length(teo)){
                eval(parse(text=c("tedata <- dataL$",teo[tei])))
                tedata <- tedata[tedata[,"aveA"]>=5.5,]
                eval(parse(text=paste0("dataL$",teo[tei],"<- tedata")))}
        
        # mRNA profiles filter by ave >= 5.5 & Rsd >= mRNArsd
        data <-dataL$BAL_T_mRNA[rowMeans(dataL$BAL_T_mRNA)>=5.5,]
        t <- sapply(1:nrow(data),function(a){sd(data[a,])/mean(data[a,])})
        data = data[t >= mRNArsd,]
        dataL$BAL_T_mRNA=data
        
        
        dataL <- lapply(dataL,t) # convert to samples in row and features in column
        
        rm(list=ls(pattern = "te"))
        names(dataL) <- omics 
        
        ## samples
        sample <- names(label)
        dataL <- lapply(dataL,function(a){a[intersect(rownames(a),sample),]})
        
        # merge three serum targeted data into one Serum_M_targeted
        sample <- intersect(intersect(rownames(dataL$Serum_M_Kynurenine),rownames(dataL$Serum_M_Biocrates)),rownames(dataL$Serum_M_Sphingolipid))
        dataL$Serum_M_Targeted <- cbind(dataL$Serum_M_Biocrates[sample,],dataL$Serum_M_Sphingolipid[sample,],dataL$Serum_M_Sphingolipid[sample,])
        
        omics <- c("HLA_typing","BAL_T_mRNA","BAL_T_miR","BAL_P_DIGE","BAL_P_iTRAQ_ALL_impute","BEC_T_miR","BEC_P_TMT_impute","EXO_T_miR","Serum_M_Non_targeted","Serum_M_Targeted","Serum_M_Oxylip","BALF_M_Oxylip")
        dataL <- dataL[omics]
        names(dataL) <- c("Blood_HLA","BAL_T_mRNA","BAL_T_miR","BAL_P_DIGE","BAL_P_iTRAQ","BEC_T_miR","BEC_P_TMT","EXO_T_miR","Serum_M_Non_targeted","Serum_M_Targeted","Serum_M_Oxylip","BALF_M_Oxylip")
        
        ## combinations of omics
        nd <- length(dataL)
        index <- vector("list") # index of omics
        for (tei in 1:nd){
                te <- combn(nd,tei);
                te2 <- lapply(1:ncol(te),function(a) {te[,a]})
                index = c(index,te2)
        }
        rm(list=ls(pattern = "te"))
        
        ## summary of samples for each omic combinations
        
        sampleL <- vector('list',length=length(index)) # sample COSMIC IDs for each combination
        for (te in 1:length(index)){
                te1 <- index[[te]]
                te2 <- lapply(te1, function(b){rownames(dataL[[b]])})
                te3 <- te2[[1]]
                for (tei in (1:length(te2))){
                        te3 <- intersect(te3,te2[[tei]])
                }
                sampleL[[te]] <- te3
        }
        sampleL <- lapply(sampleL,function(a){intersect(a,names(label))})
        rm(list=ls(pattern = "te"));rm(nd);
        
        te <- sapply(sampleL,length)
        te_min <- sapply(sampleL,function(a){min(table(label[intersect(names(label),a)]))})
        # number of smallest group for each omics combination and label group

        pindex <- cbind.data.frame(omics=1:length(te_min),Min_sample=te_min)
        pindex <- pindex[pindex[,"Min_sample"]>=minS,]
        index <- index[pindex[,"omics"]]
        sampleL <- sampleL[pindex[,"omics"]]
        # result table of sample/omics summary
        pindex[,"N_sample"] <- te[pindex[,"omics"]]
        pindex[,"N_omics"] <- sapply(index,length)
        
        pindex[,"omics"] <- 1:nrow(pindex)
        rm(list=ls(pattern = "te"));
        res <- list(dataL,sampleL,index,pindex,label)
        names(res) <- c("dataL","sampleL","index","pindex","label")
        
        # data normalization
        Dist <- vector("list", length = length(dataL))
        for (i in 2:length(dataL)) {
                view = standardNormalization(dataL[[i]])
                # combine train and test data
                Dist[[i]] = dist2(view, view)}
        library(e1071)
        Dist[[1]] <- hamming.distance(dataL$Blood_HLA)/nrow(dataL$Blood_HLA)
        names(Dist) <- names(res$dataL)
        res$dist <- Dist

        
        return(res)

}