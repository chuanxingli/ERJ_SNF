## Description: Prepare structed input data for SNF analysis by omics name and minimal sample size
## Usage: input <- generate_input(dataL=raw$dataL,label=raw$truelabel,omics=names(raw$dataL),minS=5)
# ==============================================================================
# Input:
# dataL - list of omics data, row for samples, column for features
# label - factor vector of sample groups
# omics - character vector of omics' names
# minS - minimize sample numbers in subgroups
# ==============================================================================
# Output: a list vector with 5 components
# dataL - list of omics data, row for samples, column for features
# label - factor vector of sample groups
# index - list of omics index for each multi-omics combination 
# sampleL - list of sample IDs shared by multi-omics in index
# pindex - data.frame of summary of each multi-omics combinations with "omics" (IDs in index and sampleL),"Min_sample","N_sample","N_omics"  
# dist - list of dist matrix corresponding to dataL


generate_input <- function(dataL,label,omics=names(dataL),minS=5){

        dataL <- dataL[omics]

        
        ## samples
        sample <- names(label)
        dataL <- lapply(dataL,function(a){a[intersect(rownames(a),sample),]})
        
       
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
        for (i in 1:length(dataL)) {
                view = standardNormalization(dataL[[i]])
                # combine train and test data
                Dist[[i]] = dist2(view, view)}

        # library(e1071)
        # Dist[[1]] <- hamming.distance(dataL$Blood_HLA)/nrow(dataL$Blood_HLA)
        # names(Dist) <- names(res$dataL)
        res$dist <- Dist
        
        
        return(res)
        
}