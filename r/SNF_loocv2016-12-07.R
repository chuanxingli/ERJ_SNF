SNF_loocv <- function(snfw,label,ssize,permn)
{
        
        ## leave one out cross-validation with equal sample size
        ## output is the NMI 
        
        ## leave one out cross-validation with equal sample size
        # snfw - list of weight matrix fused
        # label - factor vector of sample true label 
        # ssize - sample size inside of subgroup, could be larger than the real minimal sample size, then the output will be NA
        # number of permutation time
        set.seed(123)
        sample = rownames(snfw)
        c = min(table(label[sample]))
        if (c > ssize) {
        res = matrix(NA,ncol=permn,nrow=2)
        # get out one random sample from dataset
        res[1:2,1:permn] <- sapply(1:permn,function(a){
                test = sample(1:length(sample))[1]
        # construct balanced training data
        trainSample = sample[-test] # raw train sample set
        
        l=label[trainSample] # label of train sample

        gs = lapply(1:length(unique(l)),function(a){names(l[l==unique(l)[a]])})
        # seperate samples by groups
        
        us = lapply(1:length(gs),function(a){sample(gs[[a]],ssize)})
        
        trainSample = unlist(us)
        
        # prepare input W and groups for groupPredictW.R
        sp_sample=c(trainSample,sample[test])
        newLabel = groupPredictW(snfw[sp_sample,sp_sample],as.numeric(label[trainSample]))
        return(c(test,newLabel[sample[test]]))
        })
        rest = res[2,]
        names(rest) = sample[res[1,]]
        v = calNMI(label[names(rest)],rest)
        return(v)
        }
        else {return(NA)}
}