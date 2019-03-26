weightMatrix <- function(Dist,vk,va) 
{
        ## function to generate weight matrix for each data distance matrix in list dist by parameter K in vk, and alpha in va
        # vk=3:6;va=seq(0.3,0.8,by=0.1);Dist=input$dist;
        par <- cbind(k=rep(vk,length(va)),a=sort(rep(va,length(vk))))
        res <- lapply(c(1:nrow(par)),function(b){
                # b is the rownumber of par
                K = par[b,"k"];alpha = par[b,"a"];
                Wi = lapply(Dist,function(a){affinityMatrix(a, K, alpha)})
                names(Wi) = names(Dist);
                return(Wi)
        })
        result <- list(value=res,par=par)       
        return(result)
}

SNFweight <- function(input=input,Wi=Wi,vt=vt){
        ## This function is to generate SNF merged W weight matrix for each omics combinations
        ## for single omic, just use the original weight matrix
        input = input
        index =  input$index
        sampleL = input$sampleL
        par <- cbind(Wi$par %x% rep(1, length(vt)),t=rep(vt,nrow(Wi$par))) # combination between vt and Wi$par(k+a)
        par <- cbind(rep(1:length(index),nrow(par)),par %x% rep(1,length(index))) # combination between omics and par
        colnames(par)=c("omics","k","a","t")
        n = nrow(par)
        # code for non parallel computing 
        # value <- pblapply(1:n,function(a){
        #         wm <- Wi$value[[which(Wi$par[,"k"]==par[a,"k"]&Wi$par[,"a"]==par[a,"a"])]]
        #         # weight matrix for matched k and alpha
        #         wm <- wm[index[[par[a,"omics"]]]]
        #         wm <- lapply(wm, function(b){b[sampleL[[par[a,"omics"]]],sampleL[[par[a,"omics"]]]]})
        #         W = SNF(wm, par[a,"k"], par[a,"t"])
        #         rownames(W) <- colnames(W) <- sampleL[[par[a,"omics"]]]
        #         return(W)
        #         
        # })
        
        
        # code for parallel computing
        # Calculate the number of cores
        no_cores <- detectCores() - 3
        
        # Initiate cluster
        cl <- makeCluster(no_cores)
        clusterExport(cl, "Wi")
        clusterExport(cl, "input") 
        clusterExport(cl, "par")
        value <- parLapply(cl,1:n,function(a){
                .parSNF <- function(a,uWi=Wi,uindex=input$index,usampleL=input$sampleL,upar=par){
                        wm <- uWi$value[[which(uWi$par[,"k"]==upar[a,"k"]&uWi$par[,"a"]==upar[a,"a"])]]
                        # weight matrix for matched k and alpha
                        wm <- wm[uindex[[upar[a,"omics"]]]]
                        # weight matrix for this combination
                        wm <- lapply(wm, function(b){b[usampleL[[upar[a,"omics"]]],usampleL[[upar[a,"omics"]]]]})
                        # get sample set
                        if (length(wm)>1) {W = SNF(wm, upar[a,"k"], upar[a,"t"])} else {W=wm[[1]]}
                        # SNF merge by k,t,a and this omics
                        rownames(W) <- colnames(W) <- usampleL[[upar[a,"omics"]]]
                        return(W)
                }
                library(SNFtool)
                W=.parSNF(a);return(W)})
        
        stopCluster(cl)
        
        res <- list(value=value,par=par)
        return(res)
}



groupPredictW <- function (W,groups, method = 1) 
{
        ### This is a small edit of groupPredict.r in SNFtools
        ## function based on groupPredict.r to only input SNF weight matrix W, number of samples(train + test), and label of train samples
        # Wi = vector("list", length = length(train))
        # for (i in 1:length(train)) {
        #         view = standardNormalization(rbind(train[[i]], test[[i]]))
        #         # combine train and test data
        #         Dist1 = dist2(view, view)
        #         # distance matrix including both train and test data
        #         Wi[[i]] = affinityMatrix(Dist1, K, alpha)
        #         # weight matrix including both train and test data
        # }
        # W = SNF(Wi, K, t) # final weight matrix from SNF, also including both train and test data
        nsample = nrow(W)
        Y0 = matrix(0, nsample, max(groups)) # initial label for all samples
        for (i in 1:length(groups)) Y0[i, groups[i]] = 1
        Y = .csPrediction(W, Y0, method)
        newgroups = rep(0, nsample)
        for (i in 1:nrow(Y)) newgroups[i] = which(Y[i, ] == max(Y[i,]))
        names(newgroups)=rownames(W)
        return(newgroups)
}

.csPrediction <- function(W,Y0,method){
        ###This function implements the label propagation to predict the label(subtype) for new patients.	
        ### note method is an indicator of which semi-supervised method to use
        # method == 0 indicates to use the local and global consistency method
        # method >0 indicates to use label propagation method.
        
        alpha=0.9;
        P= W/rowSums(W)
        if(method==0){
                Y= (1-alpha)* solve( diag(dim(P)[1])- alpha*P)%*%Y0;
        } else {
                NLabel=which(rowSums(Y0)==0)[1]-1;
                Y=Y0;
                for (i in 1:1000){
                        Y=P%*%Y;
                        Y[1:NLabel,]=Y0[1:NLabel,];
                }
        }
        return(Y);
}

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