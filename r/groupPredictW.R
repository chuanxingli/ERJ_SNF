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