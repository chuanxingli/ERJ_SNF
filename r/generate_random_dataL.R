generate_random_dataL <- function(nomics,ngroup,equal = TRUE,ssize = 30,gnum = seq(100,1000,by=100)){
        ## generate a dataset with multi omic, equal or unequal sample and feature size.
        #[Usage]: raw <- generate_random_dataL(nomics=10,ngroup=4,equal = FALSE,ssize = 30,gnum = seq(100,1000,by=100))
        #[Usage]: raw.equal <- generate_random_dataL(nomics=10,ngroup=3) # equal samples for all omics
        # =========================================================================
        # Inputs:
        # nomics - number of omics we want to generate
        # ngroup - number of sample groups
        # equal - 'TRUE' as default, with the same samples in all omics. otherwise 'FALSE'. If equal = FALSE, possible numbe of samples as an integer vector seq(ssize,ngroup*sszie,by=ssize)
        # ssize - max number of samples in each groups
        # gnum - number of genes/features. An integer vector of possible gene numbers. If the length of gnum == 1, all the omics with the same feature numbers.
        # =========================================================================
        # Outputs: a list vector with two components dataL and label
        # dataL - list of omics data, row for samples, column for features
        # label - factor vector of sample groups
        
        sn = ngroup*ssize # max sample number
        truelabel <- rep(c(1:ngroup),ssize)[sample(sn)]
        names(truelabel) <- paste0("s",1:sn)
        
        set.seed(123)

        dataL <- lapply(1:nomics,function(a){
                if (equal==TRUE){
                        temps <- c(1:sn)
                } else {
                        steps = seq(ssize,sn,by=10)
                        temps <- sample(sn)[1:sample(steps)[1]]}
                if (length(gnum)==1){tempg = gnum
                } else{
                        tempg <- sample(gnum)[1]  
                }
                matrix(rnorm(tempg*length(temps)),length(temps),tempg,dimnames=list(names(truelabel)[temps],NULL))
        })
        names(dataL) <- paste0("omics",1:length(dataL))
        return(res = list(dataL=dataL,label=truelabel))
}