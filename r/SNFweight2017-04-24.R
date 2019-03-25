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



