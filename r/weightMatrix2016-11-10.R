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