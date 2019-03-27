cosmic_input <- function(par,del_sample) {
  library(gdata)
  library(limma) # Quantile normalization of data (normalizeQuantiles )
  #rn <- c("datasummary","clinic","BAL_P_DIGE","BAL_P_iTRAQ","BAL_T_miR","BAL_T_mRNA","BEC_T_miR","EXO_T_miR","Metabolic","HLA_typing")
  cosmic <- vector(mode = 'list')
  # names(cosmic) <- rn
  # rm(rn)
  
  # # import datasummary
  datasummary_raw<-read.xls(par["datasummary","fname"], sheet=par["datasummary","sheet"], nrows = par["datasummary","nrows"], na.strings = "NA")
  datasummary<-datasummary_format(datasummary_raw, par["datasummary","scol"])
  datasummary <-datasummary[is.element(datasummary$COSMIC_ID, del_sample)==FALSE,] # delet sample ID in del_sample
  rm(datasummary_raw)
  cosmic$datasummary <- datasummary
  
  # import clinical data
  subject <- "clinic"
  data<-read.xls(par[subject,"fname"], sheet=par[subject,"sheet"], na.strings = "NA")
  data[,1]<-as.character(data[,1])
  # t <- paste(as.character(data[duplicated(data[,1])==TRUE,1]),'2',sep="")
  # data[duplicated(data[,1])==TRUE,1]<-t
  rownames(data) <- data[,1]
  c <- as.integer(as.matrix(data["COSMIC_ID",]))
  colnames(data)[is.na(c)==0]<-datasummary$barcode[match(c,datasummary$COSMIC_ID)][is.na(c)==0]
  data <- data[,is.na(colnames(data))==FALSE]
  data <- data[,2:122]
  
  HLA_typing <- data[which(data$TISSUE=="HLA"),]
  col_to_keep <- !is.na(match(colnames(HLA_typing),rownames(datasummary)))
  HLA_typing <- HLA_typing[,which(col_to_keep)]
  HLA_typing <- as.matrix(HLA_typing)
  class(HLA_typing) <- "integer"
  
  cosmic$HLA_typing <- HLA_typing
  
  cosmic$clinic <- data
  rm(data,c,subject,HLA_typing,col_to_keep)
  
  # import BAL_T_mRNA
  subject <- "BAL_T_mRNA"
  data<-read.csv(par[subject,"fname"],sep=",",header=TRUE,nrows=par[subject,"nrows"])
  c <- as.integer(as.matrix(data[1,]))
  colnames(data)[is.na(c)==0]<-cosmic$datasummary$barcode[match(c,cosmic$datasummary$COSMIC_ID)][is.na(c)==0]
  data <- as.matrix(data)
  rownames(data) <- data[,1]
  data <- data[2:dim(data)[1],2:dim(data)[2]]
  geneid <- data[,is.element(colnames(data),cosmic$datasummary$barcode)==FALSE]
  data <- data[,is.element(colnames(data),cosmic$datasummary$barcode)==TRUE]
  class(data) <- "numeric"
  eval(parse(text = paste("cosmic$",subject,"<-data",sep="")))
  eval(parse(text = paste("cosmic$",subject,"_geneid<-geneid",sep="")))
  rm(data,subject,c,geneid)
  
  # import BAL_P_DIGE
  subject <- "BAL_P_DIGE"
  data<-read.xls(par[subject,"fname"], sheet=par[subject,"sheet"], na.strings = "NA")
  c <- as.integer(as.matrix(data[1,]))
  colnames(data)[is.na(c)==0]<-cosmic$datasummary$barcode[match(c,cosmic$datasummary$COSMIC_ID)][is.na(c)==0]
  data <- as.matrix(data)
  rownames(data) <- data[,1]
  data <- data[2:dim(data)[1],2:dim(data)[2]]
  class(data) <- "numeric"
  data <- log2(10^data)
  eval(parse(text = paste("cosmic$",subject,"<-data",sep="")))
  rm(data,subject,c)
  
  # import BAL_P_iTRAQ
  subject <- "BAL_P_iTRAQ_ONE"
  data<-read.xls(par[subject,"fname"], sheet=par[subject,"sheet"], na.strings = "NA",head=FALSE)
  rownames(data) <- data[,1]
  data <- data[,-1]
  c <- as.integer(as.matrix(data[1,]))
  colnames(data)[is.na(c)==0]<-cosmic$datasummary$barcode[match(c,cosmic$datasummary$COSMIC_ID)][is.na(c)==0]
  data <- data[-1,]
  data <- as.matrix(data)
  class(data) <- "numeric"
  eval(parse(text = paste("cosmic$",subject,"<-data",sep="")))
  rm(data,subject,c)
  
  subject <- "BAL_P_iTRAQ_ALL"
  data<-read.xls(par[subject,"fname"], sheet=par[subject,"sheet"], na.strings = "NA",head=FALSE)
  rownames(data) <- data[,1]
  data <- data[,-1]
  c <- as.integer(as.matrix(data[1,]))
  colnames(data)[is.na(c)==0]<-cosmic$datasummary$barcode[match(c,cosmic$datasummary$COSMIC_ID)][is.na(c)==0]
  data <- data[-1,]
  data <- as.matrix(data)
  class(data) <- "numeric"
  eval(parse(text = paste("cosmic$",subject,"<-data",sep="")))
  rm(data,subject,c)
  # run iTRAQ_impute.R to determin KNN imputation K
  # source('~/Documents/proj/copd/data/iTRAQ_impute.R')
  # cosmic$BAL_P_iTRAQ_ONE_impute <- iTRAQ_impute(cosmic$BAL_P_iTRAQ_ONE,10,"./BAL_P_iTRAQ_ONE")
  # cosmic$BAL_P_iTRAQ_ALL_impute <- iTRAQ_impute(cosmic$BAL_P_iTRAQ_ONE,10,"./BAL_P_iTRAQ_ALL")
  
  # import BAL_T_miR
  subject <- "BAL_T_miR"
  data<-read.xls(par[subject,"fname"], sheet=par[subject,"sheet"], na.strings = "NA")
  rownames(data) <- data[,"Probe.Sequence"]
  data <- data[,c(colnames(data)[par[subject,"scol"]:dim(data)[2]],"aveA")]
  c <- as.integer(as.matrix(data[1,]))
  colnames(data)[is.na(c)==0]<-cosmic$datasummary$barcode[match(c,cosmic$datasummary$COSMIC_ID)][is.na(c)==0]
  colnames(data)[dim(data)[2]]<-"aveA"
  data <- as.matrix(data)
  data <- data[2:dim(data)[1],]
  class(data) <- "numeric"
  data <- data[,is.na(colnames(data))==FALSE]
  eval(parse(text = paste("cosmic$",subject,"<-data",sep="")))
  rm(data,subject,c)
  
  # import BEC_T_miR
  subject <- "BEC_T_miR"
  data<-read.xls(par[subject,"fname"], sheet=par[subject,"sheet"], na.strings = "NA",head=FALSE)
  rownames(data) <- data[,1]
  data <- data[,-1]
  c <- as.integer(as.matrix(data[1,]))
  colnames(data)[is.na(c)==0]<-cosmic$datasummary$barcode[match(c,cosmic$datasummary$COSMIC_ID)][is.na(c)==0]
  colnames(data)[1]<-"aveA"
  data <- data[,is.na(colnames(data))==FALSE]
  data <- data[-1,]
  data <- as.matrix(data)
  class(data) <- "numeric"
  eval(parse(text = paste("cosmic$",subject,"<-data",sep="")))
  rm(data,subject,c)
  
  # import EXO_T_miR
  subject <- "EXO_T_miR"
  data<-read.xls(par[subject,"fname"], sheet=par[subject,"sheet"], na.strings = "NA")
  rownames(data) <- data[,"Probe.Sequence"]
  data <- data[,c(colnames(data)[par[subject,"scol"]:dim(data)[2]],"aveA")]
  c <- as.integer(as.matrix(data[1,]))
  colnames(data)[is.na(c)==0]<-cosmic$datasummary$barcode[match(c,cosmic$datasummary$COSMIC_ID)][is.na(c)==0]
  colnames(data)[dim(data)[2]]<-"aveA"
  data <- as.matrix(data)
  data <- data[2:dim(data)[1],]
  class(data) <- "numeric"
  data <- data[,is.na(colnames(data))==FALSE]
  eval(parse(text = paste("cosmic$",subject,"<-data",sep="")))
  rm(data,subject,c)
  
  # # import Metabolic
  subject <- "BALF_M_Oxylip"
  data<-read.xls(par[subject,"fname"], sheet=par[subject,"sheet"], na.strings = "NA",head=FALSE)
  rownames(data)<-data[,2] # rowname as SecID
  data <- data[,-4:-1] # delete first 4 columns
  c <- as.integer(as.matrix(data[1,]))
  colnames(data)[is.na(c)==0]<-cosmic$datasummary$barcode[match(c,cosmic$datasummary$COSMIC_ID)][is.na(c)==0]
  data <- data[-1,] # delete first row
  data <- data[,is.element(colnames(data),cosmic$data$barcode)==TRUE] # delete sample without barcode
  data <- data[,colSums(is.na(data))!=nrow(data)]
  data <- as.matrix(data)
  tedata <-t(data)
  te <- which(is.na(tedata),arr.ind=TRUE)
  tedata[te[1],te[2]]<-min(tedata[is.na(tedata[,8])==FALSE,8])/3
  data <- t(tedata);
  class(data) <- "numeric";
  eval(parse(text = paste("cosmic$",subject,"<-data",sep="")))
  rm(data,subject,c)
  
  subject <- "Serum_M_Oxylip"
  data<-read.xls(par[subject,"fname"], sheet=par[subject,"sheet"], na.strings = "NA",head=FALSE)
  rownames(data)<-data[,2] # rowname as SecID
  data <- data[,-4:-1] # delete first 4 columns
  c <- as.integer(as.matrix(data[1,]))
  colnames(data)[is.na(c)==0]<-cosmic$datasummary$barcode[match(c,cosmic$datasummary$COSMIC_ID)][is.na(c)==0]
  data <- data[-1,] # delete first row
  data <- data[,is.element(colnames(data),cosmic$data$barcode)==TRUE] # delete sample without barcode
  data <- data[,colSums(is.na(data))!=nrow(data)]
  data <- as.matrix(data)
  class(data) <- "numeric";
  eval(parse(text = paste("cosmic$",subject,"<-data",sep="")))
  rm(data,subject,c)
  
  subject <- "Oxylip_annot"
  data<-read.xls(par[subject,"fname"], sheet=par[subject,"sheet"], na.strings = "NA",head=TRUE)
  rownames(data)<-data[,1] # rowname as SecID
  eval(parse(text = paste("cosmic$",subject,"<-data",sep="")))
  rm(data,subject)
  
  subject <- "Serum_M_Non_targeted"
  data<-read.xls(par[subject,"fname"], sheet=par[subject,"sheet"], na.strings = "NA",head=FALSE)
  rownames(data)<-data[,1] # rowname as KIID
  data <- data[,-8:-1] # delete first 8 columns
  c <- as.integer(as.matrix(data[1,]))
  colnames(data)[is.na(c)==0]<-cosmic$datasummary$barcode[match(c,cosmic$datasummary$COSMIC_ID)][is.na(c)==0]
  data <- data[-1,] # delete first row
  data <- data[,is.element(colnames(data),cosmic$data$barcode)==TRUE] # delete sample without barcode
  data <- as.matrix(data)
  data <- log2(data)
  eval(parse(text = paste("cosmic$",subject," <- data",sep="")))
  rm(data,subject,c)
  
  subject <- "Non_targeted_annot"
  data<-read.xls(par[subject,"fname"], sheet=par[subject,"sheet"], na.strings = "NA",head=TRUE)
  rownames(data)<-data[,1] # rowname as KIID
  data <- data[,1:8]
  eval(parse(text = paste("cosmic$",subject," <- data",sep="")))
  rm(data,subject)
  # update annotation
  load("./rawdata/conversionTableHMDB2Kegg.RData")
  cosmic$Non_targeted_annot <- merge(cosmic$Non_targeted_annot,conversionTable,by.x="HMDBID",by.y="HMDB",all.x=TRUE)
  
  
  subject <- "Serum_M_Biocrates"
  data<-read.xls(par[subject,"fname"], sheet=par[subject,"sheet"], na.strings = "NA",head=TRUE)
  data <- t(data) #Transpose
  c <- as.integer(as.matrix(data[1,]))
  colnames(data)[is.na(c)==0]<-cosmic$datasummary$barcode[match(c,cosmic$datasummary$COSMIC_ID)][is.na(c)==0]
  data <- data[-2:-1,]
  class(data) <-"numeric"
  eval(parse(text = paste("cosmic$",subject," <- data",sep="")))
  rm(data,subject,c)
  # delete missing values in cosmic$Serum_M_Biocrates
  cosmic$Serum_M_Biocrates <- cosmic$Serum_M_Biocrates[, which(is.na(colnames(cosmic$Serum_M_Biocrates))==FALSE)] # delete columns with name of NA
  cosmic$Serum_M_Biocrates <- cosmic$Serum_M_Biocrates[which(rowSums(is.na(cosmic$Serum_M_Biocrates))==0),] # delete rows with missing values
  
  subject <- "Serum_M_Kynurenine"
  data<-read.xls(par[subject,"fname"], sheet=par[subject,"sheet"], na.strings = "NA",head=TRUE)
  data <- t(data) #Transpose
  c <- as.integer(as.matrix(data[1,]))
  colnames(data)[is.na(c)==0]<-cosmic$datasummary$barcode[match(c,cosmic$datasummary$COSMIC_ID)][is.na(c)==0]
  data <- data[-2:-1,]
  
  # imputation of missing value as 1/3LLOQ
  tedata <-t(data)
  te <- which(is.na(tedata),arr.ind=TRUE)
  tem <- min(tedata[is.na(tedata[,2])==FALSE,2])/3
  tedata[te[1,1],te[1,2]]<- tem
  tedata[te[2,1],te[2,2]]<-tem
  data <- t(tedata)
  
  class(data) <-"numeric";
  eval(parse(text = paste("cosmic$",subject," <- data",sep="")))
  rm(data,subject,c)
  
  subject <- "Serum_M_Sphingolipid"
  data<-read.xls(par[subject,"fname"], sheet=par[subject,"sheet"], na.strings = "NA",head=TRUE)
  data <- t(data) #Transpose
  c <- as.integer(as.matrix(data[1,]))
  colnames(data)[is.na(c)==0]<-cosmic$datasummary$barcode[match(c,cosmic$datasummary$COSMIC_ID)][is.na(c)==0]
  data <- data[-2:-1,]
  eval(parse(text = paste("cosmic$",subject," <- data",sep="")))
  rm(data,subject,c)
  
  
  # import miRNA_annot_v3.6
  subject <- "miRNA_annot_v3.6"
  data<-read.xls(par[subject,"fname"],header=TRUE)
  data <- data[is.na(data$Symbol.hsa)==FALSE,]
  rownames(data) <- data[,"Probe.Sequence"]
  eval(parse(text = paste("cosmic$",subject,"<-data",sep="")))
  rm(data,subject)
  
  # re-annotate miRNA symbol
  subject <- c("BAL_T_miR","BEC_T_miR","EXO_T_miR")
  for (i in 1:(length(subject))){
    eval(parse(text = paste("data <- cosmic$",subject[i],sep="")))
    c <- is.element(rownames(data),rownames(cosmic$miRNA_annot_v3.6))
    data <- data[c==TRUE,]
    rownames(data) <- cosmic$miRNA_annot_v3.6[rownames(data),"Symbol.hsa"]
    eval(parse(text = paste("cosmic$",subject[i],"<-data",sep="")))
  }
  rm(subject,data)
  
  # BEC_P_TMT
  subject <- "BEC_P_TMT"
  data<-read.xls(par[subject,"fname"], sheet=par[subject,"sheet"], na.strings = "NA",head=FALSE)
  rownames(data) <- data[,1]
  data <- data[,-1]
  c <- as.integer(as.matrix(data[1,]))
  colnames(data)[is.na(c)==0]<-cosmic$datasummary$barcode[match(c,cosmic$datasummary$COSMIC_ID)][is.na(c)==0]
  data <- data[-1,]
  data <- as.matrix(data)
  class(data) <- "numeric"
  eval(parse(text = paste("cosmic$",subject,"<-data",sep="")))
  rm(data,subject,c)
  # run iTRAQ_impute.R to determin KNN imputation K
  # source('~/Documents/proj/copd/data/iTRAQ_impute.R')
  # cosmic$BEC_P_TMT_impute <- iTRAQ_impute(cosmic$BEC_P_TMT,10,"./BEC_P_TMT")

  
  # merge three serum targeted data into one Serum_M_targeted
  sample <- intersect(intersect(rownames(cosmic$Serum_M_Kynurenine),rownames(cosmic$Serum_M_Biocrates)),rownames(cosmic$Serum_M_Sphingolipid))
  cosmic$Serum_M_Targeted <- cbind(cosmic$Serum_M_Biocrates[sample,],cosmic$Serum_M_Sphingolipid[sample,],cosmic$Serum_M_Sphingolipid[sample,])
  
  
  return(cosmic)
}



datasummary_format <- function (datasummary_raw,n) {
        # format datasummary_raw to datasummary
        # n is the start column of data
        t1<-1-is.na(datasummary_raw[,n:ncol(datasummary_raw)])
        t2 <- data.frame()
        t3 <- data.frame()
        for (i in 1:nrow(datasummary_raw)) {
                if (datasummary_raw$gender[[i]] == "male") {
                        g <- "M" 
                } else {
                        g <- "F"
                }
                if (datasummary_raw$smoking[[i]] == "non-smoker") s <- "N"  
                if (datasummary_raw$smoking[[i]] == "smoker") s <- "S"
                if (datasummary_raw$smoking[[i]] == "Ex-smoker") s <- "E"
                if (datasummary_raw$diagnosis[[i]] == "healthy") {
                        d <- "H" 
                } else {
                        d <- "C"
                }
                t3[i,1] <- paste (s,d,sep="")
                t2[i,1] <- paste ("COSMIC",datasummary_raw$COSMIC_ID[[i]], g, s, d, sep="_")
        }
        colnames(t2) <- "barcode"
        cgroup <- factor(t3[,1],levels=c("NH","SH","SC","EC"))
        datasummary <- data.frame()
        datasummary <- cbind(t2,cgroup,datasummary_raw[,2:n-1], t1)
        rownames(datasummary) <- t2[,1];
        rm(t1,t2, t3, d,g,i,s)
        #levels(datasummary$cgroup) <- c("NH","SH","SC","EC")
        datasummary
}

iTRAQ_impute <- function(data_pre,fk,fpath) {
        
        # data_pre <- cosmic$BAL_P_iTRAQ_ONE
        # setwd("./BAL_P_iTRAQ_ONE")
        # data_pre <- cosmic$BAL_P_iTRAQ_ALL
        # setwd("./BAL_P_iTRAQ_ALL")
        setwd(fpath)
        # histgram of percent of unavailable values
        rmis <- rowSums(is.na(data_pre))/ncol(data_pre) # percent of unavaliable values for each feature
        cmis <- colSums(is.na(data_pre))/nrow(data_pre) # percent of unavaliable values for each sample
        pdf("HistUnavailable.pdf",width=10,height=5)
        par(mfcol=c(1,2))
        hist(cmis,xlab="Unavailable values %",ylab="# Samples",main="Summary by samples")
        hist(rmis, xlab="Unavailable values %",ylab="# Features",main="Summary by features")
        dev.off()
        
        # subdata with no missing values
        subdata <- data_pre[rmis==0,] 
        cmis2 <- colSums(is.na(subdata))/nrow(subdata)
        range(cmis2) # zero, so no missing values in this subdata
        data_pre <- data_pre[rmis<=0.25,] # delete rows with missing value rate >0.25
        rmis <- rowSums(is.na(data_pre))/ncol(data_pre) # percent of unavaliable values for each feature
        
        # estimation of parameter K
        library(hydroGOF)
        library(VIM)
        subdata<-as.matrix(subdata) 
        N <- nrow(subdata)*ncol(subdata) # number of values in the data
        nr <- nrow(subdata) # number of features
        nc <- ncol(subdata) # number of samples
        m <- N*(sum(is.na(data_pre))/(nrow(data_pre)*ncol(data_pre))) # m is the total number of unavailable values
        rn = 10 # number of permutation
        rk <- seq(5,25,by=5) # sequence for K
        rc <- length(rk) # length of sequence of K
        result <- array(0,dim=c(20,rn,rc))
        subv <- as.vector(subdata)
        for (i in 1:rn) {
                print(i)
                t1 <- sample(N) # sampling
                t1 <- t1[1:m] # m values selected to be as NA
                s1 <- subv
                s1[t1] = NA
                sm <- matrix(s1,nr,nc)
                sm <- data.frame(sm)
                sim <- subv[t1]
                for (j in 1:rc) {
                        t2 <-kNN(sm, variable = colnames(sm),numFun=mean, metric = NULL,k = rk[j],trace=FALSE) # KNN impute
                        t22<-data.matrix(t2[,1:nc])
                        t3<-as.vector(t22[,1:nc])
                        obs <- t3[t1]
                        rr <- gof(sim=sim,obs=obs) # gof functions
                        for (x in 1:20) {
                                result[x,i,j] <- rr[x]
                        }
                }
        }
        
        # result is a 20-rn-rc array, result[x,i,j] is the gof value for the ith K value, jth permutation and x in gof.
        # boxplot for gof result for different K
        fyname<-rownames(rr)
        fyname2<-fyname
        fyname2[5]<-"NRMSE"
        fyname2[6]<-"PBIAS"
        filename <- "raw_data_gof.pdf"
        pdf(filename)
        par(mfrow=c(3,2),mar=c(5,5,3,1))
        for (i in c(12,14,16,17,18)){ 
                boxplot(result[i,1:rn,],boxwex = 0.5,ylim=c(0,1),xlab="K",ylab=fyname[i],xaxt="n",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
                axis(1,at=c(1:5),labels=as.character(rk),cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
        }
        dev.off()
        
        # estimation of parameter for filter features - percent of unavailable values
        p<-c(0.01,0.05,seq(0.1,0.25,by=0.05))
        rn=10
        k=fk
        result <- array(0,dim=c(20,rn,length(p)))
        
        for (j in 1:length(p)) {
                print(j)
                m <- floor(p[j]*N)
                for (i in 1:rn) {
                        t1 <- sample(N)
                        t1 <- t1[1:m]
                        s1 <- subv
                        s1[t1] = NA
                        sm <- matrix(s1,nr,nc)
                        sm <- data.frame(sm)
                        sim <- subv[t1]
                        t2 <-kNN(sm, variable = colnames(sm),numFun=mean, metric = NULL,k = k,trace=FALSE)
                        t22<-data.matrix(t2[,1:nc])
                        t3<-as.vector(t22[,1:nc])
                        obs <- t3[t1]
                        rr <- gof(sim=sim,obs=obs)
                        for (x in 1:20) {
                                result[x,i,j] <- rr[x]
                        }
                }
        }
        
        
        # boxplot for gof result for different percent of unavailable values
        filename <- "raw_data_gof_percent.pdf"
        pdf(filename)
        par(mfrow=c(3,2),mar=c(5,5,3,1))
        for (i in c(12,14,16,17,18)){ 
                boxplot(result[i,1:rn,1:5],boxwex = 0.5,ylim=c(0,1),xlab="% unavailable values",ylab=fyname[i],xaxt="n",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
                axis(1,at=c(1:5),labels=as.character(p[1:5]),cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
        }
        dev.off()
        
        # filter and imputation
        fp = 0.25
        row_to_keep <- which(rmis <= fp)
        row_to_delete<-which(rmis>fp)
        data_pre2<-data_pre[row_to_delete,]
        data_pre <- data_pre[row_to_keep,]
        #fk = 10 
        data_pre<-data.frame(data_pre)
        data <- kNN(data_pre, variable = colnames(data_pre),numFun=mean, metric = NULL,k = fk,trace=FALSE)
        data <- data[,1:ncol(data_pre)]
        
        # distribution comparison before and after imputation
        filename <- "boxplot.pdf"
        pdf(filename)
        par(mfcol=c(2,1),mar=c(5,5,3,1))
        boxplot(data_pre,xaxt="n",xlab="Samples",ylab="log2(values)",main='Before imputation')
        boxplot(data,xaxt="n",xlab="Samples",ylab="log2(values)",main="After imputation")
        dev.off()
        return(as.matrix(data))
}
# BAL_P_iTRAQ_ONE_impute <- data
# BAL_P_iTRAQ_ALL_impute <- data

# clinic bioconductor, clinical data matrix using sample for rows prepared by Vincenzo
# subject <- "clinic_bioconductor"
# data <- read.csv(par[subject,"fname"])
# c <- as.integer(as.matrix(data[,1]))
# c2 <- match(c,cosmic$datasummary$COSMIC_ID)
# data$barcode[is.na(c2)==0] <- cosmic$datasummary$barcode[c2][is.na(c2)==0]
# eval(parse(text = paste("cosmic$",subject,"<-data",sep="")))
