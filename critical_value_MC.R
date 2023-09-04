# to get where this script is saved in the local disk

list.of.packages <- c("rstudioapi")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(rstudioapi)

working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

setwd(working_dir)

source("packages and dists.R") # calling the code "packages and dists.R" which contains the distribution functions and also the packages to install

# 1. PE3

sh <- c(-0.4, -0.3, -0.2, -0.1,-0.05,0.4, 0.3,0.2) #shape parameters
ns <- c(10,15,20,25,32) #sample size 
criticals.pe3.ad <- c()
criticals.pe3.ks <- c()

for (uu in ns){
  
  n.mc<-30000        ### Number of Simulations
  AD.pe<-matrix(NA,nrow = n.mc,ncol = length(sh))   ### Initialisation of array with AD-statistic
  KS.pe<-matrix(NA,nrow = n.mc,ncol = length(sh))   ### Initialisation of array with KS-statistic
  n<-uu            ### Length of the series
  J<-1:n
  
  for(k in 1:n.mc){
    for (p in seq(1,length(sh),1)){
      X <- sort(runif(n))
      X <- qpe3(X,shape=sh[p],scale=1,location = 0) #for consistency and generalization we set the scale and location paramters to be 1 and 0, respectively.
      
      X<-sort(X)
      X<- X/sd(X,na.rm = TRUE)
      
      tryCatch({
        
        lmr2 <- lmoms(X)
        para2 <- parpe3(lmr2,checklmom = TRUE)
        
        a2=para2[[2]][['mu']]
        b2=para2[[2]][['sigma']]
        c2=para2[[2]][['gamma']]
        
        fitpe3 <- fitdist(X, distr="pe3", method="mle", start = list(shape=c2,scale=b2,location=a2))
        
        
        jpe3<-gofstat(fitpe3)
        
        AD.pe[k,p] <- unname(jpe3$ad)
        
        KS.pe[k,p] <- unname(jpe3$ks)
        
      },error=function(e){
        
        #cat(k, " ", p)
        
      })
      
      
    }
  }
  
  ll<- c()
  ll<-append(ll,AD.pe[,1])
  ll<-append(ll,AD.pe[,2])
  ll<-append(ll,AD.pe[,3])
  ll<-append(ll,AD.pe[,4])
  ll<-append(ll,AD.pe[,5])
  
  ll<-append(ll,AD.pe[,6])
  ll<-append(ll,AD.pe[,7])
  ll<-append(ll,AD.pe[,8])
  
  
  ll <- na.omit(ll)
  crit.val1<-quantile(ll,probs=c(0.95))
  
  criticals.pe3.ad <- append(criticals.pe3.ad,unname(crit.val1))
  
  ##################### KS ################
  ls<- c()
  ls<-append(ls,KS.pe[,1])
  ls<-append(ls,KS.pe[,2])
  ls<-append(ls,KS.pe[,3])
  ls<-append(ls,KS.pe[,4])
  ls<-append(ls,KS.pe[,5])
  
  ls<-append(ls,KS.pe[,6])
  ls<-append(ls,KS.pe[,7])
  ls<-append(ls,KS.pe[,8])
  
  ls <- na.omit(ls)
  
  crit.val2<-quantile(ls,probs=c(0.95))
  
  criticals.pe3.ks <- append(criticals.pe3.ks,unname(crit.val2))
  
  
}

#2. GEV 

sh <- c(-0.4, -0.3, -0.2, -0.1,-0.05,0.4, 0.3,0.2) # shape parameters

ns <- c(10,15,20,25,32) # sample size
criticals.ge.ad <- c()
criticals.ge.ks <- c()

for (uu in ns){
  
  n.mc<-30000        ### Number of Simulations
  AD.ge<-matrix(NA,nrow = n.mc,ncol = length(sh))   ### Initialisation of array with AD-statistic
  KS.ge<-matrix(NA,nrow = n.mc,ncol = length(sh))   ### Initialisation of array with AD-statistic
  n<-uu            ### Length of the series
  J<-1:n
  
  for(k in 1:n.mc){
    for (p in seq(1,length(sh),1)){
      
      X <- sort(runif(n))
      
      X<- qgext(X,shape = sh[p],scale=1,location = 0)
      
      X<-sort(X)   ### Sort the data
      
      X<- X/sd(X,na.rm = TRUE)
      
      tryCatch({
        
        lmr3 <- lmoms(X,nmom = 5)
        para3 <- pargev(lmr3,checklmom = TRUE)
        
        a3=para3[[2]][['xi']]
        b3=para3[[2]][['alpha']]
        c3=para3[[2]][['kappa']]
        
        #fitting
        fitgev <- fitdist(X, "gext",method='mle',start = list(shape=c3,scale=b3,location=a3))
        
        jext<-gofstat(fitgev)
        
        AD.ge[k,p] <- unname(jext$ad)
        
        KS.ge[k,p] <- unname(jext$ks)
        
      },error=function(e){
        
        #cat(k, " ", p)
        
      })
      
      
    }
    
  }
  
  ll<- c()
  ll<-append(ll,AD.ge[,1])
  ll<-append(ll,AD.ge[,2])
  ll<-append(ll,AD.ge[,3])
  ll<-append(ll,AD.ge[,4])
  ll<-append(ll,AD.ge[,5])
  ll<-append(ll,AD.ge[,6])
  ll<-append(ll,AD.ge[,7])
  ll<-append(ll,AD.ge[,8])
  
  
  ll <- na.omit(ll)
  
  crit.val1<-quantile(ll,probs=c(0.95))
  
  criticals.ge.ad <- append(criticals.ge.ad,unname(crit.val1))
  ##################### KS ################
  ls<- c()
  ls<-append(ls,KS.ge[,1])
  ls<-append(ls,KS.ge[,2])
  ls<-append(ls,KS.ge[,3])
  ls<-append(ls,KS.ge[,4])
  ls<-append(ls,KS.ge[,5])
  ls<-append(ls,KS.ge[,6])
  ls<-append(ls,KS.ge[,7])
  ls<-append(ls,KS.ge[,8])
  
  ls<- na.omit(ls)
  crit.val2<-quantile(ls,probs=c(0.95))
  
  criticals.ge.ks <- append(criticals.ge.ks,unname(crit.val2))
  
  
}

#3. Genlog 

sh <- c(-0.4, -0.3, -0.2, -0.1,-0.05,0.4, 0.3,0.2) #shape parameters

ns <- c(10,15,20,25,32) # sample size
criticals.gl.ad <- c()
criticals.gl.ks <- c()

for (uu in ns){
  
  n.mc <-30000        ### Number of Simulations
  AD.gl <-matrix(NA,nrow = n.mc,ncol = length(sh))   ### Initialisation of array with AD-statistic
  KS.gl <-matrix(NA,nrow = n.mc,ncol = length(sh))   ### Initialisation of array with AD-statistic
  n<-uu            ### Length of the series
  J<-1:n
  
  for(k in 1:n.mc){
    for (p in seq(1,length(sh),1)){
      
      X <- sort(runif(n))
      
      X <- qgenlog(X,shape = sh[p],scale=1,location = 0)
      
      
      X<-sort(X)   ### Sort the data
      
      tryCatch({
        
        lmr <- lmoms(X,nmom = 5)
        para <- parglo(lmr,checklmom = TRUE)
        
        a=para[[2]][['xi']]
        b=para[[2]][['alpha']]
        c=para[[2]][['kappa']]
        
        
        fitglo <- fitdist(X, "genlog",method='mle',start = list(shape=c,scale=b,location=a))
        
        jglo<-gofstat(fitglo)
        
        AD.gl[k,p] <- unname(jglo$ad)
        
        KS.gl[k,p] <- unname(jglo$ks)
        
      },error=function(e){
        
        #cat(k, " ", p)
        
      })
      
    }
    
  }
  
  ############### AD ####################
  ll<- c()
  ll<-append(ll,AD.gl[,1])
  ll<-append(ll,AD.gl[,2])
  ll<-append(ll,AD.gl[,3])
  ll<-append(ll,AD.gl[,4])
  ll<-append(ll,AD.gl[,5])
  
  ll<-append(ll,AD.gl[,6])
  ll<-append(ll,AD.gl[,7])
  ll<-append(ll,AD.gl[,8])
  
  ll <- na.omit(ll)
  
  crit.val1<-quantile(ll,probs=c(0.95))
  
  
  criticals.gl.ad <- append(criticals.gl.ad,unname(crit.val1))
  ##################### KS ################
  ##################### KS ################
  ls<- c()
  ls<-append(ls,KS.gl[,1])
  ls<-append(ls,KS.gl[,2])
  ls<-append(ls,KS.gl[,3])
  ls<-append(ls,KS.gl[,4])
  ls<-append(ls,KS.gl[,5])
  
  ls<-append(ls,KS.gl[,6])
  ls<-append(ls,KS.gl[,7])
  ls<-append(ls,KS.gl[,8])
  
  ls <- na.omit(ls)
  
  crit.val2<-quantile(ls,probs=c(0.95))
  
  criticals.gl.ks <- append(criticals.gl.ks,unname(crit.val2))
  
  
}

#printing the critical values

#1. for PE3

cat(criticals.pe3.ad)
cat(criticals.pe3.ks)

#1. for Genlog

cat(criticals.gl.ad)
cat(criticals.gl.ks)

#1. for GEV

cat(criticals.ge.ad)
cat(criticals.ge.ks)
