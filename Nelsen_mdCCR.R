## Script for mdCCR - Nelsen 4.2.12 Copula ##

rm(list=ls())

library(Matrix)
library(numDeriv)
library(doParallel)
library(foreach)

nprocs=20
print( paste( as.character(detectCores()), "cores detected" ) );
cl <- makePSOCKcluster(nprocs)
doParallel::registerDoParallel(cl)
print( paste( as.character(getDoParWorkers() ), "cores registered" ) )


## Functions for Nelson 4.2.12
gfun <- function(x){log(x/(1-x))}
gfun.inv <- function(x){exp(x)/(1+exp(x))}
hfun <- function(x){log(x/(1-x))}

## Functions for likelihood
prepare.Uim <- function(a, tm, c0){ # Computing Uim without adjusting the cdf
  
  n=nrow(a)
  a.cdf1 <- ecdf(a[a$y1<1,1])
  a.cdf2 <- ecdf(a[a$y2<1,2])
  cdf1.value <- a.cdf1(a[,1])*c0[1]#*n/(n+1)
  cdf2.value <- a.cdf2(a[,2])*c0[2]#*n/(n+1)
  
  w <- 1*(1*outer(cdf1.value, c0[1], "<") +1*outer(cdf2.value, c0[1], "<") ==2) 
  cdf1.value <- cdf1.value[which(w==1)]
  cdf2.value <- cdf2.value[which(w==1)]
  
  #counts in each category
  wt.temp <- 1*(1*outer(cdf1.value, tm, "<=") + 1*outer(cdf2.value, tm, "<=")==2) 
  wt <- cbind(wt.temp[, 1], wt.temp[, -1] - wt.temp[, -m])     
  
  invisible(wt)
}


L11 <- function(thets,a,x,tm,c0){
  
  alp1 <- thets[1]
  bet1 <- thets[2]
  wt <- prepare.Uim(a, tm, c0)
  ht <- hfun(tm)
  coef <- alp1 + bet1 * x
  
  ght <- gfun.inv(coef + ht) # Psi(tm)
  dght <- c(ght[1], diff(ght))
  
  lpr <- ifelse(dght > 0, log(dght), Inf) 
  lik <- sum(wt %*% lpr)
  
  return(lik)
  
}

L12<-function(x,c0,thets){
  coef<-thets[1] + (thets[2]*x)
  fn<- c0 - gfun.inv(coef + hfun(c0))
  lpr <- ifelse(fn > 0, log(fn), Inf)
  return(lpr)
}

L22<-function(x,c0,thets){
  coef<-thets[1] + (thets[2]*x)
  fn<- 1 - (2*c0) + gfun.inv(coef + hfun(c0))
  lpr <- ifelse(fn > 0, log(fn), Inf)
  return(lpr)
}

likfun<-function(thets){ #thets= (alpha, beta, c1^1,c1^2)
  
  if(thets[1]< -log(2) | thets[1]> 0 | thets[2]< -log(2) | thets[2]> log(2)| any(thets[-c(1:2)]<0) | any(thets[-c(1:2)]>1)){
    return(Inf)
  }
  
  r <- levels(factor(C11$Rep)) # no. of replicates per wf
  x.unique <- levels(factor(C11$x)) # no. of workflows
  
  cc<-thets[-c(1:2)]
  cc2<-thets[-c(1:2)]
  
  l11<-0
  l12<-0
  l21<-0
  l22<-0
  
  for(k in 1:length(x.unique)){
    
    if(cc[k]>0.9999){
      tm <-seq(0.01, .9999, length.out=m)
    }else{
      tm <-seq(0.01, cc[k], length.out=m)
    }
    
    for(l in 1:length(r)){
      
      x1<- as.numeric(x.unique[k])
      r1<-as.numeric(r[l])
      
      lik0<-L11(thets, C11[C11$x==x1 & C11$Rep==r1,1:2], x1, tm, c(cc[k],cc2[k]))
      l11<-l11+lik0
      
      lik1<-ifelse(n12[l,k]>0,(L12(x1,cc[k],thets))*n12[l,k],0)
      l12= l12+lik1
      
      lik2<-ifelse(n21[l,k]>0,(L12(x1,cc[k],thets))*n21[l,k],0)
      l21= l21+lik2
      
      lik3<-ifelse(e_n22[k]>0,(L22(x1,cc[k],thets))*(e_n22[k]),0)
      l22= l22+lik3
    }
  }
  return(l11+l12+l21+l22) #log (L11 L12 L21 L22)
}

## generate data from Archimedean copula from (4.2.12) in Nelsen (2006)
generatefun <- function(t, thet){(1/t-1)^thet}
generatefun.inv <- function(t, thet){(1+t^(1/thet))^(-1)}
kcfun <- function(t, thet){t*(1+1/thet)- t^2*1/thet}
kcfun.inv <- function(t, thet){1/2*(1+thet-sqrt((1+thet)^2-4*thet*t))}

datSim_Nelson4212 <- function(n, thet){
  
  v1 <- runif(n, 0, 1)
  v2 <- runif(n, 0, 1)
  t <- kcfun.inv(v2, thet)
  y1 <- generatefun.inv(v1*generatefun(t, thet), thet)
  y2 <- generatefun.inv((1-v1)*generatefun(t, thet), thet)
  
  dat <- cbind(y1, y2)
  colnames(dat) <- c("y1", "y2")
  
  return(dat)
}


## Parameters for Nelson
theta1=2
theta2=1.7
num_sim=10000
num_samples=100
boot_samples=400
m=30
num_rep=3
num_wf=2

## Combinations of c1 and c2
r<-c(.9,.8,.7,.6)
grid<-expand.grid(r,r)
colnames(grid)<-c("c1","c2")
grid=grid[grid$c1==grid$c2,]
rownames(grid)=1:nrow(grid)
grid=grid[3,]

par.est<-matrix(NA, ncol=16, nrow=nrow(grid))
par.est[,1]<- (-log (2))/theta1 # alpha
par.est[,2]<- (-log (2)/theta2)+(log (2)/theta1) #beta
colnames(par.est)<-c("alpha","beta","c1","c1.est1","c1.est2",
                     "alpha.est","alpha.se","alpha.sd",
                     "beta.est","beta.se","beta.sd",
                     "nm1","nm1.est","nm2","nm2.est","beta.power")
ptm=proc.time()
for(i in 1:nrow(grid)){
  c1=grid[i,1]
  c2=c1
  par.est[i,3]=c1
  
  param_est <- foreach( j = 1:num_samples,
                        .combine="cbind") %dopar% {
                          
                          n11<-matrix(NA,nrow=num_rep, ncol= num_wf)
                          n12p<-matrix(NA,nrow=num_rep, ncol= num_wf)
                          n12c<-matrix(NA,nrow=num_rep, ncol= num_wf)
                          n21<-matrix(NA,nrow=num_rep, ncol= num_wf)
                          n22<-matrix(NA,nrow=num_rep, ncol= num_wf)
                          nm<-matrix(NA,nrow=num_rep, ncol= num_wf)
                          
                          C11=NULL
                          
                          for(nr in 1: num_rep){
                            
                            wf1<- datSim_Nelson4212(num_sim,theta1)
                            wf2<- datSim_Nelson4212(num_sim,theta2)
                            
                            y1<-wf1[,1]
                            y2<-wf1[,2]
                            num<- nrow(wf1)
                            
                            cdf1 <- ecdf(y1)
                            cdf2 <- ecdf(y2)
                            cdf1.val <- cdf1(y1)*num/(num+1)
                            cdf2.val <- cdf2(y2)*num/(num+1)
                            
                            ## Partitioning data
                            wt1 <- 1*(1*outer(cdf1.val, c1, "<=") +1*outer(cdf2.val, c1, "<=") ==2) 
                            n11[nr,1] <- length(y1[which(wt1==1)]) #Category 11 (all complete)
                            
                            wt2 <- 1*(1*outer(cdf1.val, c1, "<=") +1*outer(cdf2.val, c1, ">")
                                      +1*outer(cdf2.val, c2, "<=") ==3) 
                            n12c[nr,1] <- length(y1[which(wt2==1)]) #Category 12 (all complete)
                            
                            wt3 <- 1*(1*outer(cdf1.val, c1, "<=") +1*outer(cdf2.val, c2, ">") ==2) 
                            n12p[nr,1] <- length(y1[which(wt3==1)]) #Category 12 (Y2 missing)
                            
                            wt4 <- 1*(1*outer(cdf1.val, c1, ">") +1*outer(cdf2.val, c1, "<=") ==2) 
                            n21[nr,1] <- length(y1[which(wt4==1)]) #Category 21 (Y1 missing)
                            
                            wt5 <- 1*(1*outer(cdf1.val, c1, ">") +1*outer(cdf2.val, c1, ">") 
                                      +1*outer(cdf2.val, c2, "<=") ==3) 
                            n22[nr,1] <- length(y1[which(wt5==1)])#Category 22 (Y1 missing)
                            
                            wt6 <- 1*(1*outer(cdf1.val, c1, ">") +1*outer(cdf2.val, c2, ">") ==2)
                            nm[nr,1] <- length(y1[which(wt6==1)]) #Category 22 (Completely missing)
                            
                            w1 <- 1*(1*outer(cdf1.val, c1, "<=") ==1) 
                            y1[which(w1==0)]=1
                            w2 <- 1*(1*outer(cdf2.val, c2, "<=") ==1) 
                            y2[which(w2==0)]=1
                            C11a=as.data.frame(cbind(y1,y2,0,nr))
                            colnames(C11a)=c("y1","y2","x","Rep")
                            
                            y1<-wf2[,1]
                            y2<-wf2[,2]
                            num<- nrow(wf2)
                            
                            cdf1 <- ecdf(y1)
                            cdf2 <- ecdf(y2)
                            cdf1.val <- cdf1(y1)*num/(num+1)
                            cdf2.val <- cdf2(y2)*num/(num+1)
                            
                            ## Partitioning data
                            wt1 <- 1*(1*outer(cdf1.val, c1, "<=") +1*outer(cdf2.val, c1, "<=") ==2) 
                            n11[nr,2] <- length(y1[which(wt1==1)]) #Category 11 (all complete)
                            
                            wt2 <- 1*(1*outer(cdf1.val, c1, "<=") +1*outer(cdf2.val, c1, ">")
                                      +1*outer(cdf2.val, c2, "<=") ==3) 
                            n12c[nr,2] <- length(y1[which(wt2==1)]) #Category 12 (all complete)
                            
                            wt3 <- 1*(1*outer(cdf1.val, c1, "<=") +1*outer(cdf2.val, c2, ">") ==2) 
                            n12p[nr,2] <- length(y1[which(wt3==1)]) #Category 12 (Y2 missing)
                            n12=n12c+n12p
                            
                            wt4 <- 1*(1*outer(cdf1.val, c1, ">") +1*outer(cdf2.val, c1, "<=") ==2) 
                            n21[nr,2] <- length(y1[which(wt4==1)]) #Category 21 (Y1 missing)
                            
                            wt5 <- 1*(1*outer(cdf1.val, c1, ">") +1*outer(cdf2.val, c1, ">") 
                                      +1*outer(cdf2.val, c2, "<=") ==3) 
                            n22[nr,2] <- length(y1[which(wt5==1)])#Category 22 (Y1 missing)
                            
                            wt6 <- 1*(1*outer(cdf1.val, c1, ">") +1*outer(cdf2.val, c2, ">") ==2)
                            nm[nr,2] <- length(y1[which(wt6==1)]) #Category 22 (Completely missing)
                            
                            w1 <- 1*(1*outer(cdf1.val, c1, "<=") ==1) 
                            y1[which(w1==0)]=1
                            w2 <- 1*(1*outer(cdf2.val, c2, "<=") ==1) 
                            y2[which(w2==0)]=1
                            C11b=as.data.frame(cbind(y1,y2,1,nr))
                            colnames(C11b)=c("y1","y2","x","Rep")
                            
                            C11=as.data.frame(rbind(C11,C11a,C11b))
                            colnames(C11)=c("y1","y2","x","Rep")
                          }
                          
                          n_obs=round(colMeans(n11+n12+n21))
                          c1_est=c(c1,c1)
                          alp_est=-0.34
                          bet_est=-0.06
                          thet=c(alp_est,bet_est,c1_est)
                          xx=unique(C11$x)
                          e_n22=rep(NA,2)
                          e_n22[1]=n_obs[1]*(1-2*c1_est[1]+ gfun.inv(alp_est + (bet_est*xx[1])+ hfun(c1_est[1])))/(2*c1_est[1]- gfun.inv(alp_est + (bet_est*xx[1]) + hfun(c1_est[1])))
                          e_n22[2]=n_obs[2]*(1-2*c1_est[2]+ gfun.inv(alp_est + (bet_est*xx[2])+ hfun(c1_est[2])))/(2*c1_est[2]- gfun.inv(alp_est + (bet_est*xx[2]) + hfun(c1_est[2])))
                          fit<- optim(thet,fn = function(thets){-likfun(thets)},hessian=TRUE)
                          alp_est=fit$par[1]
                          bet_est=fit$par[2]
                          c1_est=fit$par[3:4]
                          while(likfun(c(alp_est,bet_est,c1_est))-likfun(thet)>0.1){
                            thet=c(alp_est,bet_est,c1_est)
                            e_n22[1]=n_obs[1]*(1-2*c1_est[1]+ gfun.inv(alp_est + (bet_est*xx[1])+ hfun(c1_est[1])))/(2*c1_est[1]- gfun.inv(alp_est + (bet_est*xx[1]) + hfun(c1_est[1])))
                            e_n22[2]=n_obs[2]*(1-2*c1_est[2]+ gfun.inv(alp_est + (bet_est*xx[2])+ hfun(c1_est[2])))/(2*c1_est[2]- gfun.inv(alp_est + (bet_est*xx[2]) + hfun(c1_est[2])))
                            fit<- optim(thet,fn = function(thets){-likfun(thets)},hessian=TRUE)
                            alp_est=fit$par[1]
                            bet_est=fit$par[2]
                            c1_est=fit$par[3:4]
                          }
                          
                          a.est=alp_est
                          b.est=bet_est
                          c.est1=c1_est[1]
                          c.est2=c1_est[2]

                          nm1.est=e_n22[1]
                          nm2.est=e_n22[2]
                          
                          dat=C11
                          
                          alp.est=NULL
                          bet.est=NULL
                          c1.est=NULL
                          c2.est=NULL

                          for( k in 1:boot_samples){
                            
                            C11=NULL
                            for(nr in 1: num_rep){
                              temp=dat[dat$x==0 & dat$Rep==nr,]
                              C11a=temp[sample(1:num_sim,size=num_sim, replace=T),]
                              
                              temp=dat[dat$x==1 & dat$Rep==nr,]
                              C11b=temp[sample(1:num_sim,size=num_sim, replace=T),]
                              
                              C11=rbind(C11,C11a,C11b)
                              
                              n11[nr,]=c(nrow(C11[C11$y1<1 & C11$y2<1 & C11$x==0,]),
                                         nrow(C11[C11$y1<1 & C11$y2<1 & C11$x==1,]))
                              n12[nr,]=c(nrow(C11[C11$y1<1 & C11$y2==1 & C11$x==0,]),
                                         nrow(C11[C11$y1<1 & C11$y2==1 & C11$x==1,]))
                              n21[nr,]=c(nrow(C11[C11$y1==1 & C11$y2<1 & C11$x==0,]),
                                         nrow(C11[C11$y1==1 & C11$y2<1 & C11$x==1,]))
                            }
                            
                            C11=as.data.frame(C11)
                            n_obs=round(colMeans(n11+n12+n21))
                            
                            thet=c(alp_est,bet_est,c1_est)
                            
                            fit<- optim(thet,fn = function(thets){-likfun(thets)})
                            alp_est=fit$par[1]
                            bet_est=fit$par[2]
                            c1_est=fit$par[3:4]
                            while(likfun(c(alp_est,bet_est,c1_est))-likfun(thet)>0.01){
                              thet=c(alp_est,bet_est,c1_est)
                              e_n22[1]=n_obs[1]*(1-2*c1_est[1]+ gfun.inv(alp_est + (bet_est*xx[1])+ hfun(c1_est[1])))/(2*c1_est[1]- gfun.inv(alp_est + (bet_est*xx[1]) + hfun(c1_est[1])))
                              e_n22[2]=n_obs[2]*(1-2*c1_est[2]+ gfun.inv(alp_est + (bet_est*xx[2])+ hfun(c1_est[2])))/(2*c1_est[2]- gfun.inv(alp_est + (bet_est*xx[2]) + hfun(c1_est[2])))
                              fit<- optim(thet,fn = function(thets){-likfun(thets)})
                              alp_est=fit$par[1]
                              bet_est=fit$par[2]
                              c1_est=fit$par[3:4]
                            }
                            alp.est=c(alp.est,alp_est)
                            bet.est=c(bet.est,bet_est)
                            c1.est=c(c1.est,c1_est[1])
                            c2.est=c(c1.est,c1_est[2])
                          }
                          
                          a.ese <- sd(alp.est)
                          b.ese <- sd(bet.est)
                          
                          CI.low <- b.est - qnorm(.975) * b.ese
                          CI.up <- b.est + qnorm(.975) * b.ese
                          if(0 > CI.up | 0 < CI.low){
                            rejects=1
                          } else{
                            rejects=0
                          }

                          res1=rbind(c.est1, c.est2,
                                     a.est, a.ese,
                                     b.est, b.ese,
                                     nm[1],nm1.est,nm[2],nm2.est,
                                     rejects)
                          return(res1)
                        }
  par.est[i,4]=mean(param_est[1,])
  par.est[i,5]=mean(param_est[2,])
  par.est[i,6]=mean(param_est[3,])
  par.est[i,7]=mean(param_est[4,])
  par.est[i,8]=sd(param_est[3,])
  par.est[i,9]=mean(param_est[5,])
  par.est[i,10]=mean(param_est[6,])
  par.est[i,11]=sd(param_est[5,])
  par.est[i,12]=mean(param_est[7,])
  par.est[i,13]=mean(param_est[8,])
  par.est[i,14]=mean(param_est[9,])
  par.est[i,15]=mean(param_est[10,])
  par.est[i,16]=sum(param_est[11,])/num_samples
}
proc.time()-ptm
out=round(par.est, digits=6)
out
write.csv(out, file="Nelson_theta=2,1.7_Rep3_grid3.csv")


