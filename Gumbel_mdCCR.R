## Script for Inverted mdCCR - Gumbel Hougaard Copula ##

rm(list=ls())

library(gumbel)
library(Matrix)
library(doParallel)
library(foreach)

nprocs=20
print( paste( as.character(detectCores()), "cores detected" ) );
cl <- makePSOCKcluster(nprocs)
doParallel::registerDoParallel(cl)
print( paste( as.character(getDoParWorkers() ), "cores registered" ) )

## Functions for Gumbel-Hougaard
gfun <- function(x){log(x)}
gfun.inv <- function(x){exp(x)}
hfun <- function(x){log(x)}

## Functions for log likelihood
prepare.Uim <- function(a, tm, c0){
  
  a.cdf1 <- ecdf(a[a$y1>0,1])
  a.cdf2 <- ecdf(a[a$y2>0,2])
  cdf1.value <- (a.cdf1(a[,1])*(1-c0[1])) + c0[1]
  cdf2.value <- (a.cdf2(a[,2])*(1-c0[2])) + c0[2]
  
  w <- 1*(1*outer(cdf1.value, c0[1], ">") +1*outer(cdf2.value, c0[1], ">") ==2) 
  cdf1.value <- cdf1.value[which(w==1)]
  cdf2.value <- cdf2.value[which(w==1)]
  
  # get counts in each category
  wt.temp <- 1*(1*outer(cdf1.value, tm, ">") + 1*outer(cdf2.value, tm, ">")==2) 
  wt <- cbind(wt.temp[, -m]- wt.temp[, -1] , wt.temp[, m])     # n*m
  
  invisible(wt)
}

L11 <- function(thets,a,x,tm, c0){
  
  alp <- thets[1]
  bet1 <- thets[2]
  
  m <- length(tm)
  
  wt <- prepare.Uim(a, tm, c0) # getting Uim
  
  coef <- alp + bet1 * x 
  ht <- hfun(tm)
  
  ght <- 1 - (2*tm) + gfun.inv(coef * ht) #psi*(tm)
  
  if(any(ght<0)){
    return(Inf)
  }
  
  dght <- c(ght[-m]-ght[-1], ght[m]) # Psi*(tm-1) - Psi*(tm)
  
  lpr <- ifelse(dght > 0, log(dght), Inf) 
  lik <- sum(wt %*% lpr) #log L11
  
  return(lik)
}

L12<-function(x,c0,thets){
  coef<-thets[1] + (thets[2]*x)
  fn<- c0 - gfun.inv(coef * hfun(c0))
  lpr <- ifelse(fn > 0, log(fn), Inf)
  return(lpr) #log L12 or L21
}

L22<-function(x,c0,thets){
  coef<-thets[1] + (thets[2]*x)
  fn<- gfun.inv(coef * hfun(c0))
  lpr <- ifelse(fn > 0, log(fn), Inf)
  return(lpr) #log L22
}

likfun<-function(thets){ #thets= (alpha, beta, nm0,nm1)
  
  if(thets[1]<1 | thets[1]>2 | thets[2]< -1| thets[2]>1 | any(thets[-c(1:2)]<0) | any(thets[-c(1:2)]>1)){ # avoiding negative values for nm
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
  const<-0
  
  for(k in 1:length(x.unique)){
    
    if(cc[k]>0){
      tm <-seq(cc[k], .99, length.out=m)
    }else{
      tm <-seq(0.0001, .99, length.out=m)
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

## Parameters from Gumbel-Hougaard. Here alpha represents theta.
theta1=2
theta2=1.7
num_sim=10000
num_samples=100
boot_samples=400
m=30
num_rep=3
num_wf=2

## Combinations of c1 and c2
r<-c(.1,.2,.3,.4)
grid<-expand.grid(r,r)
colnames(grid)<-c("c1","c2")
grid=grid[grid$c1==grid$c2,]
rownames(grid)=1:nrow(grid)
grid=grid[1,]

par.est<-matrix(NA, ncol=16, nrow=nrow(grid))
par.est[,1]<-2^(1/theta1) #true alpha
par.est[,2]<-(2^(1/theta2))-(2^(1/theta1)) #true beta
colnames(par.est)<-c("alpha","beta","c1","c1.est1","c1.est2",
                     "alpha.est","alpha.se","alpha.sd",
                     "beta.est","beta.se","beta.sd",
                     "nm1","nm1.est","nm2","nm2.est","beta.power")

for(i in 1:nrow(grid)){
  c1=grid[i,1]
  c2=c1
  par.est[i,3]=c1
  
  param_est <- foreach( j = 1:num_samples,
                        .combine="cbind",
                        .packages="gumbel") %dopar% {
                          
                          n11<-matrix(NA,nrow=num_rep, ncol= num_wf)
                          n12p<-matrix(NA,nrow=num_rep, ncol= num_wf)
                          n12c<-matrix(NA,nrow=num_rep, ncol= num_wf)
                          n21<-matrix(NA,nrow=num_rep, ncol= num_wf)
                          n22<-matrix(NA,nrow=num_rep, ncol= num_wf)
                          nm<-matrix(NA,nrow=num_rep, ncol= num_wf)
                          
                          C11=NULL
                          
                          for(nr in 1: num_rep){
                            
                            #simulation from gumbel for 2 workflows
                            wf1<- rgumbel(num_sim, alpha=theta1, dim=2, method=1)
                            wf2<- rgumbel(num_sim, alpha=theta2, dim=2, method=1)
                            
                            y1<-wf1[,1]
                            y2<-wf1[,2]
                            num<- nrow(wf1)
                            
                            cdf1 <- ecdf(y1)
                            cdf2 <- ecdf(y2)
                            cdf1.val <- cdf1(y1)*num/(num+1)
                            cdf2.val <- cdf2(y2)*num/(num+1)
                            
                            
                            ## Partitioning data
                            wt1 <- 1*(1*outer(cdf1.val, c1, ">") +1*outer(cdf2.val, c1, ">") ==2) 
                            n11[nr,1] <- length(y1[which(wt1==1)]) #Category 11 (all complete)
                            
                            wt2 <- 1*(1*outer(cdf1.val, c1, ">") +1*outer(cdf2.val, c1, "<=")
                                      +1*outer(cdf2.val, c2, ">") ==3) 
                            n12c[nr,1] <- length(y1[which(wt2==1)]) #Category 12 (all complete)
                            
                            wt3 <- 1*(1*outer(cdf1.val, c1, ">") +1*outer(cdf2.val, c2, "<=") ==2) 
                            n12p[nr,1] <- length(y1[which(wt3==1)]) #Category 12 (Y2 missing)
                            
                            wt4 <- 1*(1*outer(cdf1.val, c1, "<=") +1*outer(cdf2.val, c1, ">") ==2) 
                            n21[nr,1] <- length(y1[which(wt4==1)]) #Category 21 (Y1 missing)
                            
                            wt5 <- 1*(1*outer(cdf1.val, c1, "<=") +1*outer(cdf2.val, c1, "<=") 
                                      +1*outer(cdf2.val, c2, ">") ==3) 
                            n22[nr,1] <- length(y1[which(wt5==1)])#Category 22 (Y1 missing)
                            
                            wt6 <- 1*(1*outer(cdf1.val, c1, "<=") +1*outer(cdf2.val, c2, "<=") ==2)
                            nm[nr,1] <- length(y1[which(wt6==1)]) #Category 22 (Completely missing)
                            
                            w1 <- 1*(1*outer(cdf1.val, c1, ">") ==1) 
                            y1[which(w1==0)]=0
                            w2 <- 1*(1*outer(cdf2.val, c2, ">") ==1) 
                            y2[which(w2==0)]=0
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
                            wt1 <- 1*(1*outer(cdf1.val, c1, ">") +1*outer(cdf2.val, c1, ">") ==2) 
                            n11[nr,2] <- length(y1[which(wt1==1)]) #Category 11 (all complete)
                            
                            wt2 <- 1*(1*outer(cdf1.val, c1, ">") +1*outer(cdf2.val, c1, "<=")
                                      +1*outer(cdf2.val, c2, ">") ==3) 
                            n12c[nr,2] <- length(y1[which(wt2==1)]) #Category 12 (all complete)
                            
                            wt3 <- 1*(1*outer(cdf1.val, c1, ">") +1*outer(cdf2.val, c2, "<=") ==2) 
                            n12p[nr,2] <- length(y1[which(wt3==1)]) #Category 12 (Y2 missing)
                            n12=n12c+n12p
                            
                            wt4 <- 1*(1*outer(cdf1.val, c1, "<=") +1*outer(cdf2.val, c1, ">") ==2) 
                            n21[nr,2] <- length(y1[which(wt4==1)]) #Category 21 (Y1 missing)
                            
                            wt5 <- 1*(1*outer(cdf1.val, c1, "<=") +1*outer(cdf2.val, c1, "<=") 
                                      +1*outer(cdf2.val, c2, ">") ==3) 
                            n22[nr,2] <- length(y1[which(wt5==1)])#Category 22 (Y1 missing)
                            
                            wt6 <- 1*(1*outer(cdf1.val, c1, "<=") +1*outer(cdf2.val, c2, "<=") ==2)
                            nm[nr,2] <- length(y1[which(wt6==1)]) #Category 22 (Completely missing)
                            
                            w1 <- 1*(1*outer(cdf1.val, c1, ">") ==1) 
                            y1[which(w1==0)]=0
                            w2 <- 1*(1*outer(cdf2.val, c2, ">") ==1) 
                            y2[which(w2==0)]=0
                            C11b=as.data.frame(cbind(y1,y2,1,nr))
                            colnames(C11b)=c("y1","y2","x","Rep")
                            
                            C11=as.data.frame(rbind(C11,C11a,C11b))
                            colnames(C11)=c("y1","y2","x","Rep")
                          }
                          
                          n_obs=round(colMeans(n11+n12+n21))
                          c1_est=c(c1,c1)
                          alp_est=1.4
                          bet_est=0.09
                          thet=c(alp_est,bet_est,c1_est)
                          xx=unique(C11$x)
                          e_n22=rep(NA,2)
                          e_n22[1]=n_obs[1]*(gfun.inv((alp_est + (bet_est*xx[1]))* hfun(c1_est[1])))/(1-gfun.inv((alp_est + (bet_est*xx[1]))* hfun(c1_est[1])))
                          e_n22[2]=n_obs[2]*(gfun.inv((alp_est + (bet_est*xx[2]))* hfun(c1_est[2])))/(1-gfun.inv((alp_est + (bet_est*xx[2]))* hfun(c1_est[2])))
                          fit<- optim(thet,fn = function(thets){-likfun(thets)})
                          alp_est=fit$par[1]
                          bet_est=fit$par[2]
                          c1_est=fit$par[3:4]
                          while(likfun(c(alp_est,bet_est,c1_est))-likfun(thet)>0.01){
                            thet=c(alp_est,bet_est,c1_est)
                            e_n22[1]=n_obs[1]*(gfun.inv((alp_est + (bet_est*xx[1]))* hfun(c1_est[1])))/(1-gfun.inv((alp_est + (bet_est*xx[1]))* hfun(c1_est[1])))
                            e_n22[2]=n_obs[2]*(gfun.inv((alp_est + (bet_est*xx[2]))* hfun(c1_est[2])))/(1-gfun.inv((alp_est + (bet_est*xx[2]))* hfun(c1_est[2])))
                            fit<- optim(thet,fn = function(thets){-likfun(thets)})
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
                              
                              n11[nr,]=c(nrow(C11[C11$y1>0 & C11$y2>0 & C11$x==0,]),
                                         nrow(C11[C11$y1>0 & C11$y2>0 & C11$x==1,]))
                              n12[nr,]=c(nrow(C11[C11$y1>0 & C11$y2==0 & C11$x==0,]),
                                         nrow(C11[C11$y1>0 & C11$y2==0 & C11$x==1,]))
                              n21[nr,]=c(nrow(C11[C11$y1==0 & C11$y2>0 & C11$x==0,]),
                                         nrow(C11[C11$y1==0 & C11$y2>0 & C11$x==1,]))
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
                              e_n22[1]=n_obs[1]*(gfun.inv((alp_est + (bet_est*xx[1]))* hfun(c1_est[1])))/(1-gfun.inv((alp_est + (bet_est*xx[1]))* hfun(c1_est[1])))
                              e_n22[2]=n_obs[2]*(gfun.inv((alp_est + (bet_est*xx[2]))* hfun(c1_est[2])))/(1-gfun.inv((alp_est + (bet_est*xx[2]))* hfun(c1_est[2])))
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

out=round(par.est, digits=6)
out
write.csv(out, file="Gumbel_theta=2,1.7_Rep3_grid1.csv")

