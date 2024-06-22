## Script for comparing Transplex, Smarter and C1 using mdCCR ##
rm(list=ls())

library(Matrix)
library(gumbel)
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

prepare.Uim <- function(a, tm, c0){
  
  a.cdf1 <- ecdf(a[a$y1>0,1])
  a.cdf2 <- ecdf(a[a$y2>0,2])
  cdf1.value <- (a.cdf1(a[,1])*(1-c0[1])) + c0[1]
  cdf2.value <- (a.cdf2(a[,2])*(1-c0[2])) + c0[2]
  
  ## n12c
  wt1 <- 1*(1*outer(cdf1.value, c0[1], ">") +1*outer(cdf2.value, c0[2], ">") +
              1*outer(cdf2.value, c0[1], "<=")==3) 
  
  ## counts in M categories (n111, n112, ... n11M)
  w <- 1*(1*outer(cdf1.value, c0[1], ">") +1*outer(cdf2.value, c0[1], ">") ==2) 
  cdf1.value <- cdf1.value[which(w==1)]
  cdf2.value <- cdf2.value[which(w==1)]
  
  wt.temp <- 1*(1*outer(cdf1.value, tm, ">") + 1*outer(cdf2.value, tm, ">")==2) 
  wt <- cbind(wt.temp[, -m]- wt.temp[, -1] , wt.temp[, m])
  
  return(c(colSums(wt),sum(wt1)))
}

L11 <- function(coef,tm){
  
  m <- length(tm)
  ht <- hfun(tm)
  
  ght <- 1 - (2*tm) + gfun.inv(coef * ht) # Psi(tm)
  dght <- c(ght[-m]-ght[-1], ght[m])
  
  lpr <- ifelse(dght > 0, log(dght), Inf) 
  
  return(lpr) #log P11
}

L12c<-function(coef,c0){
  fn= c0[1] - c0[2] + pgumbel(u=c0[1],v=c0[2],alpha=log(2)/log(coef),dim=2) - gfun.inv(coef * hfun(c0[1])) 
  lpr <- ifelse(fn > 0, log(fn), Inf)
  return(lpr) #log P12c
}

L12p<-function(coef,c0){
  fn= c0[2] - pgumbel(u=c0[1],v=c0[2],alpha=log(2)/log(coef),dim=2)
  lpr <- ifelse(fn > 0, log(fn), Inf)
  return(lpr) #log P12p
}

L21p<-function(coef,c0){
  fn= c0[1] - pgumbel(u=c0[1],v=c0[2],alpha=log(2)/log(coef),dim=2)
  lpr <- ifelse(fn > 0, log(fn), Inf)
  return(lpr) #log P21p
}

L22<-function(coef,c0){
  fn= pgumbel(u=c0[1],v=c0[2],alpha=log(2)/log(coef),dim=2)
  lpr <- ifelse(fn > 0, log(fn), Inf)
  return(lpr) #log 1- P(phi)
}

likfun<-function(thets){ 
  
  if(thets[1]<1 | thets[1]>2 | any(thets[2:num_wf]< -1) | any(thets[2:num_wf]>1) | 
     any(thets[-c(1:num_wf)]<0) | any(thets[-c(1:num_wf)]>1) | any(thets[(num_wf+1):(2*num_wf)] < thets[(2*num_wf+1):(3*num_wf)]) ){ # avoiding negative values for nm
    return(Inf)
  }
  
  r.unique <- levels(factor(C11$Rep)) # no. of replicates per wf
  x.unique <- levels(factor(C11$x)) # no. of workflows
  
  cc<- thets[(num_wf+1):(2*num_wf)]
  cc2<- thets[(2*num_wf+1):(3*num_wf)]
  
  l11<-0
  l12p<-0
  l12c<-0
  l21p<-0
  l22<-0
  #const<-0
  
  for(k in 1:length(x.unique)){
    
    if(cc[k]>0){
      tm <-seq(cc[k], .99, length.out=m)
    }else{
      tm <-seq(0.0001, .99, length.out=m)
    }
    
    z<- as.numeric(x.unique[k])
    x1 <- 1*(x.unique[k] == "1")
    x2 <- 1*(x.unique[k] == "2")
    
    coef<-thets[1] + (thets[2]*x1) + (thets[3]*x2)
    
    for(l in 1:length(r.unique)){
      
      lik0<-sum(L11(coef,tm)*n11[l,k,])
      l11<-l11+lik0
      
      lik1<-ifelse(n12c[l,k]>0,(L12c(coef,c(cc[k],cc2[k])))*n12c[l,k],0)
      l12c= l12c+lik1
      
      lik2<-ifelse(n12p[l,k]>0,(L12p(coef,c(cc[k],cc2[k])))*n12p[l,k],0)
      l12p= l12p+lik2
      
      lik3<-ifelse(n21p[l,k]>0,(L21p(coef,c(cc[k],cc2[k])))*n21p[l,k],0)
      l21p= l21p+lik3
      
      lik4<-ifelse(e_n22[k]>0,(L22(coef,c(cc[k],cc2[k])))*(e_n22[k]),0)
      l22= l22+lik4
      
    }
  }
  return(l11+l12c+l12p+l21p+l22) #log likelihood
}


## reading data
{
  dat1<- as.data.frame(read.csv("Transplex1.csv")[,10])
  dat2<- as.data.frame(read.csv("Transplex2.csv")[,10])
  dat3<- as.data.frame(read.csv("Transplex3.csv")[,10])
  dat4<- as.data.frame(read.csv("Smarter1.csv")[,10])
  dat5<- as.data.frame(read.csv("Smarter2.csv")[,10])
  dat6<- as.data.frame(read.csv("Smarter3.csv")[,10])
  dat7<- read.table(file = '3801.tsv', sep = '\t', header = TRUE)[,2]
  dat8<- read.table(file = '3831.tsv', sep = '\t', header = TRUE)[,2]
  dat9<- read.table(file = '3826.tsv', sep = '\t', header = TRUE)[,2]
  dat10<- read.table(file = '3866.tsv', sep = '\t', header = TRUE)[,2]
  dat11<- read.table(file = '3836.tsv', sep = '\t', header = TRUE)[,2]
  dat12<- read.table(file = '3851.tsv', sep = '\t', header = TRUE)[,2]
  
  C11=as.data.frame(rbind(cbind(rank(dat2),rank(dat1),0,1),cbind(rank(dat3),rank(dat2),0,2),
                          cbind(rank(dat3),rank(dat1),0,3),
                          cbind(rank(dat5),rank(dat4),1,1),cbind(rank(dat5),rank(dat6),1,2),
                          cbind(rank(dat6),rank(dat4),1,3),
                          cbind(rank(dat7),rank(dat8),2,1),cbind(rank(dat9),rank(dat10),2,2),
                          cbind(rank(dat12),rank(dat11),2,3)))
}
## 0 is Transplex; 1 is Smarter; 2 is C1 
colnames(C11)=c("y1","y2","x","Rep")

rscale<- function(x){(x-min(x))/(max(x)-min(x))}

num_rep=3
num_wf=3

## scaling the rank data
for(i in 0:(num_wf-1)){
  for(j in 1:num_rep){
    C11[C11$x==i & C11$Rep==j,1]=rscale(C11[C11$x==i & C11$Rep==j,1])
    C11[C11$x==i & C11$Rep==j,2]=rscale(C11[C11$x==i & C11$Rep==j,2])
  }
}

n1=matrix(NA,nrow=num_rep, ncol=num_wf)
n12p=matrix(NA,nrow=num_rep, ncol=num_wf)
n21p=matrix(NA,nrow=num_rep, ncol=num_wf)
nm=matrix(NA,nrow=num_rep, ncol=num_wf)

## Number of values in each block
for(i in 1:num_wf){
  for(j in 1:num_rep){
    n1[j,i]=nrow(C11[C11$y1 >0 & C11$y2 >0 & C11$x==(i-1) & C11$Rep==j,]) # complete
    n12p[j,i]=nrow(C11[C11$y1 >0 & C11$y2 ==0 & C11$x==(i-1) & C11$Rep==j,]) #partial
    n21p[j,i]=nrow(C11[C11$y1 == 0 & C11$y2 > 0 & C11$x==(i-1) & C11$Rep==j,]) #partial
    nm[j,i]<- nrow(C11[C11$y1 == 0 & C11$y2 == 0 & C11$x==(i-1) & C11$Rep==j,]) # missing
  }
}

m=30
n_obs=colMeans(n1+n12p+n21p)

c1_est<-colMeans((n21p+nm)/(n1+n12p+n21p+nm))
c2_est<-colMeans((n12p+nm)/(n1+n12p+n21p+nm))
alp_est=1.3
bet_est1=-0.1
bet_est2=-0.3

n12c=matrix(NA,nrow=num_rep, ncol= num_wf)
n11=array(NA, c(num_rep, num_wf, m))
e_n22=rep(NA,num_wf)


r.unique <- levels(factor(C11$Rep)) # no. of replicates per wf
x.unique <- levels(factor(C11$x)) # no. of workflows

for(k in 1:num_wf){
  z<- as.numeric(x.unique[k])
  x1 <- 1*(x.unique[k] == "1")
  x2 <- 1*(x.unique[k] == "2")
  
  coeff<-alp_est + (bet_est1*x1) + (bet_est2*x2)
  e_n22[k]=n_obs[k]*pgumbel(u=c1_est[k],v=c2_est[k],alpha=log(2)/log(coeff),dim=2)/(1-pgumbel(u=c1_est[k],v=c2_est[k],alpha=log(2)/log(coeff),dim=2))
  
  for(l in 1:num_rep){
    
    r<-as.numeric(r.unique[l])
    
    wt <- prepare.Uim(C11[C11$x==z & C11$Rep==r,1:2], seq(c1_est[k], .99, length.out=m), c(c1_est[k],c2_est[k]))
    n11[l,k,]=wt[1:m]
    n12c[l,k]=wt[m+1]
  }
}

thet=c(alp_est,bet_est1,bet_est2,c1_est,c2_est)
lik=likfun(c(alp_est,bet_est1,bet_est2,c1_est,c2_est))-sum(lfactorial(e_n22))+
  sum(lfactorial(n_obs+e_n22-1))-sum(lfactorial(n12c))-sum(lfactorial(n11))

print(round(c(alp_est,bet_est1,bet_est2,c1_est,c2_est),3))

fit<- optim(thet,fn = function(thets){-likfun(thets)})

lik=c(lik,likfun(fit$par)-sum(lfactorial(e_n22))+
        sum(lfactorial(n_obs+e_n22-1))-sum(lfactorial(n12c))-sum(lfactorial(n11)))
itr=2
while((lik[itr]-lik[itr-1])>0.01){
  ## updating parameters
  alp_est=fit$par[1]
  bet_est1=fit$par[2]
  bet_est2=fit$par[3]
  c1_est=fit$par[(num_wf+1):(2*num_wf)]
  c2_est=fit$par[(2*num_wf+1):(3*num_wf)]
  thet=c(alp_est,bet_est1,bet_est2,c1_est,c2_est)
  
  print(round(c(alp_est,bet_est1,bet_est2,c1_est,c2_est),3))
  
  ## E-step
  for(k in 1:num_wf){
    z<- as.numeric(x.unique[k])
    x1 <- 1*(x.unique[k] == "1")
    x2 <- 1*(x.unique[k] == "2")
    
    coeff<-alp_est + (bet_est1*x1) + (bet_est2*x2)
    e_n22[k]=n_obs[k]*pgumbel(u=c1_est[k],v=c2_est[k],alpha=log(2)/log(coeff),dim=2)/(1-pgumbel(u=c1_est[k],v=c2_est[k],alpha=log(2)/log(coeff),dim=2))
    
    for(l in 1:num_rep){
      
      r<-as.numeric(r.unique[l])
      
      wt <- prepare.Uim(C11[C11$x==z & C11$Rep==r,1:2], seq(c1_est[k], .99, length.out=m), c(c1_est[k],c2_est[k]))
      n11[l,k,]=wt[1:m]
      n12c[l,k]=wt[m+1]
    }
  }
  
  ## M-step
  fit<- optim(thet,fn = function(thets){-likfun(thets)})
  lik=c(lik,likfun(fit$par)-sum(lfactorial(e_n22))+
          sum(lfactorial(n_obs+e_n22-1))-sum(lfactorial(n12c))-sum(lfactorial(n11)))
  itr=itr+1
  
}

est=c(alp_est,bet_est1,bet_est2)
nm.est=e_n22
res=as.data.frame(cbind(est,c1_est,c2_est,nm.est))
res
colnames(res)=c("Coefficient","c1.est","c2.est","nm.est")
rownames(res)=c("TransPlex","SMARTer","C1")
write.csv(res,"Platform_results_mdCCR.csv")
write.csv(lik,"Platform_likelihood_mdCCR.csv")
