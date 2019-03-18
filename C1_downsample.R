## Finding sequencing depth using downsampled data
rm(list=ls())
library(Matrix)
#setwd("~/Box Sync/Qunhua/Data")

#dat1<- as.data.frame(read.csv("C1_1.csv")[,-c(2,3,4,5,8,9,13)])
#dat2<- as.data.frame(read.csv("C1_4.csv")[,-c(2,3,4,5,8,9,13)])
dat1<- read.csv("C1_100.csv")[,-c(1,2)]
dat3<- read.csv("C1_10.downsample.csv")[,-c(1,2)]
dat4<- read.csv("C1_20.downsample.csv")[,-c(1,2)]
dat5<- read.csv("C1_30.downsample.csv")[,-c(1,2)]
dat6<- read.csv("C1_40.downsample.csv")[,-c(1,2)]
dat7<- read.csv("C1_50.downsample.csv")[,-c(1,2)]
#dat8<- read.csv("C1_60.downsample.csv")[,-c(1,2)]
dat9<- read.csv("C1_75.downsample.csv")[,-c(1,2)]

C11=as.data.frame(rbind(cbind(rank(dat1[,2]),rank(dat1[,1]),0),cbind(rank(dat3[,2]),rank(dat3[,1]),1),
                        cbind(rank(dat4[,2]),rank(dat4[,1]),2),cbind(rank(dat5[,2]),rank(dat5[,1]),3),
                        cbind(rank(dat6[,2]),rank(dat6[,1]),4),cbind(rank(dat7[,2]),rank(dat7[,1]),5),
                        #cbind(rank(dat8[,2]),rank(dat8[,1]),5),
                        cbind(rank(dat9[,2]),rank(dat9[,1]),6)))
#C11=as.data.frame(rbind(cbind(rank(dat2[,4]),rank(dat1[,4]),0),cbind(rank(dat3[,2]),rank(dat3[,1]),1),
 #                       cbind(rank(dat4[,2]),rank(dat4[,1]),2),cbind(rank(dat5[,2]),rank(dat5[,1]),3),
  #                      cbind(rank(dat6[,2]),rank(dat6[,1]),4),cbind(rank(dat7[,2]),rank(dat7[,1]),5),
   #                     cbind(rank(dat8[,2]),rank(dat8[,1]),6)))
## 0 is Original; 1 is 10%; 2 is 20%; 3 is 30%; 4 is 40%, 5 is 50%; 6 is 60% and 7 is 75%
colnames(C11)=c("y1","y2","x")
s=7
## scaling the rank data
for(i in 0:(s-1)){
C11[C11$x==i,1]=(C11[C11$x==i,1]-min(C11[C11$x==i,1]))/(max(C11[C11$x==i,1])- min(C11[C11$x==i,1]))
C11[C11$x==i,2]=(C11[C11$x==i,2]-min(C11[C11$x==i,2]))/(max(C11[C11$x==i,2])- min(C11[C11$x==i,2]))
}

#plot(C11[C11$x==4,1],C11[C11$x==4,2],main="50% downsampled data",xlab="y1",ylab="y2", pch=20)
plot(C11[C11$x==0,1],C11[C11$x==0,2],main="1.2M",xlab="y1",ylab="y2", pch=20)

n1=rep(NA,s)
n12p=rep(NA,s)
n21p=rep(NA,s)
nm=rep(NA,s)

## Number of values in each block
for(i in 0:s){
n1[i]=nrow(C11[C11$y1 >0 & C11$y2 >0 & C11$x==i-1,]) # complete
n12p[i]=nrow(C11[C11$y1 >0 & C11$y2 ==0 & C11$x==i-1,]) #partial
n21p[i]=nrow(C11[C11$y1 == 0 & C11$y2 > 0 & C11$x==i-1,]) #partial
nm[i]<- nrow(C11[C11$y1 == 0 & C11$y2 == 0 & C11$x==i-1,]) # missing
}

## Functions for Gumbel-Hougaard
gfun <- function(x){log(x)}
gfun.inv <- function(x){exp(x)}
hfun <- function(x){log(x)}

## Functions for log likelihood
prepare.Uim <- function(a, tm, c0){ ##Computing Uim without adjusting the cdf
  
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

L11 <- function(a,coef,tm, c0){
  
  m <- length(tm)
  
  wt <- prepare.Uim(a, tm, c0) # getting Uim
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


L12<-function(coef,c0){
  fn<- c0 - gfun.inv(coef * hfun(c0))
  lpr <- ifelse(fn > 0, log(fn), Inf)
  return(lpr) #log L12 or L21
}

L22<-function(coef,c0){
  fn<- gfun.inv(coef * hfun(c0))
  lpr <- ifelse(fn > 0, log(fn), Inf)
  return(lpr) #log L22
}

likfun<-function(thets){ #thets= (alpha, beta1,..., nm0,nm1,...)
  
  if(thets[1]<1 | thets[1]>2 | any(thets[2:s]< -1) | any(thets[2:s]>1) | any(thets[(s+1):(2*s)]<0) ){
    # avoiding negative values for nm
    return(Inf)
  }
  
  cc<-rep(NA,s)
  cc2<-rep(NA,s)
  n11=rep(NA,s)
  n21=rep(NA,s)
  n12=rep(NA,s)
  n22=rep(NA,s)
  
  x.unique <- levels(factor(C11$x)) # getting types of workflows
  l11<-0
  l12<-0
  l21<-0
  l22<-0
  const<-0
  
  for(k in 1:length(x.unique)){
    
    # computing c1 for a particular workflow
    cc[k]<-(n21p[k]+thets[k+s])/(n1[k]+n12p[k]+n21p[k]+thets[k+s])
    cc2[k]<-(n12p[k]+thets[k+s])/(n1[k]+n12p[k]+n21p[k]+thets[k+s])
    
    if(cc[k]>0){
      tm <-seq(cc[k], .99, length.out=m)
    }else{
      tm <-seq(0.0001, .99, length.out=m)
    }
    
    z<- as.numeric(x.unique[k])
    x1 <- 1*(x.unique[k] == "1")
    x2 <- 1*(x.unique[k] == "2")
    x3 <- 1*(x.unique[k] == "3")
    x4 <- 1*(x.unique[k] == "4")
    x5 <- 1*(x.unique[k] == "5")
    x6 <- 1*(x.unique[k] == "6")
    #x7 <- 1*(x.unique[k] == "7")
    
    coef<-thets[1] + (thets[2]*x1) + (thets[3]*x2) + (thets[4]*x3) + (thets[5]*x4) + (thets[6]*x5) + 
      (thets[7]*x6) #+ (thets[8]*x7)
    
    cdf1 <- ecdf(C11[C11$y1>0 & C11$x==z,1])
    cdf2 <- ecdf(C11[C11$y2>0 & C11$x==z,2])
    cdf1.val <- (cdf1(C11[C11$x==z,1])*(1-cc[k])) + cc[k]
    cdf2.val <- (cdf2(C11[C11$x==z,2])*(1-cc2[k])) + cc2[k]
    
    wt <- 1*(1*outer(cdf1.val, cc[k], ">") +1*outer(cdf2.val, cc[k], ">") ==2) 
    n11[k]=sum(wt)
    
    n12[k]=n12p[k]+ n1[k] -n11[k]
    n22[k]=(cc[k]*(n1[k]+n12p[k]+n21p[k]+thets[k+s])) - (n12[k]+thets[k+s])
    n21[k]=n21p[k]-n22[k]
    
    lik0<-L11(C11[C11$x==z,1:2], coef, tm, c(cc[k],cc2[k]))
    l11<-l11+lik0
    
    lik1<-ifelse(n12[k]>0,(L12(coef,cc[k]))*n12[k],0)
    l12= l12+lik1
    
    lik2<-ifelse(n21[k]>0,(L12(coef,cc[k]))*n21[k],0)
    l21= l21+lik2
    
    lik3<-ifelse(n22[k]+thets[k+s]>0,(L22(coef,cc[k]))*(n22[k]+thets[k+s]),0)
    l22= l22+lik3
    
    const= const+ lfactorial(n11[k]+n12[k]+n21[k]+n22[k]+thets[k+s]) - lfactorial(n22[k]+thets[k+s])
  }
  return(const+l11+l12+l21+l22) #log (L11 L12 L21 L22)
}

m=30
fit<- optim(c(1.7,0,0,0,0,0,0,nm),
            fn = function(thets){-likfun(thets)}, hessian=TRUE)
fit$par

est <- fit$par        # estimate
#SE=nearPD(solve(fit$hessian))
#ese <- sqrt(diag(SE$mat))       # standard error
ese <- sqrt(diag(solve(fit$hessian)))

CI.low <- est[1:s] - qnorm(1-.05/2, 0, 1) * ese[1:s]
CI.up <- est[1:s] + qnorm(1-.05/2, 0, 1) * ese[1:s]
res=as.data.frame(cbind(est[1:s],ese[1:s],CI.low,CI.up))
coeff=c("alpha","beta1","beta2","beta3","beta4","beta5","beta6")#,"beta7")
res=cbind(coeff,round(res,5),nm,round(fit$par[(s+1):(2*s)],5))
colnames(res)=c("Coefficient","Estimate","SE","CI.low","CI.up","True.nm","nm.estimate")
rownames(res)=c("1.2M","0.12M","0.24M","0.36M","0.48M","0.60M","0.90M")
write.csv(res,"C1.csv")



## Plots
library(ggplot2)
df=cbind(c(0.12,0.24,0.36,0.48,0.60,0.90),res[-1,])
colnames(df)=c("seqd",colnames(res))
ggplot(df, aes(x = seqd, y = Estimate)) + ylab("Estimate of coefficient") + ylim(-0.03,0.1) +
  geom_point(size = 2) + geom_hline(yintercept=0, linetype="dashed") + xlab("Sequencing Depth") +
  geom_errorbar(aes(ymax = CI.up, ymin = CI.low), size=0.5, width=0.02) +
  scale_x_continuous(breaks=c(0.12,0.24,0.36,0.48,0.60,0.90), 
                     label=c("0.12M","0.24M","0.36M","0.48M","0.60M","0.90M"))


ggplot(df, aes(x = rownames(df), y = nm.estimate)) + ylab("Number of missing values") + 
  geom_point(aes(x = rownames(df), y = nm.estimate), size = 2) + 
  geom_point(aes(x = rownames(df), y = True.nm), shape=2, size = 2) + xlab("Sequencing Depth")

ggplot(df, aes(x = True.nm, y = nm.estimate)) + ylab("Estimate of missing values") + ylim(15800,18700) + 
  xlim(15800,18700) + geom_point(size = 2) +geom_text(aes(label=rownames(df)),hjust=0.8, vjust=-1) + 
  xlab("True number of missing values") + geom_abline(intercept=0, slope=1, linetype="dashed")
  

x.unique <- levels(factor(C11$x))
c1<-rep(NA,s)
c2<-rep(NA,s)
df=NULL

for(k in 1:s){
  c1[k]<-(n21p[k]+est[k+s])/(n1[k]+n12p[k]+n21p[k]+est[k+s])
  c2[k]<-(n12p[k]+est[k+s])/(n1[k]+n12p[k]+n21p[k]+est[k+s])
  tm <-seq(c1[k], .99, length.out=m)
  
  z<- as.numeric(x.unique[k])
  x1 <- 1*(x.unique[k] == "1")
  x2 <- 1*(x.unique[k] == "2")
  x3 <- 1*(x.unique[k] == "3")
  x4 <- 1*(x.unique[k] == "4")
  
  coef<-est[1] + (est[2]*x1) + (est[3]*x2) + (est[4]*x3) + (est[5]*x4)
  a=C11[C11$x==z,1:2]
  wt <- prepare.Uim(a, tm, c(c1[k],c2[k])) # getting Uim
  ht <- hfun(tm)
  lik<- log(gfun.inv(coef * ht))
  
  df=cbind(df,log(tm),lik)
  
  #ght <- 1 - (2*tm) + gfun.inv(coef * ht) #psi*(tm)
  #lik1<- ifelse(ght > 0, log(ght), Inf)
  
  #dght <- c(ght[-m]-ght[-1], ght[m]) # Psi*(tm-1) - Psi*(tm)
  #lik2 <- ifelse(dght > 0, log(dght), Inf) 
  
}

df=as.data.frame(df)
colnames(df)=c("t1","lik1","t2","lik2","t3","lik3","t4","lik4","t5","lik5")
ggplot(df, aes(x = t1, y = lik1)) + ylab("log Psi(t)") + xlab("log t") + 
  geom_line(aes(x = t1, y = lik1)) + geom_line(aes(x = t2, y = lik2),linetype=2) + 
  geom_line(aes(x = t3, y = lik3),linetype=3) + geom_line(aes(x = t4, y = lik4),linetype=4) + 
  geom_line(aes(x = t5, y = lik5),linetype=5)


