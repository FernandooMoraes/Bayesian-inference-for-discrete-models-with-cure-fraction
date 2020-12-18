library(TeachingDemos)
library(splines)
library(survival)
library(MCMCpack)
require(coda)
require(lattice)
require(MASS)
set.seed(2020)


time=bladder[,5]
censure=bladder[,6]
n=length(censure)


######## WEIBULL DISCRETE SURVIVAL FUNCTION WITH CURED FRACTION
surv.weidc.fc<-function(x,q,b,f){ f + (1-f)*( (q)^((x+1)^b) )}

# EMPIRICAL K-M ESTIMATOR

data<-Surv(time,censure)
km<-survfit(data~1)
plot(km,conf.int=F,ylim=c(0.1,1))

# Likelihood function

VERO<-function(p,time,censure){
  q<-p[1]
  b<-p[2]
  f<-p[3]
  if ( (q>0) && (q<1) && (b>0) && (f>0) && (f<1) )
    return (-1*(
      sum(censure*log(1-f))+sum(censure*log(q^(time^b)-q^((time+1)^b)))
      +sum((1-censure)*log(f+(1-f)*q^((time+1)^b)))))
  else return (-Inf)
}

a<-optim(c(0.9,1,0.15),VERO,time=time,censure=censure)

q.est<-a$par[1]
b.est<-a$par[2]
f.est<-a$par[3]

xx<-sort(time)
survWei2<-surv.weidc.fc(xx,q.est,b.est,f.est)
points(xx,survWei2,type="b",col=4)
posterior<-function(p,time,censure){
  q<-p[1]
  b<-p[2]
  f<-p[3]
  if ( (q>0) && (q<1) &&  (f>0) && (f<1) )
    return ((
      sum( censure*log(1-f) )
      +sum( censure * log( q^(time^b) - q^((time+1)^b) ) )
      +sum( (1-censure)*log(f + (1-f)* q^((time+1)^b) ) )
      + dbeta(q,10^-3,10^-3,log=T)
      + dgamma(f,1,1,log=T)
    ))
  else return (-Inf)
}
M<-100000
theta<-MCMCmetrop1R(posterior,c(q.est,b.est,f.est),logfun=TRUE,burnin=1000,mcmc=M,time=time,censure=censure)
q.est2<-mean(theta[,1])
b.est2<-mean(theta[,2])
f.est2<-mean(theta[,3])

# Graphics

xx<-sort(time)
survWei3<-surv.weidc.fc(xx,q.est2,b.est2,f.est2)
points(xx,survWei3,type="b",col=2)

# FBST H1:b=2 
post.mcmc<-numeric(M)
for (j in 1:M) { post.mcmc[j]<-posterior(theta[j,],time,censure) }

posterior.H02<-function(p,time,censure){
  q<-p[1]
  b<-2
  f<-p[2]
  if ( (q>0) && (q<1) && (f>0) && (f<1) )
    return (-1*(
      sum( censure*log(1-f) )
      +sum( censure * log( q^(time^b) - q^((time+1)^b) ) )
      +sum( (1-censure)*log(f + (1-f)* q^((time+1)^b) ) )
      + dbeta(q,10^-3,10^-3,log=T)
      + dgamma(f,1,1,log=T)
    ))
  else return (-Inf)
}
a.H02<-optim(c(q.est,f.est),posterior.H02,time=time,censure=censure)

# Max of posterior surv H02

max.post.H02<- -1*a.H02$value
e.value2<-sum(post.mcmc<max.post.H02)/M


# FBST H2:b=1 
posterior.H03<-function(p,time,censure){
  q<-p[1]
  b<-1
  f<-p[2]
  if ( (q>0) && (q<1) && (f>0) && (f<1) )
    return (-1*(
      sum( censure*log(1-f) )
      +sum( censure * log( q^(time^b) - q^((time+1)^b) ) )
      +sum( (1-censure)*log(f + (1-f)* q^((time+1)^b) ) )
      + dbeta(q,10^-3,10^-3,log=T)
      + dgamma(f,1,1,log=T)
    ))
  else return (-Inf)
}
a.H03<-optim(c(q.est,f.est),posterior.H03,time=time,censure=censure)

# Max of posterior surv H02

max.post.H03<- -1*a.H03$value
e.value3<-sum(post.mcmc<max.post.H03)/M


#Graphics of estimatives

data<-Surv(time,censure)
km<-survfit(data~1)
plot(km,conf.int=F,xlab="t",ylab="S(t)")
### q
q.full<-mean(theta[,1])
q.H02<-a.H02$par[1]
q.H03<-a.H03$par[1]

### b
b.full<-mean(theta[,2])
b.H02<-2
b.H03<-1

### f
f.full<-mean(theta[,3])
f.H02<-a.H02$par[2]
f.H03<-a.H03$par[2]

xxx<-seq(0,max(time))
surv.full<-surv.weidc.fc(xxx,q.full,b.full,f.full)
surv.H02<-surv.weidc.fc(xxx,q.H02,b.H02,f.H02)
surv.H03<-surv.weidc.fc(xxx,q.H03,b.H03,f.H03)


pdf('ap_2.pdf', width = 12, height = 4)

par(mfrow=c(1,3))
plot(km,conf.int=F,xlab="t",ylab="S(t)",main=c("Complete Model","Discrete Weibull "),ylim=c(0.4,1))
points(xxx,surv.full,type="b",col=1,pch=16)
plot(km,conf.int=F,xlab="t",ylab="S(t)",main=c("Model under H01","Discrete Rayleigh Exp"),ylim=c(0.4,1))
points(xxx,surv.H02,type="b",col=2,pch=16)
plot(km,conf.int=F,xlab="t",ylab="S(t)",main=c("Model under H02","Discrete Exp Geomeometry"),ylim=c(0.4,1))
points(xxx,surv.H03,type="b",col=3,pch=16)

dev.off()

#e- values

e.value2;e.value3

##### Estimatives And HPD #####

#q

q.est<-mean(theta[,1]);q.est
emp.hpd(theta[,1])

# b

b.est<-mean(theta[,2]);b.est
emp.hpd(theta[,2])

# f

f.est<-mean(theta[,3]);f.est
emp.hpd(theta[,3])  

