#' Four scenarios for simulating dose interval data
#'
#' there are 4 scenarios: linear one-sided, non-linear one-sided, linear two-sided and non-linear two-sided. Each scenario has 1 continuous treatment case and an ordinal treatment case.
#'
#' @param size the number of observations to be drawn
#' @param ncov the number of covariates to be generated
#' @param seed random seed
#'
#' @examples
#'


#' @export
Scenario1.continuous <- function(size,ncov,seed,aL=-1,aU=1){
  set.seed(seed)
  X = matrix(runif(size*ncov,-1,1),ncol=ncov)
  #A = runif(size,aL,aU)
  A = rtruncnorm(size,a=aL,b=aU,mean=mean(aL+aU),sd=0.5)
  propensity=1/dnorm(A,mean(aL+aU),0.5)
  c =  2*(0.15*X[,1]+0.15*X[,2]+0.15*X[,3]+0.15*X[,4])
  r=10
  mu =(0*exp(c*r)+5*exp(r*A))/(exp(r*c)+exp(r*A))
  Y = rnorm(length(mu),mu,1)#+0.15*X[,5]# further noise
  S=2.5
  label=mu>S
  plot(A,Y)
  abline(h=S)
  alpha=matrix(rep(c(0.5,0.5),size),ncol=2,byrow=TRUE)
  traininfo = list(X=X,A=A,Y=Y,mu=mu,S=S,alpha=alpha,propensity=propensity/mean(propensity),opt=c,label=label)
  return(traininfo)
}
#' @export
Scenario1.ordinal <- function(size,ncov,seed,K=10){
  set.seed(seed)
  X = matrix(runif(size*ncov,-1,1),ncol=ncov)
  c =(0.15*X[,1]+0.15*X[,2]+0.15*X[,3]+0.15*X[,4])
  c= pmin(aU,pmax(aL,c))
  levels=seq(aL-0.0001,aU+0.0001,length.out = K+1)
  c=as.numeric(cut(c,breaks=levels,labels=F))
  A = sample(1:K,size,replace = T)
  r=0.5
  mu =(0*exp(c*r)+5*exp(r*A))/(exp(r*c)+exp(r*A))
  Y = rnorm(length(mu),mu,1)+0.15*X[,5]
  S=2.5
  label=mu>S
  plot(A,Y)
  abline(h=S)
  alpha=matrix(rep(c(0.5,0.5),size),ncol=2,byrow=TRUE)
  traininfo = list(X=X,A=A,Y=Y,mu=mu,S=S,alpha=alpha,propensity=rep(1,size),opt=c)
  return(traininfo)
}
#' @export
Scenario2.continuous <- function(size,ncov,seed,aL=-1,aU=1){
  set.seed(seed)
  X = matrix(runif(size*ncov,-1,1),ncol=ncov)
  #A = runif(size,aL,aU)
  A = rtruncnorm(size,a=aL,b=aU,mean=mean(aL+aU),sd=0.5)
  propensity=1/dnorm(A,mean(aL+aU),0.5)
  c =  2*(0.15*I(X[,1]>0.5)-0.15*I(X[,2]<(-0.5))+0.15*sin(pi*X[,3])-0.15*X[,4]^3+0.02)
  r=10
  mu =(0*exp(c*r)+5*exp(r*A))/(exp(r*c)+exp(r*A))
  Y = rnorm(length(mu),mu,1)#+0.15*log(abs(X[,5])+1)-0.058
  S=2.5
  label=mu>S
  plot(A,Y)
  abline(h=S)
  alpha=matrix(rep(c(0.5,0.5),size),ncol=2,byrow=TRUE)
  traininfo = list(X=X,A=A,Y=Y,mu=mu,S=S,alpha=alpha,propensity=propensity/mean(propensity),opt=c,label=label)
  return(traininfo)
}
#' @export
Scenario2.ordinal <- function(size,ncov,seed,K=10){
  set.seed(seed)
  X = matrix(runif(size*ncov,-1,1),ncol=ncov)
  levels=seq(aL-0.0001,aU+0.0001,length.out = K+1)
  c =  (0.15*I(X[,1]>0.5)-0.15*I(X[,2]<(-0.5))+0.15*sin(pi*X[,3])-0.15*X[,4]^3+0.02)
  c=as.numeric(cut(c,breaks=levels,labels=F))
  A = sample(1:K,size,replace = T)
  r=0.5
  mu =(0*exp(c*r)+5*exp(r*A))/(exp(r*c)+exp(r*A))
  Y = rnorm(length(mu),mu,1)+0.15*log(abs(X[,5])+1)-0.058
  S=2.5
  label=mu>S
  plot(A,Y)
  abline(h=S)
  alpha=matrix(rep(c(0.5,0.5),size),ncol=2,byrow=TRUE)
  traininfo = list(X=X,A=A,Y=Y,mu=mu,S=S,alpha=alpha,propensity=rep(1,size),opt=c,label=label)
  return(traininfo)
}
#' @export
Scenario3.continuous <- function(size,ncov,seed,aL=-1,aU=1){
  set.seed(seed)
  X = matrix(runif(size*ncov,-1,1),ncol=ncov)
  c=2*(0.15*X[,1]+0.15*X[,2]+0.15*X[,3]+0.15*X[,4])
  A = runif(size,aL,aU)
  #mu =-2*(A-c)^2*(abs(A-c)<0.5)-(abs(A-c)>=0.5)*0.5
  mu =(1-abs(c-A))
  Y = rnorm(length(mu),mu,0.1)#+0.15*X[,5]
  S=0.5
  len=0.5
  label=mu>S
  plot(A,Y)
  abline(h=S)
  alpha=matrix(rep(c(0.5,0.5),size),ncol=2,byrow=TRUE)
  traininfo = list(X=X,A=A,Y=Y,mu=mu,S=S,alpha=alpha,propensity=rep(1,size),opt_L=pmax(c-len,aL),opt_R=pmin(c+len,aU),label=label)
  return(traininfo)
}
#' @export
Scenario3.ordinal <- function(size,ncov,seed,K=20){
  set.seed(seed)
  levels=seq(aL-0.0001,aU+0.0001,length.out = K+1)
  X = matrix(runif(size*ncov,-1,1),ncol=ncov)
  c=(0.15*X[,1]+0.15*X[,2]+0.15*X[,3]+0.15*X[,4])
  c=pmax(pmin(aU,c),aL)
  c=as.numeric(cut(c,breaks=levels,labels=F))
  A = sample(1:K,size,replace = T)
  mu =(K-abs(c-A))
  Y = rnorm(length(mu),mu,1)+0.15*X[,5]
  S=quantile(mu,0.5)
  len=K-S
  label=mu>S
  plot(A,Y)
  abline(h=S)
  alpha=matrix(rep(c(0.5,0.5),size),ncol=2,byrow=TRUE)
  traininfo = list(X=X,A=A,Y=Y,mu=mu,S=S,alpha=alpha,propensity=rep(1,size),K=K,opt_L=pmax(c-len,0),opt_R=pmin(c+len,K),label=label)
  return(traininfo)
}
#' @export
Scenario4.continuous <- function(size,ncov,seed,aL=-1,aU=1){
  set.seed(seed)
  X = matrix(runif(size*ncov,-1,1),ncol=ncov)
  A = runif(size,aL,aU)
  c=2*(0.15*I(X[,1]>0.5)-0.15*I(X[,2]<(-0.5))+0.15*sin(pi*X[,3])-0.15*X[,4]^3)
  #mu =-2*(A-c)^2*(abs(A-c)<0.5)-(abs(A-c)>=0.5)*0.5+0.15*log(abs(X[,5])+1)+5
  mu =(1-abs(c-A))
  Y = rnorm(length(mu),mu,0.1)#+0.15*log(abs(X[,5])+1)-0.058
  S=0.5
  len=0.5
  label=mu>S
  plot(A,Y)
  abline(h=S)
  alpha=matrix(rep(c(0.5,0.5),size),ncol=2,byrow=TRUE)
  traininfo = list(X=X,A=A,Y=Y,mu=mu,S=S,alpha=alpha,propensity=rep(1,size),opt_L=pmax(c-len,aL),opt_R=pmin(c+len,aU),label=label)
  return(traininfo)
}
#' @export
Scenario4.ordinal <- function(size,ncov,seed,K=20){
  set.seed(seed)
  levels=seq(aL-0.0001,aU+0.0001,length.out = K+1)
  X = matrix(runif(size*ncov,-1,1),ncol=ncov)
  c=(0.15*I(X[,1]>0.5)-0.15*I(X[,2]<(-0.5))+0.15*sin(pi*X[,3])-0.15*X[,4]^3)
  c=pmax(pmin(aU,c),aL)
  c=as.numeric(cut(c,breaks=levels,labels=F))
  A = sample(1:K,size,replace = T)
  mu =(K-abs(c-A))
  Y = rnorm(length(mu),mu,1)+0.15*log(abs(X[,5])+1)-0.058
  S=quantile(mu,0.5)
  len=K-S
  label=mu>S
  plot(A,Y)
  abline(h=S)
  alpha=matrix(rep(c(0.5,0.5),size),ncol=2,byrow=TRUE)
  traininfo = list(X=X,A=A,Y=Y,mu=mu,S=S,alpha=alpha,propensity=rep(1,size),K=K,opt_L=pmax(c-len,0),opt_R=pmin(c+len,K),label=label)
  return(traininfo)
}
