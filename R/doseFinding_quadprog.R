#library(CVXR)
#library(LowRankQP)
#library(quadprog)

#' @export
doseFind<-function(train, test=NULL,eps=NULL,aL=NULL,aU=NULL,step=1,thres=NULL,maxiter=20,lambda=0){
  if (is.null(thres)) thres=quantile(train$Y,0.8)
  if (is.null(aL)) aL=min(train$A)-0.5*sd(train$A)
  if (is.null(aU)) aU=max(train$A)+0.5*sd(train$A)
  p=dim(train$X)[2]
  n=dim(train$X)[1]

  index=train$Y>=thres
  beta <- Variable(p)
  alpha<-Variable(1)
  obj <- Minimize(sum_entries(abs(train$A[index] -alpha- train$X[index,] %*% beta)))
  prob<- Problem(obj)
  solv <- solve(prob)

  beta0= solv$getValue(beta)
  alpha0=solv$getValue(alpha)
  pred0=alpha0+train$X%*% beta0
  for (i in 1:maxiter){
    solv=fit.dose(train,pred0=pred0,eps=eps)
    if (sum((solv$beta-beta0)^2)<1e-5){
      cat('Converged at iteration: ',i,'-th\n')
      break;
    } else {
      beta0=solv$beta
      alpha0=solv$alpha
    }

    pred0=alpha0+train$X%*% beta0
    dosedif=pred0-train$A
    index=(abs(dosedif)<0.1)
    cat('iteration:',i,' sample size:',sum(index),'\n')
    print(quantile(pred0))
  }
  pred_dose=predict.doseCVX(test$X,solv)
  pred_dose=pmin(pmax(pred_dose,3),12)
  return(list(pred=pred_dose,solv=solv,beta=solv$beta,alpha=solv$alpha))
}

#' @export
fit.dose<-function(train,pred0=NULL,beta0=NULL,alpha0=NULL,eps,aL=NULL,aU=NULL){
  
  if (is.null(pred0)) pred0=alpha0+train$X%*% beta0
  p=dim(train$X)[2]
  Q=(abs(train$A-pred0)<eps)
  Q1=train$A<pred0-eps
  Q2=train$A>pred0+eps
  beta <- Variable(p)
  alpha<-Variable(1)
  pred=alpha+train$X %*% beta
  res=train$A - pred
  part1=sum_entries(Q*abs(res))
  part2=sum_entries((Q1)*2*pos(train$A-pred))
  part3=sum_entries((Q2)*2*pos(pred-train$A))
  partAll=sum_entries(1/train$propensity*train$Y*(part1+part2+part3))
  obj <- Minimize(partAll)
  prob<- Problem(obj)
  solv <- solve(prob)
  class(solv)='doseCVX'
  solv$beta=solv$getValue(beta)
  solv$alpha=solv$getValue(alpha)
  return(solv)
}

#' @export
predict.doseCVX<-function(X,solv){
  return(solv$alpha+X%*%solv$beta)
}


##############################################################################
## one-sided dosing interval
##############################################################################
#' @export
lineOne<-function(b,s1,s2,s3,s4,a,fit){
  res=a-fit-b
  return(sum(s1*pmax(-res,0)+s2*pmax(res,0)))
}
#' @export
fit.Int1<-function(train,eps,XX=NULL,CC=NULL,pred0=NULL,rbfkernel=NULL,
                   Intcep0=NULL,beta0=NULL,lambda=0.0001,type='PDI',lower=TRUE,method='completeLinear',...){
  # if (FALSE){
  #   alpha=c(0.82,0.18)
  #   eps=0.3
  #   lambda=0.0001
  #   n=length(train$A)
  #   pred0=6
  #
  # }
  n=length(train$A)
  p=dim(train$X)[2]

  X=apply(train$X,2,function(x) (x-mean(x))/sd(x))
  if (is.null(CC)){
    if (method=='completeLinear'){
      XX=X%*%t(X)
    } else {
      rbfkernel <- rbfdot(sigma = sigma)
      XX=as.matrix(as.data.frame(kernelMatrix(rbfkernel,X)))
    }
    CC=cbind(rbind(XX,-XX),rbind(-XX,XX))
  }

  d=matrix(c(-train$A,train$A),ncol = 1)
  Aeq=matrix(c(rep(1,n),rep(-1,n)),ncol = 1)

  if (type=='PDI'){
    H1=(train$Y>train$S)*(train$alpha[,2])*train$propensity/lambda
    H2=(train$Y<train$S)*(train$alpha[,1])*train$propensity/lambda
  } else{
    H1=pmax(train$Y-train$S,0)*(train$alpha[,2])*train$propensity/lambda
    H2=pmax(train$S-train$Y,0)*(train$alpha[,1])*train$propensity/lambda
  }

  Q1=train$A<(pred0-eps)
  Q2=train$A>(pred0+eps)
  W1=train$A>(pred0-eps)
  W2=train$A<(pred0+eps)

  #Q1=Q1*0
  #Q2=Q2*0
  if (lower){
    ub=c(H1*W1+H2*Q2,H1*Q1+H2*W2)
  } else{
    ub=c(H2*W1+H1*Q2,H2*Q1+H1*W2)
  }

  zero=(ub==0)

  ubMat=diag(2*n)

  Amat <- cbind(Aeq[!zero,],-ubMat[!zero,!zero],ubMat[!zero,!zero])
  bvec=matrix(c(0,-ub[!zero],rep(0,sum(!zero))),ncol=1)
  solv=solve.QP(Dmat=CC[!zero,!zero]+diag(0.1,sum(!zero)), dvec=d[!zero,], Amat=Amat, bvec=bvec, meq=1, factorized=FALSE)
  alfa=rep(0,2*n)
  alfa[!zero]=solv$solution
  if (method=='completeLinear'){
    beta=t((alfa[-c(1:n)]-alfa[1:n])%*%X)
    pred=X%*%beta
  } else {
    beta=(alfa[-c(1:n)]-alfa[1:n])
    pred=XX%*%beta
  }

  opt=optimize(interval = c(min(train$A),max(train$A)),f=lineOne,s1=H1*W1+H2*Q2,s2=H1*Q1+H2*W2,a=train$A,fit=pred)
  intcept=opt$minimum
  
  
  output=list(beta=beta,X=X,intcept=intcept,pred=pred+intcept,solv=solv,method=method,rbfkernel=rbfkernel)
  class(output)='Quad1'
  return(output)
}


#' @export
predict.Quad1<-function(fit,X){
  X=apply(X,2,function(x) (x-mean(x))/sd(x))
  if (fit$method=='completeLinear'){
    pred=X%*%fit$beta+fit$intcept
  } else {
    XX=as.matrix(as.data.frame(kernelMatrix(fit$rbfkernel,X,fit$X)))
    pred=XX%*%fit$beta+fit$intcept
  }
  return(pred)
}
#' @export
predict.RDCDI1<-function(obj,test){
  fit.model=obj$fit.model
  lower=obj$lower
  type=obj$type
  n1=dim(test$X)[1]
  pred0=obj$pred0

  misclass<-NULL->region->pred
  if (isTRUE(all.equal(fit.model,NA))){
    fit.model=pred0$fit
    pred=predict.init(fit.model,test)
  } else{
    pred=predict.Quad1(fit.model,test$X)
  }

  if (!is.null(test$A)){
    if (lower){
      region=pred<test$A
    } else {
      region=pred>test$A
    }
    if ((!is.null(test$S))&(!is.null(test$Y))){
      if (type=="EDI"){
        misclass=(sum(test$alpha[,1]*(region)*pmax((test$S-test$Y),0))+sum(test$alpha[,2]*(!region)*pmax((test$Y-test$S),0)))/n1
      } else if (type=="PDI"){
        misclass=(sum(test$alpha[,1]*(region)*(test$Y<test$S))+sum(test$alpha[,2]*(!region)*(test$Y>test$S)))/n1
      }
    }
  }
  return(list(misclass=misclass,region=region,pred=pred))
}
#' @export
RDCDI1<-function(train,type='PDI',method=NULL,pred0=NULL,
                 maxiter=25,lower=TRUE,eps=NULL,aL=NULL,aU=NULL,
                 lambda=0.0001,...){
  n=length(train$A)
  p=dim(train$X)[2]
  sigma=1/dim(train$X)[2]
  
  if (is.null(aL))  aL<-min(train$A)-0.5*sd(train$A)
  if (is.null(aU))  aU<-max(train$A)+0.5*sd(train$A)
  if (is.null(eps)) eps=0.1*(aU-aL)
  if (is.null(pred0)){
    if (length(table(train$A))>4){
      init0=initiate.continuous(train,train,side=1)
      pred0=init0$pred
      Intcep0=drop(coef(init0$fit))[1]
      beta0=drop(coef(init0$fit))[-1]
    } else {
      init0=initiate.binary(train,train,side=1)
      pred0=init0$pred
      Intcep0=drop(coef(init0$fit))[1]
      beta0=drop(coef(init0$fit))[-1]
    }
  } else {
    Intcep0=median(pred0)
    if (length(pred0)==1) pred0=rep(pred0,n)
    if (method=='completeLinear'){
      beta0=rep(0,p)
    } else {
      beta0=rep(0,n)
    }
  }
  
  X=apply(train$X,2,function(x) (x-mean(x))/sd(x))
  if (method=='completeLinear'){
    XX=X%*%t(X)
  } else {
    rbfkernel <- rbfdot(sigma = sigma)
    XX=as.matrix(as.data.frame(kernelMatrix(rbfkernel,X)))
  }
  CC=cbind(rbind(XX,-XX),rbind(-XX,XX))
  for (i in 1:maxiter){
    before<-Sys.time()
    fit=fit.Int1(train,XX=XX,CC=CC,eps=eps,pred0=pred0,lambda=lambda,type=type,lower=lower,method=method,rbfkernel=rbfkernel)
    after<-Sys.time()
    after-before

    fit0=fit
    Intcep0=fit$intcept
    pred0=drop(predict(fit,train$X))
    if (sum(abs(fit$beta-beta0))<1e-15){
      beta0=fit$beta
      cat('Converged at iteration: ',i-1,'-th\n')
      break;
    }
    beta0=fit$beta
    print(quantile(pred0))
  }
  #1-sum((pred0-train$opt)^2)/sum(train$opt^2)
  if (lower){
    region=pred0<train$A
  } else{
    region=pred0>train$A
  }
  if (type=="EDI"){
    if (lower){
      misclass=(sum(train$alpha[,1]*(region)*pmax((train$S-train$Y),0))+sum(train$alpha[,2]*(!region)*pmax((train$Y-train$S),0)))/n
    } else{
      misclass=(sum(train$alpha[,1]*(!regopm)*pmax((train$S-train$Y),0))+sum(train$alpha[,2]*(region)*pmax((train$Y-train$S),0)))/n
    }
  } else if (type=="PDI"){
    if (lower){
      misclass=(sum(train$alpha[,1]*(region)*(train$Y<train$S))+sum(train$alpha[,2]*(!region)*(train$Y>train$S)))/n
    } else{
      misclass=(sum(train$alpha[,1]*(!region)*(train$Y<train$S))+sum(train$alpha[,2]*(region)*(train$Y>train$S)))/n
    }
  }
  output=list(beta=fit$beta,intcept=fit$intcept,pred=pred0,solv=fit$solv,region=region,type=type,lower=lower,numiter=i,lambda=lambda,fit.model=fit,misclass=misclass)
  class(output)='RDCDI1'
  return(output)
}
#
# fit_L=RDCDI1(train=train,alpha=c(0.82,0.18),type='PDI',pred0=6,
#              maxiter=5,lower=TRUE,eps=0.2,aL=NULL,aU=NULL,
#              lambda=0.000001)
#
# pred_L=predict(fit_L,test$X)
#
# fit_R=RDCDI1(train=train,alpha=c(0.82,0.18),type='PDI',pred0=8,
#              maxiter=5,lower=FALSE,eps=0.2,aL=NULL,aU=NULL,
#              lambda=0.000001)
#
#
# pred_R=predict(fit_R,test$X)
#
# plotdata=data.frame(observed_A1c=test$A,recommended_A1c=pred_dose,
#                     lowerBound=pred_L,upperBound=pred_R)
# colnames(plotdata)[1]='observed_A1c'
# plotdata<- melt(plotdata)
# p<-ggplot(plotdata,aes(x=value, fill=variable,colour=variable)) +
#   scale_x_continuous(breaks = sort(c(seq(2,12,0.5))))+
#   geom_density(alpha=0.25)+
#   xlab('A1c')
# p

#####################################################################
## two-sided version
#####################################################################


# lineTwo<-function(b,s1,s2,s3,s4,s5,s6,s7,s8,a,fitL,fitR,eps){
#   resL1=a-fitL-b[1]
#   resL2=a-fitL-b[1]+eps
#   resR1=a-fitR-b[2]
#   resR2=a-fitR-b[2]-eps
#   return(sum(s1*pmax(-resL1,0)+s2*pmax(resL1,0)+s3*pmax(-resL2,0)+s4*pmax(resL2,0)+s5*pmax(resR1,0)+s6*pmax(-resR1,0)+s7*pmax(resR2,0)+s8*pmax(-resR2,0)))
# }
# fit.Int2<-function(train,eps,lambda=0.0001,type='PDI',CC=NULL,XX=NULL,lower=TRUE,
#                    pred0L=NULL,Intcep0L=NULL,beta0L=NULL,pred0R=NULL,Intcep0R=NULL,beta0R=NULL,...){
#   if (is.null(pred0L)) pred0L=Intcep0L+train$X %*% beta0L
#   if (is.null(pred0R)) pred0R=Intcep0R+train$X %*% beta0R
#   # library(quadprog)
#   # train=subsetData(train,sample(1:length(train$A),1000))
#   # alpha=c(0.82,0.18)
#   # eps=0.3
#   # lambda=0.000001
#   # pred0L=6.5
#   # pred0R=8
#   n=length(train$A)
#   p=dim(train$X)[2]
# 
#   if (is.null(X)) X=apply(train$X,2,function(x) (x-mean(x))/sd(x))
#   if (is.null(CC)){
#     XX=X%*%t(X)
#     col1=rbind(XX,-XX,XX,-XX,XX,matrix(0,nrow=4*n,ncol=n))
#     col2=rbind(matrix(0,nrow=4*n,ncol=n),XX,-XX,XX,-XX,XX)
#     CC=cbind(col1,-col1,col1,-col1,col1+col2,-col2,col2,-col2,col2)
#   }
# 
#   Q1L=train$A<(pred0L-eps)
#   Q2L=train$A<(pred0L)
#   W1L=train$A>(pred0L-eps)
#   W2L=train$A>(pred0L)
# 
#   Q1R=train$A>(pred0R+eps)
#   Q2R=train$A>(pred0R)
#   W1R=train$A<(pred0R+eps)
#   W2R=train$A<(pred0R)
# 
#   #Q1L=Q1L*0;Q2L=Q2L*0;Q1R=Q1R*0;Q2R=Q2R*0
#   if (type=='PDI'){
#     H1=(train$Y>train$S)*(train$alpha[,2])*train$propensity/lambda
#     H2=(train$Y<train$S)*(train$alpha[,1])*train$propensity/lambda
#   } else {
#     H1=pmax(train$Y-train$S,0)*(train$alpha[,2])*train$propensity/lambda
#     H2=pmax(train$S-train$Y,0)*(train$alpha[,1])*train$propensity/lambda
#   }
# 
# 
# 
#   ub=c(H1*W1L,H1*Q1L,H2*W2L,H2*Q2L,rep(100000,n),H1*W1R,H1*Q1R,H2*W2R,H2*Q2R)
#   zero0<-(ub==0)->zero1
#   zero1[(4*n+1):(5*n)]=TRUE
# 
# 
# 
#   d=matrix(c(-train$A,train$A,-train$A-rep(eps,n),train$A+rep(eps,n),rep(0,n),-train$A,train$A,-train$A-rep(eps,n),train$A+rep(eps,n)),ncol = 1)
#   Aeq=cbind(c(rep(1,n),rep(-1,n),rep(1,n),rep(-1,n),rep(1,n),rep(0,4*n)),c(rep(0,4*n),rep(-1,n),rep(1,n),rep(-1,n),rep(1,n),rep(-1,n)))
# 
#   ubMat=diag(9*n)
# 
#   Amat <- cbind(Aeq[!zero0,],-ubMat[!zero0,!zero1],ubMat[!zero0,!zero0])
#   bvec=matrix(c(rep(0,2),-ub[!zero1],rep(0,sum(!zero0))),ncol=1)
#   solv=solve.QP(Dmat=CC[!zero0,!zero0]+diag(0.00000001,sum(!zero0)), dvec=d[!zero0,], Amat=Amat, bvec=bvec, meq=2, factorized=FALSE)
#   alfa=rep(0,9*n)
#   alfa[!zero0]=solv$solution
#   betaL=t((alfa[(n+1):(2*n)]-alfa[(1:n)]+alfa[(3*n+1):(4*n)]-alfa[(2*n+1):(3*n)]-alfa[(4*n+1):(5*n)])%*%X)
#   betaR=t((alfa[(6*n+1):(7*n)]-alfa[(5*n+1):(6*n)]+alfa[(8*n+1):(9*n)]-alfa[(7*n+1):(8*n)]+alfa[(4*n+1):(5*n)])%*%X)
# 
#   #lineTwo(c(6,8),s1=H1*W1L+H2*Q2L,s2=H1*Q1L+H2*W2L,s3=H1*W1R+H2*Q2R,s4=H1*Q1R+H2*W2R,a=train$A,fitL=X%*%betaL,fitR=X%*%betaR)
#   fit=optim(par=c(quantile(train$A,0.2),quantile(train$A,0.8)),fn=lineTwo,eps=eps,s1=H1*W1L,s2=H1*Q1L,s3=H2*W2L,s4=H2*Q2L,s5=H1*W1R,s6=H1*Q1R,s7=H2*W2R,s8=H2*Q2R,a=train$A,fitL=X%*%betaL,fitR=X%*%betaR)
# 
# 
#   # mu=ub[-c((2*n+1):(3*n))]-alfa[-c((2*n+1):(3*n))]
#   # aa=c(train$A-X%*%betaL,X%*%betaL-train$A,train$A-X%*%betaR,X%*%betaR-train$A)
#   # aa[which((mu!=0)&(abs(alfa)[-c((2*n+1):(3*n))]>0.0000001))]
# 
#   output=list()
#   output$betaL=betaL
#   output$betaR=betaR
#   output$intcepL=fit$par[1]
#   output$intcepR=fit$par[2]
#   output$pred_L=X%*%betaL+fit$par[1]
#   output$pred_R=X%*%betaR+fit$par[2]
#   output$fit=fit
#   #output$solv=solv
#   class(output)='Quad2'
#   return(output)
# }
#' @export
lineTwo<-function(b,s1,s2,s3,s4,a,fitL,fitR){
  resL=a-fitL-b[1]
  resR=a-fitR-b[2]
  return(sum(s1*pmax(-resL,0)+s2*pmax(resL,0)+s3*pmax(-resR,0)+s4*pmax(resR,0)+10*pmax(fitL+b[1]-fitR-b[2],0)))
}
#' @export
fit.Int2<-function(train,eps,lambda=0.0001,type='PDI',CC=NULL,XX=NULL,method='completeLinear',rbfkernel=NULL,
                   pred0L=NULL,Intcep0L=NULL,beta0L=NULL,pred0R=NULL,Intcep0R=NULL,beta0R=NULL,...){
  if (is.null(pred0L)) pred0L=Intcep0L+train$X %*% beta0L
  if (is.null(pred0R)) pred0R=Intcep0R+train$X %*% beta0R
  #library(quadprog)
  #train=subsetData(train,sample(1:length(train$A),1000))
  #alpha=c(0.82,0.18)
  #eps=0.3
  #lambda=0.000001
  n=length(train$A)
  sigma=1/dim(train$X)[2]
  #pred0L=6.5
  #pred0R=8
  X=apply(train$X,2,function(x) (x-mean(x))/sd(x))
  if (is.null(CC)){
    if (method=='completeLinear'){
      XX=X%*%t(X)
    } else {
      rbfkernel <- rbfdot(sigma = sigma)
      XX=as.matrix(as.data.frame(kernelMatrix(rbfkernel,X)))
    }
    col1=rbind(XX,-XX,XX,matrix(0,nrow=2*n,ncol=n))
    col2=rbind(matrix(0,nrow=2*n,ncol=n),XX,-XX,XX)
    CC=cbind(col1,-col1,col1+col2,-col2,col2)
  }
  

  Q1L=train$A<(pred0L-eps)
  Q2L=train$A<(pred0L+eps)
  W1L=train$A>(pred0L-eps)
  W2L=train$A>(pred0L+eps)

  Q1R=train$A>(pred0R+eps)
  Q2R=train$A>(pred0R-eps)
  W1R=train$A<(pred0R+eps)
  W2R=train$A<(pred0R-eps)

  #Q1L=Q1L*0;Q2L=Q2L*0;Q1R=Q1R*0;Q2R=Q2R*0
  if (type=='PDI'){
    H1=(train$Y>train$S)*(train$alpha[,2])*train$propensity/lambda
    H2=(train$Y<train$S)*(train$alpha[,1])*train$propensity/lambda
  } else{
    H1=pmax(train$Y-train$S,0)*(train$alpha[,2])*train$propensity/lambda
    H2=pmax(train$S-train$Y,0)*(train$alpha[,1])*train$propensity/lambda
  }

  ub=c(H1*W1L+H2*W2L,H1*Q1L+H2*Q2L,rep(1e05,n),H1*Q1R+H2*Q2R,H1*W1R+H2*W2R)
  zero0<-(ub==0)->zero1
  zero1[(2*n+1):(3*n)]=TRUE



  d=matrix(c(-train$A,train$A,rep(2*eps,n),-train$A,train$A),ncol = 1)
  Aeq=cbind(c(rep(1,n),rep(-1,n),rep(1,n),rep(0,2*n)),c(rep(0,2*n),rep(-1,n),rep(1,n),rep(-1,n)))

  ubMat=diag(5*n)

  Amat <- cbind(Aeq[!zero0,],-ubMat[!zero0,!zero1],ubMat[!zero0,!zero0])
  bvec=matrix(c(rep(0,2),-ub[!zero1],rep(0,sum(!zero0))),ncol=1)
  solv=solve.QP(Dmat=CC[!zero0,!zero0]+diag(0.1,sum(!zero0)), dvec=d[!zero0,], Amat=Amat, bvec=bvec, meq=2, factorized=FALSE)
  alfa=rep(0,5*n)
  alfa[!zero0]=solv$solution


  if (method=='completeLinear'){
    betaL=t((alfa[(n+1):(2*n)]-alfa[(1:n)]-alfa[(2*n+1):(3*n)])%*%X)
    betaR=t((alfa[(4*n+1):(5*n)]-alfa[(3*n+1):(4*n)]+alfa[(2*n+1):(3*n)])%*%X)
    pred_L=X%*%betaL
    pred_R=X%*%betaR
  } else {
    betaL=(alfa[(n+1):(2*n)]-alfa[(1:n)]-alfa[(2*n+1):(3*n)])
    betaR=(alfa[(4*n+1):(5*n)]-alfa[(3*n+1):(4*n)]+alfa[(2*n+1):(3*n)])
    pred_L=XX%*%betaL
    pred_R=XX%*%betaR
  }


  #lineTwo(c(6,8),s1=H1*W1L+H2*Q2L,s2=H1*Q1L+H2*W2L,s3=H1*W1R+H2*Q2R,s4=H1*Q1R+H2*W2R,a=train$A,fitL=X%*%betaL,fitR=X%*%betaR)
  fit=optim(par=c(quantile(train$A,0.2),quantile(train$A,0.8)),fn=lineTwo,s1=H1*W1L+H2*W2L,s2=H1*Q1L+H2*Q2L,s3=H1*Q1R+H2*Q2R,s4=H1*W1R+H2*W2R,a=train$A,fitL=pred_L,fitR=pred_R)

  # mu1L=ub[c((1):(n))]-alfa[c((1):(n))]
  # mu2L=ub[c((n+1):(2*n))]-alfa[c((n+1):(2*n))]
  # index1L=!zero0[c((1):(n))]
  # bLub=(-X%*%betaL-train$A)[index1L&(abs(mu1L)<0.00000001)]
  # index2L=!zero0[c((1+n):(2*n))]
  # bLlb=(train$A-X%*%betaL)[index2L&(abs(mu2L)<0.00000001)]
  #
  # ##
  # ##  # this is supposed to be the right one
  # bLlb=(train$A-X%*%betaL)[((alfa[c((1):(n))]>1e5)&index1L)|((mu2L>1e5)&index2L)]
  # quantile(bLlb)
  # bLub=(X%*%betaL-train$A)[((alfa[c((n+1):(2*n))]>1e5)&index2L)|((mu1L>1e5)&index1L)]
  # quantile(bLub)
  #
  # mu1L=ub[c((1):(n))]-alfa[c((1):(n))]
  # mu2L=ub[c((n+1):(2*n))]-alfa[c((n+1):(2*n))]
  # index1L=!zero0[c((1):(n))]
  # bLub=(-X%*%betaL-train$A)[index1L&abs(mu2L)==0]
  # bLub=(-X%*%betaL-train$A)[index1L&(abs(alfa[c((1):(n))])<1e-10)]
  # bLub=(-X%*%betaL-train$A)[index1L&(alfa[c((1):(n))]>1e-2)&(mu1L>1e-2)]
  #
  # index2L=!zero0[c((1+n):(2*n))]
  # bLlb=(train$A-X%*%betaL)[index2L&abs(mu2L)==0]
  # index2L=!zero0[c((1+n):(2*n))]
  # bLlb=(train$A-X%*%betaL)[index2L&(alfa[c((n+1):(2*n))]>1e-2)&(mu2L>1e-2)]
  # bLlb=(train$A-X%*%betaL)[index2L&(alfa[c((n+1):(2*n))]>1e-2)&(mu2L>1e-2)]
  #
  # index=(abs(mu)<0.00001)&(ub[-c((2*n+1):(3*n))]!=0)
  # aa=c(train$A-X%*%betaL,X%*%betaL-train$A,train$A-X%*%betaR,X%*%betaR-train$A)
  # aa[index]
  # aa[which((mu!=0)&(abs(alfa)[-c((2*n+1):(3*n))]>0.0000001))]

  output=list()
  output$betaL=betaL
  output$betaR=betaR
  output$solv=solv
  output$intcepL=fit$par[1]
  output$intcepR=fit$par[2]
  output$pred_L=pred_L+fit$par[1]
  output$pred_R=pred_R+fit$par[2]
  output$rbfkernel=rbfkernel
  output$method=method
  output$X=X
  #output$solv=solv
  class(output)='Quad2'
  return(output)
}
#' @export
RDCDI2<-function(train,method='completeLinear',type='PDI',pred0L=NULL,pred0R=NULL,
                 maxiter=25,eps=NULL,aL=NULL,aU=NULL,lower=TRUE,
                 lambda=0.0001,...){

  n=length(train$A)
  p=dim(train$X)[2]
  
  rbfkernel=NULL

  X=apply(train$X,2,function(x) (x-mean(x))/sd(x))
  if (method=='completeLinear'){
    XX=X%*%t(X)
  } else {
    rbfkernel <- rbfdot(sigma = 1/dim(train$X)[2])
    XX=as.matrix(as.data.frame(kernelMatrix(rbfkernel,X)))
  }

  col1=rbind(XX,-XX,XX,matrix(0,nrow=2*n,ncol=n))
  col2=rbind(matrix(0,nrow=2*n,ncol=n),XX,-XX,XX)
  CC=cbind(col1,-col1,col1+col2,-col2,col2)

  if (is.null(aL))  aL<-min(train$A)-0.5*sd(train$A)
  if (is.null(aU))  aU<-max(train$A)+0.5*sd(train$A)
  if (is.null(eps)) eps=0.2*sd(train$A)

  if (is.null(pred0L)|is.null(pred0R)){
    if (length(table(train$A))>4){
      init0=initiate.continuous(train,train,side=2)
    } else {
      init0=initiate.binary(train,train,side=2)
    }
    pred0L=fit$pred_L
    pred0R=fit$pred_R
    Intcep0L=drop(coef(init0$fit$fit_L))[1]
    beta0L=drop(coef(init0$fit$fit_L))[-1]
    Intcep0R=drop(coef(init0$fit$fit_R))[1]
    beta0R=drop(coef(init0$fit$fit_R))[-1]
  } else {
    if (method=='completeLinear'){
      beta0L=rep(0,p)
      beta0R=rep(0,p)
    } else {
      beta0L=rep(0,n)
      beta0R=rep(0,n)
    }
    Intcep0L=median(pred0L)
    Intcep0R=median(pred0R)
    if (length(pred0L)==1) pred0L=rep(pred0L,n)
    if (length(pred0R)==1) pred0R=rep(pred0R,n)
  }
  #pred0L=rep(-0.5,n);pred0R=rep(0.5,n)
  for (i in 1:maxiter){
    before=Sys.time()
    fit=fit.Int2(train,eps=eps,lambda=lambda,type=type,pred0L=pred0L,pred0R=pred0R,rbfkernel=rbfkernel,CC=CC,XX=XX,method=method)
    after=Sys.time()
    print(after-before)

    pred0L=fit$pred_L
    pred0R=fit$pred_R
    if (sum(abs(fit$betaL-beta0L)+abs(fit$betaR-beta0R))<1e-15){
      cat('Converged at iteration: ',i-1,'-th\n')
      beta0L=fit$betaL
      beta0R=fit$betaR
      break;
    }
    beta0L=fit$betaL
    beta0R=fit$betaR

    cat('two-sided iteration: ',i,'\n')
    print(quantile(pred0L))
    print(quantile(pred0R))
  }


  region=(pred0L<train$A)&(pred0R>train$A)

  if (type=="EDI"){
    misclass=(sum(train$alpha[,1]*(region)*pmax((train$S-train$Y),0))+sum(train$alpha[,2]*(!region)*pmax((train$Y-train$S),0)))/n
  } else if (type=="PDI"){
    misclass=(sum(train$alpha[,1]*(region)*(train$Y<train$S))+sum(train$alpha[,2]*(!region)*(train$Y>train$S)))/n
  }
  output=list(method=method,betaL=fit$betaL,betaR=fit$betaR,intcepL=fit$intcepL,intcepR=fit$intcepR,pred_L=fit$pred_L,pred_R=fit$pred_R,region=region,type=type,lower=lower,numiter=i,lambda=lambda,fit.model=fit,misclass=misclass)
  class(output)='RDCDI2'
  return(output)
}
#' @export
predict.Quad2<-function(fit,X){
  X=apply(X,2,function(x) (x-mean(x))/sd(x))
  if (fit$method=='completeLinear'){
    pred_L=fit$intcepL+X %*% fit$betaL
    pred_R=fit$intcepR+X %*% fit$betaR
  } else {
    XX=as.matrix(as.data.frame(kernelMatrix(fit$rbfkernel,X,fit$X)))
    pred_L=fit$intcepL+XX %*% fit$betaL
    pred_R=fit$intcepR+XX %*% fit$betaR
  }
  return(list(pred_L=pred_L,pred_R=pred_R))
}
#' @export
predict.RDCDI2<-function(obj,test){
  fit.model=obj$fit.model
  lower=obj$lower
  type=obj$type
  n1=dim(test$X)[1]
  pred0=obj$pred0

  misclass<-NULL->region->pred
  if (isTRUE(all.equal(fit.model,NA))){
    fit.model=pred0$fit
    pred=predict.init(fit.model,test)
  } else{
    pred=predict.Quad2(fit.model,test$X)
    pred_L=pred$pred_L
    pred_R=pred$pred_R
  }

  if (!is.null(test$A)){
    region=(test$A<pred_R)&(test$A>pred_L)
    if ((!is.null(test$S))&(!is.null(test$Y))){
      if (type=="EDI"){
        misclass=(sum(test$alpha[,1]*(region)*pmax((test$S-test$Y),0))+sum(test$alpha[,2]*(!region)*pmax((test$Y-test$S),0)))/n1
      } else if (type=="PDI"){
        misclass=(sum(test$alpha[,1]*(region)*(test$Y<test$S))+sum(test$alpha[,2]*(!region)*(test$Y>test$S)))/n1
      }
    }
  }
  return(list(misclass=misclass,region=region,pred_L=pred_L,pred_R=pred_R))
}

# fit_two=RDCDI2(train=train,alpha=c(0.85,0.15),type='PDI',pred0=c(5.5,7.6),
#                maxiter=5,eps=0.2,aL=NULL,aU=NULL,
#                lambda=0.000000001)
#
# pred_L=predict(fit_two,test$X)$pred_L
# pred_R=predict(fit_two,test$X)$pred_R
#
# plotdata=data.frame(observed_A1c=test$A,recommended_A1c=pred_dose,
#                     lowerBound=pred_L,upperBound=pred_R)
# colnames(plotdata)[1]='observed_A1c'
# plotdata<- melt(plotdata)
# p<-ggplot(plotdata,aes(x=value, fill=variable,colour=variable)) +
#   scale_x_continuous(breaks = sort(c(seq(2,12,0.5))))+
#   geom_density(alpha=0.25)+
#   xlab('A1c')
# p



###############################################################





















