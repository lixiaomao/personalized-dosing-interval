#' Fits dose interval using DC-algorithm
#'
#' This function is the wrapper of fitting functions for doseInt
#' 
#' @param train the training data set.
#' @param pred0 the initial prediction for training data set, if NULL, default initialization is used
#' @param nfolds the number of folds used in cross-validation.
#' @param alpha the guaranteed probability in PDI and the relative loss in EDI. Default is c(0.5,0.5).
#' @param type indicates whether PDI or EDI is to be calculated.
#' @param two.sided indicator of whether two-sided interval is considered.
#' @param family specifis the methods. 'continuous' leads to DC-algorithm,'ordinal' handles the ordinal treatments, while 'as.ordinal' cuts continuous treatment into several ordinal levels and apply the ordinal algorithm.
#' @param method specified methods used for DC and ordinal approaches.
#' @param maxiter maximal iterations for DC algorithm
#' @param step stepsize for DC algorithm
#' @param trace the level of detail to be printed: 1,2,3,4
#' @param lower True or False, Whether lower boundary is considered when applied to one-sided interval
#' @param breaks ways 'as.ordinal' cuts continuous treatment into several ordinal levels, default 'quantile'. User can also choose 'uniform' or specify customized breaks
#' @param Eps parameter for the relaxation
#' @param Lambda  penalty for the quantile regression
#' @param Cost  cost for the support vector machine
#' @param Embed.mtry  proportion of mtry
#'    
#' @examples
#' 
#' 
#' 


DCDI1<-function(train,alpha,method='rq',type='PDI',pred0=NULL,
                maxiter=20,step=1,trace=0,lower=TRUE, eps=NULL,aL=NULL,aU=NULL,
                lambda=0.0000001,cost=1,embed.mtry = 1/2,svmKernel=NULL,...){
  n=length(train$A)
  fit.method=switch(method,'rq'=fit.rq.model,'svmLinear'=fit.svmL.model,'svmRadial'=fit.svmK.model,'RLT'=fit.rlt.model)
  
  para=list(...)
  if (is.null(aL))  aL<-min(train$A)-0.5*sd(train$A)
  if (is.null(aU))  aU<-max(train$A)+0.5*sd(train$A)
  if (is.null(eps)) eps<-para$eps else eps=0.1*(aU-aL)
  if (is.null(pred0)){
    if (length(table(train$A))>4){
      pred0=initiate.continuous(train,train,side=1)
    } else {
      pred0=initiate.binary(train,train,side=1)
    }
  }
  pred=pred0$pred
  predOld=pred
  
  W1=((pred-train$A)<eps)&((pred-train$A)>0)
  W2=((train$A-pred)<eps)&((train$A-pred)>0)
  id=which(W1|W2)
  
  indexOld1=which(W1)
  indexOld2=which(W2)
  indexNew1<-NULL->indexNew2
  
  region=(train$A>predOld)
  if (type=="EDI"){
    if (lower){
      misclass0=(alpha[1]*sum((region)*pmax((train$S-train$Y),0))+alpha[2]*sum((!region)*pmax((train$Y-train$S),0)))/n
    } else{
      misclass0=(alpha[1]*sum((!region)*pmax((train$S-train$Y),0))+alpha[2]*sum((region)*pmax((train$Y-train$S),0)))/n
    }
  } else if (type=="PDI"){
    if (lower){
      misclass0=(alpha[1]*sum((region)*(train$Y<train$S))+alpha[2]*sum((!region)*(train$Y>train$S)))/n
    } else{
      misclass0=(alpha[1]*sum((!region)*(train$Y<train$S))+alpha[2]*sum((region)*(train$Y>train$S)))/n
    }
  }
  if (trace) cat('Initial Error:',misclass0,'\n')
  misclass_old=misclass0
  misclass=misclass_old
  
  fit.model_old=pred0$fit
  fit.model=NA
  i=1
  for (i in 1:maxiter){
    if (trace) cat("\nIteration:",i,'\n')
    if (type=="EDI"){
      if (lower){
        weights=train$propensity*(alpha[2]*W1*pmax((train$Y-train$S),0)+alpha[1]*W2*pmax((train$S-train$Y),0))
      } else {
        weights=train$propensity*(alpha[2]*W2*pmax((train$Y-train$S),0)+alpha[1]*W1*pmax((train$S-train$Y),0))
      }
    } else if(type=="PDI"){
      if (lower){
        weights=train$propensity*(alpha[2]*W1*(train$Y>train$S)+alpha[1]*W2*(train$S>train$Y))
      } else{
        weights=train$propensity*(alpha[2]*W2*(train$Y>train$S)+alpha[1]*W1*(train$S>train$Y))
      }
    }
    if (trace) cat('Percentage used: ',sum(weights>0),'(',round(mean(weights>0)*100),'%) \n')
    if (trace) print(quantile(train$A[weights!=0]))
    fit.model=fit.method(train,weights=weights,eps=eps,lambda=lambda,cost=cost,embed.mtry = embed.mtry, pred=predOld,side="left")
    
    if (isTRUE(all.equal(fit.model,NA))){
      cat("Fitting failed at interation: ",i,'\n')
      break;
    }
    pred=predict.owl(fit.model,train)
    pred=step*as.vector(pred)+(1-step)*as.vector(predOld)
    predOld=pred
    
    W1=((pred-train$A)<eps)&((pred-train$A)>0)
    W2=((train$A-pred)<eps)&((train$A-pred)>0)
    
    indexNew1=which(W1)
    indexNew2=which(W2)
    
    region=pred<train$A
    if (type=="EDI"){
      if (lower){
        misclass=(alpha[1]*sum((region)*pmax((train$S-train$Y),0))+alpha[2]*sum((!region)*pmax((train$Y-train$S),0)))/n
      } else{
        misclass=(alpha[1]*sum((!regopm)*pmax((train$S-train$Y),0))+alpha[2]*sum((region)*pmax((train$Y-train$S),0)))/n
      }
    } else if (type=="PDI"){
      if (lower){
        misclass=(alpha[1]*sum((region)*(train$Y<train$S))+alpha[2]*sum((!region)*(train$Y>train$S)))/n
      } else{
        misclass=(alpha[1]*sum((!region)*(train$Y<train$S))+alpha[2]*sum((region)*(train$Y>train$S)))/n
      }
    }
    if (trace>2) print(paste("Training Error:",misclass))
    if ((misclass>misclass_old*1.2)){
      cat("Fitting led to larger Error, stopped at interation: ",i,'\n')
      break
    }
    
    if (trace>2){
      plot(train$Y~train$A,main=paste("iteration:",i),cex=1.5,pch=21)
      points(train$Y[region]~train$A[region],col="red",bg='red2',pch=21,cex=1.5)
      points(train$Y[!region]~train$A[!region],col="blue",bg='royalblue3',pch=21,cex=1.5)
      points(train$Y[W1]~train$A[W1],col="orange",cex=1.5,pch=21)
      points(train$Y[W2]~train$A[W2],col="green",cex=1.5,pch=21)
      Sys.sleep(0.5)
    }
    
    if (identical(indexOld1,indexNew1)&identical(indexOld2,indexOld2)){
      if (trace) print(paste("Converge after",i,"steps"))
      break
    }
    indexOld1=indexNew1
    indexOld2=indexNew2
  } # end iteration
  if (trace) print(paste("Iteration stopped:",i))
  obj=list(region=region,method=method,type=type,alpha=alpha,lower=lower,numiter=i,lambda=lambda,fit.model=fit.model,W1=W1,W2=W2,misclass=misclass,pred=pred,weights=weights)
  class(obj)='DCDI1'
  return(obj)
}

predict.DCDI1<-function(obj,test){
  fit.model=obj$fit.model
  alpha=obj$alpha
  lower=obj$lower
  type=obj$type
  n1=dim(test$X)[1]
  
  misclass<-NULL->region->pred
  if (isTRUE(all.equal(fit.model,NA))){
    fit.model=pred0$fit
    pred=predict.init(fit.model,test)
  } else{
    pred=predict.owl(fit.model,test)
  }
  if (!is.null(test$A)){
    region=pred<test$A
    if ((!is.null(test$S))&(!is.null(test$Y))){
      if (type=="EDI"){
        if (lower){
          misclass=(alpha[1]*sum((region)*pmax((test$S-test$Y),0))+alpha[2]*sum((!region)*pmax((test$Y-test$S),0)))/n1
        } else{
          misclass=(alpha[1]*sum((!region)*pmax((test$S-test$Y),0))+alpha[2]*sum((region)*pmax((test$Y-test$S),0)))/n1
        }
      } else if (type=="PDI"){
        if (lower){
          misclass=(alpha[1]*sum((region)*(test$Y<test$S))+alpha[2]*sum((!region)*(test$Y>test$S)))/n1
        } else {
          misclass=(alpha[1]*sum((!region)*(test$Y<test$S))+alpha[2]*sum((region)*(test$Y>test$S)))/n1
        }
      }
    }
  } 
  return(list(misclass=misclass,region=region,pred=pred))
}

DCDI2<-function(train,alpha,method='rq',type='PDI',pred0=NULL,
                maxiter=20,step=1,trace=0, eps=NULL,aL=NULL,aU=NULL,
                lambda=0.0000001,cost=1,embed.mtry = 1/2,svmKernel=NULL,...){
  n=length(train$A)
  fit.method=switch(method,'rq'=fit.rq.model,'svmLinear'=fit.svmL.model,'svmRadial'=fit.svmK.model,'RLT'=fit.rlt.model)
  
  para=list(...)
  if (is.null(aL)) aL<-para$aL else aL<-min(train$A)-0.5*sd(train$A)
  if (is.null(aU)) aU<-para$aU else aU<-max(train$A)+0.5*sd(train$A)
  if (is.null(eps)) eps<-para$eps else eps=0.1*(aU-aL)
  if (is.null(pred0)){
    if (length(table(train$A)>4)){
      pred0=initiate.continuous(train,train,side=2)
    } else {
      pred0=initiate.binary(train,train,side=2)
    }
  }
  
  pred_L=pred0$pred_L
  pred_R=pred0$pred_R
  predOld_L=pred_L
  predOld_R=pred_R
  
  W1_L=((pred_L-train$A)<eps)&((pred_L-train$A)>(0))
  W2_L=((train$A-pred_L)<eps)&((train$A-pred_L)>(0))
  id_L=which(W1_L|W2_L)
  W1_R=((pred_R-train$A)<eps)&((pred_R-train$A)>(0))
  W2_R=((train$A-pred_R)<eps)&((train$A-pred_R)>(0))
  id_R=which(W1_R|W2_R)
  
  region=(train$A>pred_L)&(train$A<pred_R)
  if (type=="EDI"){
    misclass0=(alpha[1]*sum(region*pmax((train$S-train$Y),0))+alpha[2]*sum((!region)*pmax((train$Y-train$S),0)))/n
  } else if (type=="PDI"){
    misclass0=(alpha[1]*sum(region*(train$Y<train$S))+alpha[2]*sum((!region)*(train$Y>train$S)))/n
  }
  if (trace) cat('Initial Error:',misclass0,'\n')
  misclass_old=misclass0
  misclass=misclass_old
  
  indexOld1_L=which(W1_L); 
  indexOld2_L=which(W2_L)
  indexNew1_L<-NULL->indexNew2_L
  indexOld1_R=which(W1_R); 
  indexOld2_R=which(W2_R)
  indexNew1_R<-NULL->indexNew2_R
  
  fit.model_L_old=pred0$fit$fit_L
  fit.model_R_old=pred0$fit$fit_R
  fit.model_L<-NA->fit.model_R
  i=1
  flag_L<-TRUE->flag_R
  for (i in 1:maxiter){
    if (trace) cat("\nIteration:",i,'\n')
    if (type=="EDI"){
      weights_L=train$propensity*(alpha[2]*W1_L*pmax((train$Y-train$S),0)+alpha[1]*W2_L*pmax((train$S-train$Y),0))
    } else if(type=="PDI"){
      weights_L=train$propensity*(alpha[2]*W1_L*(train$Y>train$S)+alpha[1]*W2_L*(train$S>train$Y))
    }
    if (type=="EDI"){
      weights_R=train$propensity*(alpha[2]*W2_R*pmax((train$Y-train$S),0)+alpha[1]*W1_R*pmax((train$S-train$Y),0))
    } else if(type=="PDI"){
      weights_R=train$propensity*(alpha[2]*W2_R*(train$Y>train$S)+alpha[1]*W1_R*(train$S>train$Y))
    }
    
    if (flag_L){
      index_L=(train$A<pred_R)
      weights_L=weights_L[index_L]
      train_L=list(X=train$X[index_L,],A=train$A[index_L],Y=train$Y[index_L],S=train$S)
      if (trace) cat('Left boundary observations used: ',sum(weights_L>0),'(',round(mean(weights_L>0)*100),'%) \n')
      fit.model_L=fit.method(train_L,weights=weights_L,eps=eps,lambda=lambda,cost=cost,embed.mtry = embed.mtry)
      if (isTRUE(all.equal(fit.model_L,NA))){
        cat("Left fitting failed at interation: ",i,'\n')
        flag_L=FALSE
      } else {
        pred_L=predict.owl(fit.model_L,train)
        pred_L=step*as.vector(pred_L)+(1-step)*as.vector(predOld_L)
      }
    }
    
    if (flag_R){
      index_R=(train$A>pred_L)
      weights_R=weights_R[index_R]
      if (trace) cat('Right boundary observations used: ',sum(weights_R>0),'(',round(mean(weights_R>0)*100),'%) \n')
      train_R=list(X=train$X[index_R,],A=train$A[index_R],Y=train$Y[index_R],S=train$S)
      fit.model_R=fit.method(train_R,weights=weights_R,eps=eps,lambda=lambda,cost=cost,embed.mtry = embed.mtry )
      if (isTRUE(all.equal(fit.model_R,NA))){
        cat("Right fitting failed at interation: ",i,'\n')
        flag_R=FALSE;
      } else{
        pred_R=predict.owl(fit.model_R,train)
        pred_R=step*as.vector(pred_R)+(1-step)*as.vector(predOld_R)
      }
    }
    
    if (!(flag_L|flag_R)){
      print(paste("Both fitting failed at iteration:",i))
      break;
    }
    W1_L=((pred_L-train$A)<eps)&((pred_L-train$A)>0)
    W2_L=((train$A-pred_L)<eps)&((train$A-pred_L)>0)
    indexNew1_L=which(W1_L)
    indexNew2_L=which(W2_L)
    W1_R=((pred_R-train$A)<eps)&((pred_R-train$A)>0)
    W2_R=((train$A-pred_R)<eps)&((train$A-pred_R)>0)
    indexNew1_R=which(W1_R)
    indexNew2_R=which(W2_R)
    
    region=(train$A>pred_L)&(train$A<pred_R)
    if (type=="EDI"){
      misclass=(alpha[1]*sum(region*pmax((train$S-train$Y),0))+alpha[2]*sum((!region)*pmax((train$Y-train$S),0)))/n
    } else if (type=="PDI"){
      misclass=(alpha[1]*sum(region*(train$Y<train$S))+alpha[2]*sum((!region)*(train$Y>train$S)))/n
    }
    
    if (misclass>misclass_old*1.2){
      cat("Fitting led to larger Error, stopped at interation:",i,'\n')
      pred_L=predOld_L
      pred_R=predOld_R
      fit.model_L=fit.model_L_old
      fit.model_R=fit.model_R_old
      break
    } else{
      misclass_old=misclass
      predOld_L=pred_L
      predOld_R=pred_R
      fit.model_L_old=fit.model_L
      fit.model_R_old=fit.model_R
    }
    if (trace>2) print(paste("Training Error:",misclass))
    if (trace>2){
      plot(train$Y~train$A,main=paste("iteration:",i),cex=1.5,pch=21)
      points(train$Y[region]~train$A[region],cex=1.5,col="red",pch=21,bg="red2")
      points(train$Y[!region]~train$A[!region],cex=1.5,col="blue",pch=21,bg="royalblue3")
      points(train$Y[W1_L]~train$A[W1_L],col="orange",cex=1.5,pch=21)
      points(train$Y[W2_L]~train$A[W2_L],col="green",cex=1.5,pch=21)
      points(train$Y[W1_R]~train$A[W1_R],col="orange",cex=1.5,pch=21)
      points(train$Y[W2_R]~train$A[W2_R],col="green",cex=1.5,pch=21)
      Sys.sleep(0.5)
    }
    if (identical(indexOld1_L,indexNew1_L)&identical(indexOld2_L,indexOld2_L)&identical(indexOld1_R,indexNew1_R)&identical(indexOld2_R,indexOld2_R)){
      break # converged
    }
    indexOld1_L=indexNew1_L
    indexOld2_L=indexNew2_L
    indexOld1_R=indexNew1_R
    indexOld2_R=indexNew2_R
  } # end of iterations
  
  if (trace) print(paste("Iteration stopped:",i))
  if (isTRUE(all.equal(fit.model_L,NA))){
    if (trace) cat("No improvement: left boundary \n")
  }
  if (isTRUE(all.equal(fit.model_R,NA))){
    if (trace) cat("No improvement: right boundary \n")
  }
  obj=list(region=region,method=method,type=type,alpha=alpha,numiter=i,lambda=lambda,fit.model_L=fit.model_L,fit.model_R=fit.model_R,W1_L=W1_L,W2_L=W2_L,W1_R=W1_R,W2_R=W2_R,misclass=misclass,pred_L=pred_L,pred_R=pred_R,weights_L=weights_L,weights_R=weights_R)
  class(obj)='DCDI2'
  return(obj)
}

predict.DCDI2<-function(obj,test){
  fit.model_L=obj$fit.model_L
  fit.model_R=obj$fit.model_R
  alpha=obj$alpha
  lower=obj$lower
  n1=dim(test$X)[1]
  type=obj$type
  pred
  region<-NULL->misclass->pred_R->pred_L
  if (isTRUE(all.equal(fit.model_L,NA))){
    fit.model_L=pred0$fit$fit_L
    pred_L=predict.init(fit.model_L,test)
  } else{
    pred_L=predict.owl(fit.model_L,test)
  }
  if (isTRUE(all.equal(fit.model_R,NA))){
    fit.model_R=pred0$fit$fit_R
    pred_R=predict.init(fit.model_R,test)
  } else{
    pred_R=predict.owl(fit.model_R,test)
  }
  if (!is.null(test$A)){
    region=(test$A>pred_L)&(test$A<pred_test_R)
    if ((!is.null(test$Y))&(!is.null(test$S))){
      if (type=="EDI"){
        misclass=(alpha[1]*sum((region)*pmax((test$S-test$Y),0))+alpha[2]*sum((!region)*pmax((test$Y-test$S),0)))/n1
      } else if (type=="PDI"){
        misclass=(alpha[1]*sum((region)*(test$Y<test$S))+alpha[2]*sum((!region)*(test$Y>test$S)))/n1
      }
    }
  }
  return(list(region=region,misclass=misclass,pred_L=pred_L,pred_R=pred_R))
}