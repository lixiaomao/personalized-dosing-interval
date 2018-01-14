#' Fits dose interval using ordinal approach
#'
#' This function is the wrapper of fitting functions for doseInt
#' 
#' @param train the training data set.
#' @param test the testing data set. If test is NULL, use train as test.
#' @param pred0 the initial prediction for training data set, if NULL, default initialization is used
#' @param nfolds the number of folds used in cross-validation.
#' @param alpha the guaranteed probability in PDI and the relative loss in EDI. Default is c(0.5,0.5).
#' @param type indicates whether PDI or EDI is to be calculated.
#' @param two.sided indicator of whether two-sided interval is considered.
#' @param family specifis the methods. 'continuous' leads to DC-algorithm,'ordinal' handles the ordinal treatments, while 'as.ordinal' cuts continuous treatment into several ordinal levels and apply the ordinal algorithm.
#' @param method specified methods used for DC and ordinal approaches.
#' @param trace the level of detail to be printed: 1,2,3,4
#' @param lower True or False, Whether lower boundary is considered when applied to one-sided interval
#' @param K level of ordinal treatment, default 10
#' @param global if global data are used for ordinal approach, default False
#' @param margin how large is the margin if global=False, default 0
#' @param breaks ways 'as.ordinal' cuts continuous treatment into several ordinal levels, default 'quantile'. User can also choose 'uniform' or specify customized breaks
#' @param Cost  cost for the support vector machine
#'    
#' @examples
#' 
#' 
#' 

find.ordinal<-function(x){
  if (all(x!=0)){
    return(NA)
  } else {
    return(which(x==0))
  }
}
predict.ordinal<-function(fit,X,K){
  Pred=NULL
  for (k in 1:(K-1)){
    tmp=matrix(0,ncol=K-1,nrow=nrow(X))
    tmp[,k]=1
    pred=predict(fit,cbind(X,tmp))
    Pred=cbind(Pred,as.numeric(pred)-1)
  }
  total=apply(Pred,1,sum)+1
  return(list(labels=total,Prediction=Pred))
}
fit.ordinal<-function(train,TEST,two.sided=FALSE,K,continuous=TRUE,global,margin=0,alpha,breaks='quantile',kernel='linear'){
  values=NULL
  n=length(train$A)
  if (continuous){
    if (length(breaks)>1){
      levels=breaks
      for (i in 1:K) values=c(values,mean(levels[i:(i+1)]))
    } else if (breaks=='uniform'){
      levels=seq(aL,aU,length.out = K+1)
      for (i in 1:K) values=c(values,mean(levels[i:(i+1)]))
    } else if (breaks=='quantile'){
      levels=quantile(train$A,probs=seq(0,1,length.out = K+1))
      for (i in 1:K) values=c(values,mean(with(train,A[(A<levels[i+1])&(A>=levels[i])])))
    }
    levels[1]=min(train$A)-0.00001
    levels[K+1]=max(train$A)+0.00001
    train_levels=cut(train$A,breaks=levels,labels=F)
    test_levels=cut(TEST$A,breaks=levels,labels=F)
  } else{
    train_levels=train$A
    test_levels=TEST$A
    levels<-rep(NA,K+1)->values
  }
  Sub=list()
  Sub$S=train$S
  fits=list()
  for (k in 1:(K-1)){
    if (global){
      sub_id=rep(TRUE,n)
      neg_id=train_levels<=k
      pos_id=train_levels>k
    } else {
      sub_id=train_levels%in%c((k-margin):(k+1+margin))
      neg_id=train_levels%in%c((k-margin):k)
      pos_id=train_levels%in%c((k+1):(k+1+margin))
    }
    tmp=train$A
    tmp[neg_id]=0
    tmp[pos_id]=1
    Sub$A=rbind(Sub$A,matrix(tmp[sub_id],ncol=1))
    tmp=matrix(0,ncol=K-1,nrow=sum(sub_id))
    tmp[,k]=1
    Sub$X=rbind(Sub$X,cbind(train$X[sub_id,],tmp))
    Sub$Y=c(Sub$Y,train$Y[sub_id])
    Sub$propensity=c(Sub$propensity,train$propensity[sub_id])
    # sum((Sub$Y<Sub$S)*(Sub$A==1))
    # sum((Sub$Y>Sub$S)*(Sub$A==0))
    # sum((Sub$Y<Sub$S)*(Sub$A==0))
    # sum((Sub$Y>Sub$S)*(Sub$A==1))
  }
  rownames(Sub$X)=NULL
  Sub$weights=Sub$propensity*(alpha[1]*(Sub$Y<Sub$S)*(Sub$A==1)+alpha[2]*(Sub$Y>Sub$S)*(Sub$A==0))
  svmfit=svm(x=Sub$X,y=factor(Sub$A),w=Sub$weights+0.01,kernel ="linear",probability=TRUE)
  pred=predict.ordinal(svmfit,TEST$X,K)
  tmp=table(pred$labels)
  
  if (!two.sided) {
    output=list(labels=pred$labels,Prediction=pred$Prediction,fit=svmfit,values=values,value=values[pred$labels])
    return(output)
  } else {
    tmp=lapply(split(pred$Prediction,seq(NROW(pred$Prediction))),find.ordinal)
    lb=sapply(tmp,function(x){return(min(x))})
    ub=sapply(tmp,function(x){return(max(x))})
  }
  output=list(labels=list(lb=lb,ub=ub,interval=tmp),Prediction=pred$Prediction,fit=svmfit,values=values,value=list(lvalue=values[lb],uvalue=values[ub]))
  return(output)
}



predict.ordinal.class<-function(fit,X){
  Pred=NULL
  for (k in 1:K){
    tmp=matrix(0,ncol=K,nrow=nrow(X))
    tmp[,k]=1
    pred=predict(fit,cbind(X,tmp))
    Pred=cbind(Pred,as.numeric(pred)-1)
  }
  total=pmax(pmin(apply(Pred,1,sum)+1,K),1)
  return(list(labels=total,Prediction=Pred))
}
class.ordinal<-function(train,TEST,K,two.sided=FALSE,continuous=TRUE,alpha,breaks='quantile',kernel='linear'){
  values=NULL
  n=length(train$A)
  if (continuous){
    if (length(breaks)>1){
    } else if (breaks=='uniform'){
      levels=seq(aL,aU,length.out = K+1)
      for (i in 1:K) values=c(values,mean(levels[i:(i+1)]))
    } else if (breaks=='quantile'){
      levels=quantile(train$A,probs=seq(0,1,length.out = K+1))
      for (i in 1:K) values=c(values,mean(with(train,A[(A<levels[i+1])&(A>=levels[i])])))
    }
    levels[1]=min(train$A)-0.001
    levels[K+1]=max(train$A)+0.001
    train_levels=cut(train$A,breaks=levels,labels=F)
    test_levels=cut(TEST$A,breaks=levels,labels=F)
  } else {
    train_levels=train$A
    test_levels=TEST$A
  }
  
  Sub=list()
  Sub$S=train$S
  fits=list()
  for (k in 1:K){
    sub_id=(train_levels==k)
    tmp=train$Y[sub_id]
    tmp[train$Y[sub_id]>train$S]=0
    tmp[train$Y[sub_id]<train$S]=1
    Sub$Y=rbind(Sub$Y,matrix(tmp,ncol=1))
    tmp=matrix(0,ncol=K,nrow=sum(sub_id))
    tmp[,k]=1
    Sub$X=rbind(Sub$X,cbind(train$X[sub_id,],tmp))
    Sub$propensity=c(Sub$propensity,train$propensity[sub_id])
  }
  Sub$weights=Sub$propensity*(alpha[1]*(Sub$Y==1)+alpha[2]*(Sub$Y==0))
  svmfit=svm(x=Sub$X,y=factor(Sub$Y),w=Sub$weights+0.01,kernel ="linear",probability=TRUE)
  pred=predict.ordinal.class(svmfit,TEST$X)
  pred$labels=pred$labels
  tmp=table(pred$labels)
  
  if (!two.sided) {
    output=list(labels=pred$labels,Prediction=pred$Prediction,fit=svmfit)
    if (continuous){
      output$levels=levels
      output$values=values
      output$value=values[pred$labels]
    }
    return(output)
  } else {
    tmp=lapply(split(pred$Prediction,seq(NROW(pred$Prediction))),find.ordinal)
    lb=sapply(tmp,function(x){return(min(x))})
    ub=sapply(tmp,function(x){return(max(x))})
  }
  output=list(labels=list(lb=lb,ub=ub,interval=tmp),Prediction=pred$Prediction,fit=svmfit)
  if (continuous){
    output$levels=levels
    output$values=values
    output$value=list(lvalue=values[lb],uvalue=values[ub])
  }
  return(output)
}