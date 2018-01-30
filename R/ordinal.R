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
#'   fit1=ordinalDI(train1,alpha=c(0.5,0.5),K=20,continuous = T)
#'   predNew1=predict(fit1,test1)
#'   1-mean((predNew1$value-test1$opt)^2)/var(test1$opt)
#'
#'   fit3=ordinalDI(train3,two.sided=T,alpha=c(0.5,0.5),K=20,continuous=T)
#'   predNew3=predict(fit3,test3)
#'   cor(predNew3$value_L,test3$opt_L)^2
#'   cor(predNew3$value_R,test3$opt_R)^2
#'

ordinalDI<-function(train,type='PDI',two.sided=FALSE,cost=1,K=10,continuous=FALSE,
                    global=F,margin=0,alpha=c(0.5,0.5),breaks='quantile',lower=TRUE,method='svmLinear',...){
  levels=NULL
  n=length(train$A)
  if(!is.null(train$K)) K=train$K
  if (continuous){
    if (length(breaks)>1){
      breaks=breaks
      for (i in 1:K) levels=c(levels,mean(breaks[i:(i+1)]))
    } else if (breaks=='uniform'){
      breaks=seq(aL,aU,length.out = K+1)
      for (i in 1:K) levels=c(levels,mean(breaks[i:(i+1)]))
    } else if (breaks=='quantile'){
      breaks=quantile(train$A,probs=seq(0,1,length.out = K+1))
      for (i in 1:K) levels=c(levels,mean(with(train,A[(A<breaks[i+1])&(A>=breaks[i])])))
    }
    breaks[1]=aL
    breaks[K+1]=aU
    train_breaks=cut(train$A,breaks=breaks,labels=F)
  } else{
    train_breaks=train$A
    breaks<-rep(NA,K+1)->levels
  }
  Sub=list()
  Sub$S=train$S
  fits=list()
  for (k in 1:(K-1)){
    if (global){
      sub_id=rep(TRUE,n)
      neg_id=train_breaks<=k
      pos_id=train_breaks>k
    } else {
      sub_id=train_breaks%in%c((k-margin):(k+1+margin))
      neg_id=train_breaks%in%c((k-margin):k)
      pos_id=train_breaks%in%c((k+1):(k+1+margin))
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
  }
  if (type=='PDI'){
    Sub$weights=Sub$propensity*(alpha[1]*(Sub$Y<Sub$S)*(Sub$A==1)+alpha[2]*(Sub$Y>Sub$S)*(Sub$A==0))
  } else if (type=='EDI'){
    Sub$weights=Sub$propensity*(alpha[1]*pmax(0,Sub$S-Sub$Y)*(Sub$A==1)+alpha[2]*pmax(0,Sub$Y-Sub$S)*(Sub$A==0))
  }

  kernel=switch(method,'svmLinear'='linear','svmRadial'='radial')
  index=Sub$weights>0
  svmfit=svm(x=Sub$X[index,],y=factor(Sub$A)[index],w=Sub$weights[index],cost=cost,kernel =kernel,probability=TRUE)
  svmpred=predict.ordinal(svmfit,train$X,K)

  if (!two.sided) {
    svmpred$pred[is.na(svmpred$pred)]=median(svmpred$pred,na.rm = TRUE)
    output=list(two.sided=FALSE,pred=svmpred$pred,value=levels[svmpred$pred])
    class(output)='ordinalDI1'

    if ((!is.null(train$Y))&(!is.null(train$S))&(!is.null(train$A))){
      if (lower) region=output$value<=train$A else region=output$value>=train$A
      misclass=mean(train$propensity*(alpha[1]*(train$Y<train$S)*(region)+alpha[2]*(train$Y>train$S)*(!region)))
    }
  } else {
    tmp=lapply(split(svmpred$Prediction,seq(NROW(svmpred$Prediction))),find.ordinal)
    lb=sapply(tmp,function(x){return(min(x))})
    lb[is.na(lb)]=median(lb,na.rm = TRUE)
    ub=sapply(tmp,function(x){return(max(x))})
    ub[is.na(ub)]=median(ub,na.rm = TRUE)
    ub=pmin(K,ub+1)
    output=list(K=K,two.sided=TRUE,pred_L=lb,pred_R=ub,value_L=levels[lb],value_R=levels[ub])
    class(output)='ordinalDI2'

    if ((!is.null(train$Y))&(!is.null(train$S))&(!is.null(train$A))){
      region=(output$value_L<=train$A)&(output$value_R>=train$A)
      misclass=mean(train$propensity*(alpha[1]*(train$Y<train$S)*(region)+alpha[2]*(train$Y>train$S)*(!region)))
    }
  }
  output$misclass=misclass
  output$alpha=alpha
  output$K=K
  output$Prediction=svmpred$Prediction
  output$fit=svmfit
  output$levels=levels
  output$region=region
  return(output)
}


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
    out=predict(fit,cbind(X,tmp))
    Pred=cbind(Pred,as.numeric(out)-1)
  }
  pred=apply(Pred,1,sum)+1
  return(list(pred=pred,Prediction=Pred))
}
predict.ordinalDI1<-function(obj,test,K,lower=TRUE){
  two.sided=obj$two.sided
  K=obj$K
  levels=obj$levels
  fit=obj$fit
  alpha=obj$alpha

  Prediction=NULL
  for (k in 1:(K-1)){
    tmp=matrix(0,ncol=K-1,nrow=nrow(test$X))
    tmp[,k]=1
    out=predict(fit,cbind(test$X,tmp))
    Prediction=cbind(Prediction,as.numeric(out)-1)
  }
  pred=apply(Prediction,1,sum)+1
  pred[is.na(pred)]=median(pred,na.rm = TRUE)
  output=list(pred=pred,Prediction=Prediction,value=levels[pred])

  if ((!is.null(test$Y))&(!is.null(test$S))&(!is.null(test$A))){
    if (lower) region=output$value<=test$A
    misclass=mean(test$propensity*(alpha[1]*(test$Y<test$S)*(region)+alpha[2]*(test$Y>test$S)*(!region)))
    output$misclass=misclass
    output$region=region
  }
  return(output)
}

predict.ordinalDI2<-function(obj,test,K){
  two.sided=obj$two.sided
  K=obj$K
  levels=obj$levels
  fit=obj$fit
  alpha=obj$alpha

  Prediction=matrix(0,nrow=nrow(test$X),ncol=K-1)
  for (k in 1:(K-1)){
    tmp=matrix(0,ncol=K-1,nrow=nrow(test$X))
    tmp[,k]=1
    out=predict(fit,cbind(test$X,tmp))
    Prediction[,k]=as.numeric(out)-1
  }
  pred=apply(Prediction,1,sum)+1
  tmp=lapply(split(Prediction,seq(NROW(Prediction))),find.ordinal)
  lb=sapply(tmp,function(x){return(min(x))})
  lb[is.na(lb)]=median(lb,na.rm = TRUE)
  ub=sapply(tmp,function(x){return(max(x))})
  ub[is.na(ub)]=median(ub,na.rm = TRUE)
  output=list(pred_L=lb,pred_R=ub,Prediction=Prediction,value_L=levels[lb],value_R=levels[ub])

  if ((!is.null(test$Y))&(!is.null(test$S))&(!is.null(test$A))){
    region=(output$value_L<=test$A)&(output$value_R>=test$A)
    misclass=mean(test$propensity*(alpha[1]*(test$Y<test$S)*(region)+alpha[2]*(test$Y>test$S)*(!region)))
    output$misclass=misclass
    output$region=region
  }
  return(output)
}
