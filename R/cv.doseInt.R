#' Conduct cross-validation for dose interval
#'
#' This function is the wrapper of cross-validation for doseInt
#'
#' @param train the training data set.
#' @param test the testing data set. If test is NULL, use train as test.
#' @param pred0 the initial prediction for training data set, if NULL, default initialization is used
#' @param nfolds the number of folds used in cross-validation.
#' @param type indicates whether PDI or EDI is to be calculated.
#' @param two.sided indicator of whether two-sided interval is considered.
#' @param family specifis the methods. 'continuous' leads to DC-algorithm,'ordinal' handles the ordinal treatments, while 'as.ordinal' cuts continuous treatment into several ordinal levels and apply the ordinal algorithm.
#' @param method specified methods used for DC and ordinal approaches.
#' @param maxiter maximal iterations for DC algorithm
#' @param step stepsize for DC algorithm
#' @param trace the level of detail to be printed: 1,2,3,4
#' @param lower True or False, Whether lower boundary is considered when applied to one-sided interval
#' @param K level of ordinal treatment, default 10
#' @param global if global data are used for ordinal approach, default False
#' @param margin how large is the margin if global=False, default 0
#' @param breaks ways 'as.ordinal' cuts continuous treatment into several ordinal levels, default 'quantile'. User can also choose 'uniform' or specify customized breaks
#' @param Eps parameters for the relaxation
#' @param Lambda candidate penalties for the quantile regression
#' @param Cost candidate costs for the support vector machine
#' @param Embed.mtry candidate proportions of mtry
#'
#' @examples
#' train1=Scenario1.continuous(2000,5,1)
#' test1=Scenario1.continuous(10000,5,1)
#'
#' cv1=cv.doseInt(train1,test1,pred0=NULL,nfolds=2,alpha=c(0.5,0.5),type='PDI',two.sided=FALSE,
#'                family='continuous',
#'                method='svmLinear',
#'                maxiter=20,step=1,trace=0,lower=TRUE,Lambda = 2^(-(1:5)),Cost=c(1,7,12,19))
#' pred1=predict(cv1$fit,test1)
#'
#' train2=Scenario2.continuous(2000,6,2)
#' test2=Scenario2.continuous(10000,6,2)
#'
#' cv2=cv.doseInt(train2,test2,pred0=NULL,nfolds=2,alpha=c(0.5,0.5),type='PDI',two.sided=FALSE,
#'                family='as.ordinal',
#'                method='svmLinear',K=20,
#'                maxiter=2,step=1,trace=3,lower=TRUE,Lambda = 2^(-(1:5)),Cost=c(0.0001,0.001,0.01,0.01,0.05,0.1,0.5,1))
#' pred2=predict(cv2$fit,test2)
#'
#'
#' train4=Scenario4.continuous(2000,5,4)
#' test4=Scenario4.continuous(10000,5,4)
#'
#' cv4=cv.doseInt(train4,test4,pred0=NULL,nfolds=2,alpha=c(0.5,0.5),type='PDI',two.sided=TRUE,
#'                family='as.ordinal',
#'                method='svmLinear',K=20,
#'               maxiter=2,step=1,trace=3,lower=TRUE,Lambda = 2^(-(1:5)),Cost=c(0.0001,0.001,0.01,0.01,0.05,0.1,0.5,1))
#' pred4=predict(cv4$fit,test4)
#' @export

cv.doseInt<-function(train,test=NULL,pred0=NULL,nfolds=5,type='PDI',two.sided=FALSE,
                     family=c('continuous','ordinal','as.ordinal'),
                     method=c('rq','svmLinear','svmRadial','RLT','tree','completeLinear','completeKernel'),
                     maxiter=20,step=1,trace=0,lower=TRUE,
                     K=10,global=F,margin=0,breaks='quantile',aL=NULL,aU=NULL,
                     Eps=NULL,Lambda=c(2^(-(0:10)),0.00000001),Cost=c(0.00001,0.001,0.1,seq(1,15,2)),Embed.mtry=c(0.2,0.3,0.4,0.5),
                     ...){
  #if (is.null(Eps)) Eps=seq(0.1*sd(train$A),0.5*sd(train$A),length.out = 10)
  
  if (is.null(test)) test=train
  NULL->pred->pred_L->pred_R
  n=length(train$A)
  p=dim(train$X)[2]
  nsize=n/nfolds

  if (is.null(Eps)) Eps=0.1*(max(train$A)-min(train$A))
  grid.cv = switch(method,'rq'=expand.grid(Eps = Eps, Lambda = Lambda),
                   'svmLinear'=expand.grid(Eps = Eps, Cost = Cost),
                   'svmRadial'=expand.grid(Eps = Eps, Cost = Cost),
                   'RLT'=expand.grid(Eps = Eps, Embed.mtry = Embed.mtry),
                   'completeLinear'=expand.grid(Eps = Eps, Lambda = Lambda),
                   'completeKernel'=expand.grid(Eps = Eps, Lambda = Lambda))

  Rec=rep(0,nrow(grid.cv))
  i=1;k=1
  for(i in 1:nfolds){
    print(paste("=== Fold ",i,"===",sep=" "))
    index = seq(((i-1)*nsize+1),(i*nsize))
    train0 = subsetData(train,-index)
    test0 = subsetData(train,index)
    for(k in 1:dim(grid.cv)[1]){
      print(grid.cv[k,])
      #tryCatch({
        #fit.Int1(train=train0, eps = eps, alpha , pred0 = 0, lambda = lambda,
                 #type = type, lower = lower)
        #aa=RDCDI1(train=train0, alpha=alpha, type = "PDI", method = method, pred0 = 0,
         #      maxiter = 25, lower = lower, eps = eps, lambda = lambda)
        fit=doseInt(train=train0,pred0=pred0,
                    type=type,
                    family=family,
                    two.sided=two.sided, aL=aL,aU=aU,
                    method=method, eps=grid.cv[k,1],
                    maxiter=maxiter,step=step,trace=trace,lower=lower,lambda=grid.cv[k,2],cost=grid.cv[k,2],embed.mtry = grid.cv[k,2],
                    K=K,global=global,margin=margin)
        tmp=predict(fit,test0)
        Rec[k]= Rec[k]+tmp$misclass
      #},error=function(e){print(e);Rec[k]= Rec[k]+NA;})
    }
  }
  Rec=Rec/nfolds
  Rec=cbind(Rec,grid.cv)
  names(Rec)=c("Misclass","Eps","Para")
  id=which.min(Rec$Misclass)
  print(Rec[id,])
  misclass_min=Rec$Misclass[id]

  # independent test
  fit=doseInt(train=train,pred0=pred0,
              type=type,
              family=family,
              two.sided=two.sided,
              method=method, eps=Rec[id,2],aL=aL,aU=aU,
              maxiter=maxiter,step=step,trace=trace,lower=lower,lambda=Rec[id,3],cost=Rec[id,3],embed.mtry = Rec[id,3],
              K=K,global=global,margin=margin)
  tmp=predict(fit,test)
  measures=measurement(tmp,test,two.sided=two.sided,family=family,lower=lower)
  output=list(misclass_min=misclass_min,measures=measures,fit=fit,Rec=Rec,bestPara=Rec[id,3])
  return(output)
}

#' @export
measurement=function(obj,test,two.sided,family,lower=NULL){

  n=length(test$A)
  NA->correlation->Rsquare->pred->pred_L->pred_R->misclass->cost->percentage
  if(!two.sided){
    pred=obj$pred
    if (lower){
      region=(test$A>pred)
    } else {
      region=(test$A<pred)
    }
    if (!is.null(test$opt)){
      if (family=='continuous'){
        Rsquare=1-mean((pred-test$opt)^2)/var(test$opt)
        correlation=cor(drop(pred),test$opt)
      } else {
        value=obj$value
        Rsquare=1-mean((value-test$opt)^2)/var(test$opt)
        correlation=cor(value,test$opt)
      }
    }
  } else {
    pred_L=obj$pred_L
    pred_R=obj$pred_R
    region=(test$A>pred_L)&(test$A<pred_R)
    if ((!is.null(test$opt_L))&(!is.null(test$opt_R))){
      if (family=='continuous'){
        Rsquare=1-mean((pred_L-test$opt_L)^2)/var(test$opt_L)
        Rsquare=Rsquare+1-mean((pred_R-test$opt_R)^2)/var(test$opt_R)
        Rsquare=Rsquare/2
        correlation=0.5*cor(drop(pred_L),test$opt_L)+0.5*cor(drop(pred_R),test$opt_R)
      } else {
        value_L=obj$value_L
        value_R=obj$value_R
        Rsquare=1-mean((value_L-test$opt_L)^2)/var(test$opt_L)
        Rsquare=Rsquare+1-mean((value_R-test$opt_R)^2)/var(test$opt_R)
        Rsquare=Rsquare/2
        correlation=0.5*cor(value_L,test$opt_L)+0.5*cor(value_R,test$opt_R)
      }
    }
  }
  misclass=(sum(test$alpha[,1]*(region)*(test$Y<test$S))+sum(test$alpha[,2]*(!region)*(test$Y>test$S)))/length(test$A)
  cost=(sum(test$alpha[,1]*(region)*pmax((test$S-test$Y),0))+sum(test$alpha[,2]*(!region)*pmax((test$Y-test$S),0)))/length(test$A)
  percentage=sum(region*(test$Y>test$S))/sum(region)

  return(list(correlation=correlation,Rsquare=Rsquare,pred=pred,pred_L=pred_L,pred_R=pred_R,misclass=misclass,cost=cost,percentage=percentage))
}
