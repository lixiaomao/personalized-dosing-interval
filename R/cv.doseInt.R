#' Conduct cross-validation for dose interval
#'
#' This function is the wrapper of cross-validation for doseInt
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
#' 
#' 
#' 

cv.doseInt<-function(train,test=NULL,pred0=NULL,nfolds=5,alpha=c(0.5,0.5),type='PDI',two.sided=FALSE,
                     family=c('continuous','ordinal','as.ordinal'),
                     method=c('rq','svmLinear','svmRadial','RLT','tree'),
                     maxiter=20,step=1,trace=0,lower=TRUE,
                     K=10,global=F,margin=0,breaks='quantile',
                     Eps=c(0.1,0.2,0.3),Lambda=c(2^(-(0:10)),0.00000001),Cost=c(0.00001,0.001,0.1,seq(1,15,2)),Embed.mtry=c(0.2,0.3,0.4,0.5),
                     ...){
  #if (is.null(Eps)) Eps=seq(0.1*sd(train$A),0.5*sd(train$A),length.out = 10)
  if (is.null(test)) test=train
  NULL->pred->pred_L->pred_R
  n=length(train$A)
  p=dim(train$X)[2]
  nsize=n/nfolds
  
  Eps=0.1
  grid.cv = switch(method,'rq'=expand.grid(Eps = Eps, Lambda = Lambda),
                   'svmLinear'=expand.grid(Eps = Eps, Cost = Cost),
                   'svmRadial'=expand.grid(Eps = Eps, Cost = Cost),
                   'RLT'=expand.grid(Eps = Eps, Embed.mtry = Embed.mtry))
  
  Rec=rep(0,nrow(grid.cv))
  i=1;k=1
  for(i in 1:nfolds){
    print(paste("=== Fold ",i,"===",sep=" "))
    index = seq(((i-1)*nsize+1),(i*nsize))
    train0 = subsetData(train,-index)
    test0 = subsetData(train,index)
    for(k in 1:dim(grid.cv)[1]){
      tryCatch({
        tmp=doseInt(train=train0,pred0=pred0,
                    alpha=alpha,type=type,
                    family=family,
                    two.sided=two.sided,
                    method=method, eps=grid.cv[k,1],
                    maxiter=maxiter,step=step,trace=trace,lower=lower,lambda=grid.cv[k,2],cost=grid.cv[k,2],embed.mtry = grid.cv[k,2],
                    K=K,global=global,margin=margin)[[1]]
        Rec[k]= Rec[k]+tmp$misclass
      },error=function(e){print("Error happens in cv.owlk");return(NA)})
    }
  }
  Rec=Rec/nfolds
  Rec=cbind(Rec,grid.cv)
  names(Rec)=c("Misclass","Eps","Para")
  id=which.min(Rec$Misclass)
  print(Rec[id,])
  misclass_min=Rec$Misclass[id]
  
  # independent test
  fit=doseInt(train=train,
              alpha=alpha,type=type,
              family=family,
              two.sided=two.sided,
              method=method, eps=Rec[id,2],
              maxiter=maxiter,step=step,trace=trace,lower=lower,lambda=Rec[id,3],cost=Rec[id,3],embed.mtry = Rec[id,3],
              K=K,global=global,margin=margin)[[1]]
  tmp=predict(fit,test)
  misclass=tmp$misclass
  region=tmp$region
  percentage=sum(region*(test$Y>test$S))/sum(region)
  
  output=list(misclass_min=misclass_min,misclass=misclass,fit=fit,Rec=Rec,bestPara=Rec[id,3],percentage=percentage)
  if(!two.sided){
    output$pred=tmp$pred
    if (!is.null(test$opt)){
      if (family=='continuous'){
        output$Rsquare=1-mean((output$pred-test$opt)^2)/var(test$opt)
        output$correlation=cor(output$pred,test$opt)
      } else {
        output$value=tmp$value
        output$Rsquare=1-mean((output$value-test$opt)^2)/var(test$opt)
        output$correlation=cor(output$value,test$opt)
      }
    }
    
  } else {
    output$pred_L=tmp$pred_L
    output$pred_R=tmp$pred_R
    if ((!is.null(test$opt_L))&(!is.null(test$opt_R))){
      if (family=='continuous'){
        Rsquare=1-mean((output$pred_L-test$opt_L)^2)/var(test$opt_L)
        Rsquare=Rsquare+1-mean((output$pred_R-test$opt_R)^2)/var(test$opt_R)
        output$Rsquare=Rsquare/2
        output$correlation=0.5*cor(output$pred_L,test$opt_L)+0.5*cor(output$pred_R,test$opt_R)
      } else {
        output$value_L=tmp$value_L
        output$value_R=tmp$value_R
        Rsquare=1-mean((output$value_L-test$opt_L)^2)/var(test$opt_L)
        Rsquare=Rsquare+1-mean((output$value_R-test$opt_R)^2)/var(test$opt_R)
        output$Rsquare=Rsquare/2
        output$correlation=0.5*cor(output$value_L,test$opt_L)+0.5*cor(output$value_R,test$opt_R)
      }
    }
  }
  return(output)
}