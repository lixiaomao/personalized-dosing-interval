#' Fits dose interval
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
#' @param maxiter maximal iterations for DC algorithm
#' @param step stepsize for DC algorithm
#' @param trace the level of detail to be printed: 1,2,3,4
#' @param lower True or False, Whether lower boundary is considered when applied to one-sided interval
#' @param K level of ordinal treatment, default 10
#' @param global if global data are used for ordinal approach, default False
#' @param margin how large is the margin if global=False, default 0
#' @param breaks ways 'as.ordinal' cuts continuous treatment into several ordinal levels, default 'quantile'. User can also choose 'uniform' or specify customized breaks
#' @param Eps parameter for the relaxation
#' @param Lambda  penalty for the quantile regression
#' @param Cost  cost for the support vector machine
#' @param Embed.mtry  proportion of mtry
#'
#' @examples
#' train1=Scenario2.continuous(2000,50,1)
#' test1=Scenario2.continuous(10000,50,1)
#' fit=doseInt(train=train1,
#'             alpha=c(0.5,0.5),type='PDI',
#'             family='continuous',
#'             two.sided=FALSE,
#'             method=c('rq','svmRadial'),
#'             maxiter=20,step=1,trace=3,lower=TRUE,lambda=0.000001,
#'             K=10,global=F,margin=0)
#' pred=predict(fit,test1)
#' @export
#'


doseInt<-function(y=NULL,x=NULL,a=NULL,s=NULL,propensity=NULL,train=NULL,
                  alpha=c(0.5,0.5),type='PDI',
                  family=c('continuous','ordinal','as.ordinal'),
                  two.sided=FALSE,
                  method=c('rq','svmLinear','svmRadial','RLT','tree'),
                  pred0=NULL,
                  maxiter=20,step=1,trace=0,lower=TRUE,lambda=NULL,cost=1,embed.mtry = 1/2,
                  K=10,global=F,margin=0,breaks='quantile',
                  ...){
  if (!is.null(family)){
    family=match.arg(family)
  } else if (length(unique(a))>20) {
    family='continuous'
  } else {
    family='ordinal'
  }
  if (is.null(train)) train=constructData(y,x,a,s,propensity)

  output=list()
  if (!is.null(method)) method=match.arg(method,several.ok = T)

  if (family=='continuous'){
    for (m in seq_along(method)){
      if (!two.sided){
        output[[length(output)+1]]=DCDI1(train,alpha,method=method[m],type=type,pred0=pred0,
                                         maxiter=maxiter,step=step,trace=trace,lower=lower,lambda=lambda,cost=cost,embed.mtry=embed.mtry)

      } else{
        output[[length(output)+1]]=DCDI2(train,alpha,method=method[m],type=type,pred0=pred0,
                                         maxiter=maxiter,step=step,trace=trace,lambda=lambda,cost=cost,embed.mtry=embed.mtry)
      }
    }
  } else if (family=='ordinal'|family=='as.ordinal'){
    if (is.null(method)) method='svmRadial'
    for (m in seq_along(method)){
      if (!method[m]%in%c('svmLinear','svmRadial')){
        cat('Not able to use ',method[m],' for ordinal treatment\n')
        next
      }
      if (family=='as.ordinal') continuous=TRUE else continuous=FALSE
      output[[length(output)+1]]=ordinalDI(train,type=type,continuous=continuous,two.sided=two.sided,lower=lower,
                                           K=K,cost=cost,global=global,margin=margin,alpha=alpha,breaks=breaks,method=method[m])
    }
  }
  names(output)=method
  class(output)='doseInt'
  return(output)
}

predict.doseInt<-function(x,test){
  lst=list()
  for (i in seq_along(x)) lst[[i]]=predict(x[[i]],test)
  return(lst)
}
