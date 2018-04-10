
# initiate.continuous<-function(train,test=NULL,two.sided,lambda=0.1,...){
#   if (is.null(test)) test=train
#   n=length(train$A)
#   mytrain=cbind(train$Y,train$A,train$X,1:length(train$A))
#   colnames(mytrain)=c("Y","A",paste("X",1:dim(train$X)[2],sep=""),"id")
#   mytrain=as.data.frame(mytrain)
#   mytest=as.data.frame(test$X)
#   colnames(mytest)=c(paste("X",1:dim(train$X)[2],sep=""))
#   if (two.sided==FALSE){
#     index1=(train$Y<quantile(train$Y,min(sum(train$Y<train$S)/n+0.15,1)))&(train$Y>(quantile(train$Y,max(sum(train$Y<train$S)/n-0.15,0))))
#     id2=mytrain$id[index1][(train$A[index1]>quantile(train$A[index1],0.35))&(train$A[index1]<quantile(train$A[index1],0.95))]
#     lmfit=rq(A~.,data=mytrain[id2,c("A",paste("X",1:dim(train$X)[2],sep=""))],method="lasso",lambda=lambda)
#     pred=predict(lmfit,mytest)
#     return(list(fit=lmfit,pred=pred))
#   } else if (two.sided==TRUE){
#     index1=(train$Y<quantile(train$Y,min(sum(train$Y<train$S)/n+0.15,1)))&(train$Y>(quantile(train$Y,max(sum(train$Y<train$S)/n-0.15,0))))
#     id1=mytrain$id[index1][(train$A[index1]>quantile(train$A[index1],0.05))&(train$A[index1]<quantile(train$A[index1],0.45))]
#     id2=mytrain$id[index1][(train$A[index1]>quantile(train$A[index1],0.55))&(train$A[index1]<quantile(train$A[index1],0.95))]
#     model1=rq(A~.,data=mytrain[id1,!colnames(mytrain)%in%c('Y','id')],method='lasso',lambda=lambda)
#     model2=rq(A~.,data=mytrain[id2,!colnames(mytrain)%in%c('Y','id')],method='lasso',lambda=lambda)
#     h1=predict(model1,mytest)
#     h2=predict(model2,mytest)
#     fit=list(fit_L=model1,fit_R=model2)
#     return(list(fit=fit,pred_L=h1,pred_R=h2))
#   }
# }
# 
# initiate.binary<-function(train,test=NULL,two.sided,lambda=0.1,...){
#   if (is.null(test)) test=train
#   n=length(train$A)
#   mytrain=as.data.frame(cbind(train$A,train$X))
#   colnames(mytrain)=c("A",paste("X",1:dim(train$X)[2],sep=""))
#   mytest=as.data.frame(train$X)
#   colnames(mytest)=c(paste("X",1:dim(test$X)[2],sep=""))
#   if (two.sided==FALSE){
#     index=(train$A<quantile(train$A,0.65))&(train$A>quantile(train$A,0.35))
#     lmfit=rq(A~.,data=mytrain[index,],method="lasso",lambda=lambda)
#     pred=predict(lmfit,mytest)
#     return(list(pred_test=pred,fit=lmfit))
#   } else if (two.sided==TRUE){
#     id1=(mytrain$A>quantile(mytrain$A,0.2))&(mytrain$A<quantile(mytrain$A,0.35))
#     id2=(mytrain$A>quantile(mytrain$A,0.55))&(mytrain$A<quantile(mytrain$A,0.7))
#     model1=rq(A~.,data=mytrain[id1,],method="lasso",lambda=lambda)
#     model2=rq(A~.,data=mytrain[id2,],method="lasso",lambda=lambda)
#     h1=predict(model1,mytest)
#     h2=predict(model2,mytest)
#     fit=list(fit_L=model1,fit_R=model2)
#     return(list(fit=fit,pred_L=h1,pred_R=h2))
#   }
# }
#' @export
initiate.continuous<-function(train,test=NULL,pred_L=NULL,pred_R=NULL,pred=NULL,two.sided,lower=NULL,type,...){
  if ((!is.null(pred))&(!two.sided)){
    output=list(pred=median(pred),two.sided=FALSE,lower=lower,type=type)
    class(output)='init'
    return(output)
  }
  if ((!is.null(pred_L))&two.sided){
    output=list(pred_L=median(pred_L),pred_R=median(pred_R),two.sided=TRUE,lower=lower,type=type)
    class(output)='init'
    return(output)
  }
  fn<-function(par,train,two.sided,mylower,type){
    if (two.sided==FALSE){
      region=par[1]<train$A
      if (type=="EDI"){
        if (lower){
          misclass=(sum(train$alpha[,1]*(region)*pmax((train$S-train$Y),0))+sum(train$alpha[,2]*(!region)*pmax((train$Y-train$S),0)))/length(train$A)
        } else{
          misclass=(sum(train$alpha[,1]*(!region)*pmax((train$S-train$Y),0))+sum(train$alpha[,2]*(region)*pmax((train$Y-train$S),0)))/length(train$A)
        }
      } else if (type=="PDI"){
        if (mylower){
          misclass=(sum(train$alpha[,1]*(region)*(train$Y<train$S))+sum(train$alpha[,2]*(!region)*(train$Y>train$S)))/length(train$A)
        } else {
          misclass=(sum(train$alpha[,1]*(!region)*(train$Y<train$S))+sum(train$alpha[,2]*(region)*(train$Y>train$S)))/length(train$A)
        }
      }
    } else if (two.sided==TRUE){
      region=(par[2]<train$A)&(par[3]>train$A)
      if (type=="EDI"){
        misclass=(sum(train$alpha[,1]*(region)*pmax((train$S-train$Y),0))+sum(train$alpha[,2]*(!region)*pmax((train$Y-train$S),0)))/length(train$A)
      } else if (type=="PDI"){
        misclass=(sum(train$alpha[,1]*(region)*(train$Y<train$S))+sum(train$alpha[,2]*(!region)*(train$Y>train$S)))/length(train$A)
      }
    }
    return(misclass+999*max(par[2]-par[3],0))
  }
  fit=optim(par=quantile(train$A,c(0.5,0.01,0.99)),fn=fn,method="BFGS",lower=rep(quantile(train$A,0.1),3),upper=rep(quantile(train$A,0.9),3),
            train=train,mylower=lower,type=type,two.sided=two.sided)
  output=list(pred=fit$par[1],pred_L=fit$par[2],pred_R=fit$par[3],two.sided=two.sided,lower=lower,type=type)
  class(output)='init'
  return(output)
}

#' @export
initiate.binary=initiate.continuous
########################################################
# function for linear region, linear function
#' @export
fit.rq.model<-function(train,weights,lambda=NULL,...){
  n=length(train$A)
  p=dim(train$X)[2]
  nzero=which(weights>0)
  weights[nzero]=weights[nzero]/mean(weights[nzero])
  if (length(train$A[nzero])<=p){
    cat(paste0('influntial sample:',length(train$A[nzero]),' parameters to estimate:',p+1))
    return(NA);
  }
  rq.model=rq.lasso.fit(x=train$X[nzero,]*weights[nzero],y=train$A[nzero]*weights[nzero],tau=.5,lambda=lambda)
  return(rq.model)
}
#' @export
fit.rlt.model<-function(train,weights,...){
  n=length(train$A)
  p=dim(train$X)[2]
  nzero=which(weights>0)
  weights[nzero]=weights[nzero]/mean(weights[nzero])
  if (length(train$A[nzero])<=p & FALSE){
    cat(paste0('influntial sample:',length(train$A[nzero]),' parameters to estimate:',p+1))
    return(NA);
  }
  RLT.fit = RLT(x=as.matrix(train$X[nzero,]), y=train$A[nzero], model = "regression", use.cores = 2, ntrees = 200, split.gen = "random",
                nsplit = 1, resample.prob = 1, track.obs = FALSE, importance = TRUE,
                replacement = TRUE, reinforcement = TRUE, combsplit = 5, embed.ntrees = 25,
                muting = -1, muting.percent = 0.9, embed.mtry = 2/3, subject.weight=weights[nzero])
  return(RLT.fit)
}


#' @export
fit.svmL.model<-function(train,weights,cost=1,...){
  n=length(train$A)
  svm.model=NA
  nonzero=which(weights>0)
  tryCatch({
    svm.model  = svm(x = train$X[nonzero,], y = (train$A[nonzero]), w= weights[nonzero], type="eps-regression",kernel = "linear",cost=cost,epsilon = 0,scale=TRUE)
  },error=function(e){svm.model=NA},finally={})
  return(svm.model)
}
#' @export
fit.svmK.model<-function(train,weights,cost=1,...){
  n=length(train$A)
  svm.model=NA
  nonzero=which(weights>0)
  tryCatch({
    svm.model  = svm(x = train$X[nonzero,], y = (train$A[nonzero]), w= weights[nonzero], type="eps-regression",kernel = "radial",cost=cost,epsilon = 0,scale=TRUE)
  },error=function(e){svm.model=NA},finally={})
  return(svm.model)
}
#' @export
fit.LK.model<-function(train,weights,...){
  n=length(train$A)
  svm.model=NA
  para=list(...)
  if (!is.null(para$eps)) eps<-para$eps else eps=0.05*(aU-aL)
  if (!is.null(para$gamma)) gamma<-para$gamma else gamma=NULL
  gamma=min(gamma,0.999999)
  nzero=which(weights>0)
  tryCatch({
    svm.model  = ksvm(train$X[nzero,],train$A[nzero],C=-log(gamma),epsilon=0,type="eps-svr",kernel = "rbfdot",kpar="automatic")
    #print(svm.model@kernelf@kpar)
  },error=function(e){svm.model=NA},finally={})
  return(svm.model)
}



##############################################
#' @export
predict.owl <- function(model,test){
  if (sum(class(model)=="svm")>0){
    pred=predict(model,test$X)
  } else if (sum(class(model)=="rq.pen")>0){
    pred = as.vector(predict(model,test$X))
  } else if (sum(grep('rq',model$call)) ==1){
    pred = model$coefficients[1] + test$X %*% model$coefficients[-1]
  } else if(sum(grep('RLT',model$call))==0){
    pred = predict(model,test$X)$Prediction
  } else{
    pred = predict(model,test$X)
  }
  return(pred)
}

# predict.init<-function(fit,test){
#   newdata=as.data.frame(test$X)
#   colnames(newdata)=c(paste("X",1:dim(test$X)[2],sep=""))
#   pred=predict(fit,newdata)
#   return(pred)
# }
#' @export
predict.init<-function(fit,test){
  if (!fit$two.sided){
    pred=rep(fit$pred,length(test$A))
  } else{
    pred=list(pred_L=rep(fit$pred_L,length(test$A)),pred_R=rep(fit$pred_R,length(test$A)))
  }
  return(pred)
}



# fit=initiate.continuous(train,test=NULL,pred=NULL,two.sided=two.sided,lower=lower,type=type)
# pp=predict(fit,test)
