#' @export
constructData<-function(y,x,a,s,alpha=NULL,propensity=NULL,K=20,Extra=NULL,grid=NULL,stablizer=0.01,noise=0,numerator=FALSE,type='multinomial'){
  Data=list()
  Data$X=as.matrix(x)
  Data$Y=drop(y)
  Data$A=matrix(a,ncol=1)
  Data$S=s
  if (is.null(alpha)){
    alpha=c(0.5,0.5)
  }
  if (is.vector(alpha)){
    alpha=matrix(rep(alpha,length(a)),ncol=2,byrow=TRUE)
  }
  Data$alpha=alpha
  if (is.null(propensity)) {
    if (noise>0) a=a+runif(length(a),-noise,noise)
    if (length(unique(a))>20){
      if (is.null(grid)){
        grid=seq(min(a),max(a),length.out = K)
        grid[1]=-Inf;grid[length(grid)]=Inf
      }
      a=cut(a,breaks=grid,labels = FALSE,include.lowest = TRUE)
    } else{
      a=as.character(a)
    }
    A_prob=table(a)/length(a)
    dat=as.data.frame(cbind(a,x))
    names(dat)=c('a',paste0('x',1:dim(x)[2]))
    if (type=='multinomial'){

      glmfit=nnet::multinom(a~.,data=dat)
      propensity=rep(0,length(a))
      for (i in seq_along(propensity)){
        if (numerator){
          propensity[i]=A_prob[a[i]]/(glmfit$fitted.values[i,a[i]]+stablizer)
        } else{
          propensity[i]=1/(glmfit$fitted.values[i,a[i]]+stablizer)
        }
        
      }
    } else if (type=='ordinal'){
      fit=orm(a~.,data=dat)
      cumprob=predict(fit,type='fitted')
      cumprob=cbind(cumprob,rep(0,dim(cumprob)[1]))
      propensity=rep(0,length(a))
      for (i in seq_along(propensity)){
        tmp=ifelse(a[i]>1,cumprob[i,a[i]-1]-cumprob[i,a[i]],1-cumprob[i,a[i]])
        if (numerator){
          propensity[i]=A_prob[a[i]]/(tmp+stablizer)
        } else{
          propensity[i]=1/(tmp+stablizer)
        }
      }
    } else {
      print('type not supported')
    }
  }
  Data$propensity=propensity#/mean(propensity)
  if (!is.null(Extra)){
    for (i in seq_along(Extra)){
      Data[names(Extra)[i]]=Extra[i]
    }
  }
  return(Data)
}

#' @export
subsetData<-function(Data,index){
  NewData=list()
  NewData$X=Data$X[index,]
  NewData$alpha=Data$alpha[index,]
  if (!is.null(Data$A)) NewData$A=matrix(Data$A[index],ncol=1)
  for (i in seq_along(Data)){
    if (!names(Data)[i]%in%c('X','A','alpha')) {
      if (length(Data[[i]])==length(Data$A)){
        NewData[[length(NewData)+1]]=Data[[i]][index]
      } else {
        NewData[[length(NewData)+1]]=Data[[i]]
      }
      names(NewData)[length(NewData)]=names(Data)[i]
    }
  }
  return(NewData)
}

#' @export
plot.FDI<-function(Data,pos,size=1,main=NULL,saving=FALSE,...){
  if (saving){
    jpeg(filename = paste0(main,".jpg"),width = 700, height = 700, quality = 99)
  }
  Dat1=data.frame(treatment=jitter(Data$A[which(!pos==1)]),response=jitter(Data$Y[which(!pos==1)]),label_pred=rep('in',sum(!pos==1)))
  Dat2=data.frame(treatment=jitter(Data$A[which(pos==1)]),response=jitter(Data$Y[which(pos==1)]),label_pred=rep('out',sum(pos==1)))
  Dat=rbind(Dat1,Dat2)
  p<-ggplot(data=Dat,aes(x=treatment,y=response,group=label_pred))+geom_point(aes(colour = label_pred),size=1.5,alpha =0.4)+geom_hline(yintercept=Data$S)+ggtitle(main)+theme(plot.title = element_text(hjust = 0.5))
  ggsave(filename=paste0(main,".jpg"), device="jpeg",width = 40, height = 20, units = "cm",plot = p, dpi =1000)
  #plot(Data$A,Data$Y,main=main,pch=".")
  #abline(h=Data$S)
  #points(jitter(Data$A[which(!pos==1)]),jitter(Data$Y[which(!pos==1)]),col="blue",pch=19,cex=size,alpha=0.5)
  #points(jitter(Data$A[which(pos==1)]),jitter(Data$Y[which(pos==1)]),col="red",pch=19,cex=size,alpha=0.5)
  if (saving){
    dev.off()
  }
}


#' @export
predict_indirect<-function(fitted,dta,classification=TRUE){
  if ('cv.glmnet'%in%class(fitted)|'glmnet'%in%class(fitted)) dta=design.lasso(list(A=dta[,1],X=dta[,-1]))

  if (classification){
    if ('ksvm'%in%class(fitted)){
      pred=predict(fitted,dta,type='probabilities')[,2]
    } else if ('randomForest'%in%class(fitted)){
      colnames(dta)=NULL
      pred=predict(fitted,as.matrix(dta),type='prob')[,2]
    } else if ('svm'%in%class(fitted)){
      pred=predict(fitted,as.matrix(dta),probability=TRUE)
      pred=attr(pred, "probabilities")[,1]
    } else if ('cv.glmnet'%in%class(fitted) ){
      pred=predict(fitted,dta,type='response')
    } else{
      pred=predict(fitted,dta,type='prob')
    }
  } else{
    if ('glm'%in%class(fitted)){
      colnames(dta)=NULL
      dta=as.data.frame(dta)
      pred=predict(fitted,dta)
    } else if (('glmnet'%in%class(fitted))){
      pred=predict(fitted,dta)
    } else if ('randomForest'%in%class(fitted)){
      colnames(dta)=NULL
      pred=predict(fitted,as.matrix(dta))
    } else {
      pred=predict(fitted,as.matrix(dta))
    }
  }
  return(pred)
}

#' @export
huntPoint<-function(pred,constant,two.sided,cand,lb,ub){
  NA->pred_L->pred_R
  allup=(sum(pred<=constant)==0)
  alldown=(sum(pred>=constant)==0)
  if (allup|alldown){
    if (allup){
      return(c((lb+ub)/2,lb,ub,2))
    } else {
      return(c(NA,ub,lb,-2))
    }
  }
  diff=abs(pred-constant)
  if (!two.sided){
    tmp=median(cand[which.min(diff)],na.rm = TRUE)
    return(c(tmp,NA,NA,1))
  } else if (two.sided){
    mid=which.max(pred)
    if(mid==1){
      tmp=median(cand[which.min(diff)],na.rm = TRUE)
      pred_L=lb;pred_R=tmp;
    } else if(mid==length(cand)){
      tmp=median(cand[which.min(diff)],na.rm = TRUE)
      pred_L=tmp;pred_R=ub;
    } else {
      predL=(pred[1:mid])
      AL=cand[1:mid]
      predR=(pred[(mid+1):length(pred)])
      AR=cand[(mid+1):length(pred)]
      indexL=which.min(abs(predL-constant))
      indexR=which.min(abs(predR-constant))
      if (sum(indexL)>0){
        tmpL=median(AL[indexL])
      } else{
        tmpL=lb
      }
      if (sum(indexR)>0){
        tmpR=median(AR[indexR])
      } else{
        tmpR=ub
      }
      pred_L=tmpL
      pred_R=tmpR
      pred=mean(c(pred_L,pred_R))
    }
    return(c((pred_L+pred_R)/2,pred_L,pred_R,1))
  }
}

#' @export
searchInterval=function(fitted,test,two.sided=FALSE,lower=TRUE,classification=TRUE, num.plot=5,
                        dt=0.00001,left=NA,right=NA,aL=NULL,aU=NULL){
  if (is.null(aL)) aL=min(test$A)-0.5*sd(test$A)
  if (is.null(aU)) aU=max(test$A)+0.5*sd(test$A)
  if (classification){
    Constant=test$alpha[,1]
  } else {
    Constant=rep(test$S,length(test$Y))
  }
  lb=aL;ub=aU;
  cand=seq(lb,ub,0.1)
  pred_L<-rep(NA,length(test$A))->pred_R->pred
  newdata=cbind(test$A,test$X)
  Rsquare=NA
  correlation=NA

  Dta=cbind(rep(cand,length(test$A)),apply(newdata[,-1],2,function(x) rep(x,each=length(cand))))
  Pred<-predict_indirect(fitted,Dta,classification = classification)->Pred0

  res=matrix(NA,ncol=4,nrow=length(test$A))
  for (i in 1:length(test$A)){
    res[i,]=huntPoint(pred=Pred0[(i-1)*length(cand)+1:length(cand)],constant=Constant[i],two.sided=two.sided,cand=cand,lb=lb,ub=ub)
    if (i<num.plot){
      plot(cand,Pred0[(i-1)*length(cand)+1:length(cand)],'l')
      abline(h=Constant[i])
      abline(v=res[i,1],col='red')
      abline(v=res[i,2],col='red')
      abline(v=res[i,3],col='blue')
    }
  }

  pred=res[,1];pred_L=res[,2];pred_R=res[,3];found=res[,4]

  if (!two.sided){
    pred[is.na(pred)]=median(pred,na.rm = TRUE)
    if (lower){
      region=test$A>pred
    } else {
      region=test$A<pred
    }
    if (!is.null(test$opt)){
      Rsquare=1-mean((pred-test$opt)^2)/var(test$opt)
      correlation=cor(pred,test$opt)
    }
  } else{
    pred_L[is.na(pred_L)]=median(pred_L,na.rm = TRUE)
    pred_R[is.na(pred_R)]=median(pred_R,na.rm = TRUE)
    region=(test$A<pred_R)&(test$A>pred_L)
    if ((!is.null(test$opt_R))&(!is.null(test$opt_L))){
      Rsquare=1-mean((pred_L-test$opt_L)^2)/var(test$opt_L)
      Rsquare=Rsquare+1-mean((pred_R-test$opt_R)^2)/var(test$opt_R)
      Rsquare=Rsquare/2
      correlation=0.5*cor(pred_L,test$opt_L)+0.5*cor(pred_R,test$opt_R)
    }
  }
  cost=(sum(test$alpha[,1]*(region)*pmax((test$S-test$Y),0))+sum(test$alpha[,2]*(!region)*pmax((test$Y-test$S),0)))/length(test$A)
  misclass=(sum(test$alpha[,1]*(region)*(test$Y<test$S))+sum(test$alpha[,2]*(!region)*(test$Y>test$S)))/length(test$A)

  percentage=sum(region*(test$Y>test$S))/sum(region)

  output=list(pred=pred,pred_L=pred_L,pred_R=pred_R,found=found,region=region,misclass=misclass,cost=cost,correlation=correlation,Rsquare=Rsquare,percentage=percentage)
  return(output)
}
