constructData<-function(y,x,a,s,propensity=NULL,K=20){
  Data=list()
  Data$X=as.matrix(x)
  Data$Y=drop(y)
  Data$A=matrix(a,ncol=1)
  Data$S=s
  
  if (is.null(propensity)) {
    if (length(unique(a))>20){
      grid=seq(min(a),max(a),length.out = K)
      grid[1]=-Inf;grid[length(grid)]=Inf
      a=cut(a,breaks=grid)
    } else{
      a=as.character(a)
    } 
    dat=as.data.frame(cbind(a,x))
    names(dat)=c('a',paste0('x',1:dim(x)[2]))
    glmfit=nnet::multinom(a~.,data=dat)
    propensity=rep(0,length(a))
    for (i in seq_along(propensity)){
      propensity[i]=1/glmfit$fitted.values[i,a[i]]
    }
  }
  Data$propensity=propensity/mean(propensity)
  return(Data)
}
subsetData<-function(Data,index){
  NewData=list()
  NewData$X=Data$X[index,]
  if (!is.null(Data$A)) NewData$A=matrix(Data$A[index],ncol=1)
  for (i in seq_along(Data)){
    if (!names(Data)[i]%in%c('X','A')) {
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

valueInterval<-function(x,region,type=c('mean','median','var','misclass','mse')){
  
}



predict_indirect=function(fitted,test,two.sided,dt=0.01,left=NA,right=NA,type='PDI',lower=TRUE,responseType='response'){
  lb=aL;ub=aU
  n1=length(test$A)
  cand=seq(lb,ub,0.05)
  pred_L<-rep(NA,length(test$A))->pred_R->pred
  found=rep(1,length(test$A))
  newdata=cbind(test$A,test$X)
  for (i in 1:length(test$A)){
    dta=cbind(cand,matrix(rep(newdata[i,-1],length(cand)),nrow=length(cand),byrow = TRUE))
    if ('cv.glmnet'%in%class(fitted)|'glmnet'%in%class(fitted)) dta=design.lasso(list(A=dta[,1],X=dta[,-1]))
    if ('glm'%in%class(fitted)){
      colnames(dta)[1]='treatment'
      dta=as.data.frame(dta)
      pred=predict(fitted,dta,type=responseType)
    } else if (('glmnet'%in%class(fitted))){
      pred=predict(fitted,dta,type=responseType)
    } else if ('randomForest'%in%class(fitted)){
      colnames(dta)=NULL
      pred=predict(fitted,as.matrix(dta),type=responseType)
    } else {
      pred=predict(fitted,as.matrix(dta))
      if ((responseType =='prob')&('svm'%in%class(fitted))) pred=predict(fitted,as.matrix(dta),probability=TRUE)
    }
    if (sum(pred>test$S)==0|sum(pred<test$S)==0){pred_L[i]=left;pred_R[i]=right;found[i]=0;next}
    if (!two.sided){
      tmp=median(cand[abs(pred-test$S)<dt],na.rm = TRUE)
      pred_L[i]=tmp
    } else if (two.sided){
      mid=which.max(pred)
      if(mid==1){
        tmp=median(cand[abs(pred-test$S)<dt],na.rm = TRUE)
        pred_L[i]=lb;pred_R[i]=mid;
        found[i]=1;
        next
      }
      if(mid==length(cand)){
        tmp=median(cand[abs(pred-test$S)<dt],na.rm = TRUE)
        pred_L[i]=mid;pred_R[i]=ub;
        found[i]=1;
        next
      }
      predL=(pred[1:mid])
      AL=cand[1:mid]
      predR=(pred[(mid+1):length(pred)])
      AR=cand[(mid+1):length(pred)]
      indexL=abs(predL-test$S)<dt
      indexR=abs(predR-test$S)<dt
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
      pred_L[i]=tmpL
      pred_R[i]=tmpR
    }
  }
  if (!two.sided){
    pred_L[is.na(pred_L)]=median(pred_L,na.rm = TRUE)
    region=test$A>pred_L
    Rsquare=1-mean((pred_L-test$opt)^2)/var(test$opt)
    correlation=cor(pred_L,test$opt)
  } else{
    pred_L[is.na(pred_L)]=median(pred_L,na.rm = TRUE)
    pred_R[is.na(pred_R)]=median(pred_R,na.rm = TRUE)
    region=(test$A<pred_R)&(test$A>pred_L)
    Rsquare=1-mean((pred_L-test$opt_L)^2)/var(test$opt_L)
    Rsquare=Rsquare+1-mean((pred_R-test$opt_R)^2)/var(test$opt_R)
    Rsquare=Rsquare/2
    correlation=0.5*cor(pred_L,test$opt_L)+0.5*cor(pred_R,test$opt_R)
  }
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
  percentage=sum(region*(test$Y>test$S))/sum(region)
  
  output=list(pred=pred_L,pred_L=pred_L,pred_R=pred_R,found=found,region=region,misclass=misclass,correlation=correlation,Rsquare=Rsquare,percentage=percentage)
  return(output)
}

