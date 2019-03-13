main <-
function(y,x,S,G,is.graph="IVGC",nlambda=50,nalpha=50,nfold=5,epsilon=1e-6){

  p=ncol(x)
  N=nrow(x)
  q=ncol(S)
  limit=100
  df=1
  J=p

  ### estimate Gamma
  Gamma_hat=matrix(0,q+1,p)
  for(i in 1:p)
  {
    cvfit=cv.glmnet(S,x[,i],family="gaussian")
    coe=coef(cvfit,s="lambda.min")
    Gamma_hat[,i]=as.matrix(coe)
  }

  ### obtain the fitted X
  S1=cbind(1,S)
  X0_hat=S1%*%Gamma_hat

  ### standardize X0_hat and center Y0
  X0_hat.sd=apply(X0_hat,2,sd)
  X0_hat.sd=X0_hat.sd*sqrt((N-1)/N)
  X0_hat.final=matrix(0,N,p)
  for(i in 1:J) X0_hat.final[,i]=(X0_hat[,i]-mean(X0_hat[,i]))/X0_hat.sd[i]
  Y0.final=y-mean(y)
  if(is.graph != "IV")
  {
    ### set up tuning parameters
    alpha0=rep(0,nalpha)
    lambda0=matrix(0,nalpha,nlambda)
    
    for(i in 1:nalpha) alpha0[i]=1/nalpha*i
    
    corr_xy=rep(0,p)
    for(i in 1:p) corr_xy[i]=sum(X0_hat.final[,i]*Y0.final)
    
    for(i in 1:nalpha)
    {
      lambdaMax=max(abs(corr_xy))/alpha0[i]/N
      lambdaMin=epsilon*lambdaMax
      d=(log(lambdaMax)-log(lambdaMin))/nlambda
      for(j in 1:nlambda) lambda0[i,j]=lambdaMin*exp((j+1)*d)
    }
    
    ### k(=5) fold cross-validation step to select the tuning parameter
    num=N/nfold
    index=sample(1:N)
    RSS.gc=matrix(0,nfold,nlambda*nalpha)
    
    for(ff in 1:nfold)
    {
      ind2=index[((ff-1)*num+1):(ff*num)]
      ind1=setdiff(1:N,ind2)
      
      ### training data
      n=N-num
      X=X0_hat[ind1,]
      Y=y[ind1]
      
      X.sd=apply(X,2,sd)
      X.sd=X.sd*sqrt((n-1)/n)
      X.final=matrix(0,n,p)
      for(i in 1:J) X.final[,i]=(X[,i]-mean(X[,i]))/X.sd[i]
      Y.final=Y-mean(Y)
      
      #### graph constrained estimation
      vec=c(n,J,df,nlambda,nalpha,limit)
      
      gammaSolution=matrix(0,nalpha*nlambda,df*J)
      fit<-.C("GraphConstrainedEstimation",as.double(gammaSolution),as.integer(vec),as.double(Y.final),as.double(X.final),as.double(G),as.double(alpha0),as.double(lambda0))
      gammaSolution=matrix(fit[[1]],nalpha*nlambda,df*J)
      
      A=diag((1/X.sd))
      gammaSolution1=gammaSolution%*%A
      
      mm=apply(X,2,mean)
      vv=mm/X.sd
      cons=mean(Y)-gammaSolution%*%vv
      
      gammaSolution.gc=cbind(cons,gammaSolution1)
      
      ### testing data
      n=num
      X=X0_hat[ind2,]
      X1=cbind(1,X)
      Y=y[ind2]
      
      ### graph constrained
      Y_hat=X1%*%t(gammaSolution.gc)
      for(i in 1:(nalpha*nlambda)) RSS.gc[ff,i]=sum((Y-Y_hat[,i])^2)/n
    }
    
    ###### use the complete data
    ### graph constrained estimation
    err=apply(RSS.gc,2,mean)
    ind=which(err==min(err))[1]
    
    ind.alpha=as.integer(ind/nlambda)+1
    ind.lambda=ind-(ind.alpha-1)*nlambda
    
    vec=c(N,J,df,nlambda,1,limit)
    
    gammaSolution=matrix(0,nlambda,df*J)
    fit<-.C("GraphConstrainedEstimation",as.double(gammaSolution),as.integer(vec),as.double(Y0.final),as.double(X0_hat.final),as.double(G),as.double(alpha0[ind.alpha]),as.double(lambda0[ind.alpha,]))
    gammaSolution=matrix(fit[[1]],nlambda,df*J)
    
    A=diag((1/X0_hat.sd))
    gammaSolution1=gammaSolution%*%A
    
    mm=apply(X0_hat,2,mean)
    vv=mm/X0_hat.sd
    cons=mean(y)-gammaSolution%*%vv
    gammaSolution.gc=cbind(cons,gammaSolution1)
    lambda.gc=gammaSolution.gc[ind.lambda,]
  }

  #### estimation without incorporating the graph structure
  if(is.graph != "IVGC"){
    cvfit=cv.glmnet(X0_hat,y,family="gaussian")
    coe=coef(cvfit,s="lambda.min")
    lambda.i=as.vector(coe)
  }
  if(is.graph == "IVGC")  fitivgc = list(beta_ivgc=lambda.gc, Gamma_hat=Gamma_hat)
  else if(is.graph == "IV")  fitivgc = list(beta_iv=lambda.i, Gamma_hat=Gamma_hat)
  else  fitivgc = list(beta_ivgc=lambda.gc,beta_iv=lambda.i, Gamma_hat=Gamma_hat)
  
  return(fitivgc)
}
