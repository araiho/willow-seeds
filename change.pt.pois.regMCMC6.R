change.pt.pois.regMCMC6 <- function(n.mcmc,X,y.seed,opt.ht.strt,hshape1,hshape2,opt.ht.tune,beta,beta.mns,beta.tune,beta.var,N,N.tune,N.shape,N.rate,x.strt){
	
	###
	### Subroutines and Libraries 
	###
	
	library(mvtnorm)
	library(mgcv)
	library(LaplacesDemon)
	
	###
	### Storage
	###

	opt.ht.save=rep(0,n.mcmc)
	N.save=rep(0,n.mcmc)
	beta.save=matrix(0,n.mcmc,2)
  	x.save=rep(0,n.mcmc)
  
	mse.y=rep(0,n.mcmc)
	mse.ypred=rep(0,n.mcmc)
	msediffsave=rep(0,n.mcmc)
	ypredsave=matrix(0,length(y.seed),n.mcmc)
	msesave=rep(0,n.mcmc)
	
	###
	### Initial Conditions
	###
	
	opt.ht = opt.ht.strt	
	W=matrix(c(.6,.3,.3),1,3)
  	x = cbind(1,rep(x.strt,length(y.seed[,1])))
  	sum.vec=as.vector(by(X[,3]>opt.ht,X[,1],sum))
  	p=as.vector(rdirichlet(1,sum.vec+1))
   
  	beta = beta
	
	accept.opt.ht = 0
	accept.beta = 0
	accept.N = 0
  	accept.x = 0
  	accept.p = 0
	
	###
	### Begin MCMC
	###
		
	for(k in 1:n.mcmc){
    if(k%%1000==0) cat(k," ");flush.console()
    
    ###    
    ### Sample Opt.ht
    ###
    
    opt.ht.star = rnorm(1,opt.ht,opt.ht.tune)
    
    if(opt.ht.star < max(X[,3]) & opt.ht.star > 1){
      
    sum.vec.star=as.vector(by(X[,3]>opt.ht.star,X[,1],sum))
    
    mh1.ht = sum(ddirichlet(p,sum.vec.star+1,log=TRUE) + dlnorm(opt.ht.star,hshape1,hshape2,log=TRUE)) 
    mh2.ht = sum(ddirichlet(p,sum.vec+1,log=TRUE) + dlnorm(opt.ht,hshape1,hshape2,log=TRUE))
    mh.ht = exp(mh1.ht - mh2.ht)
    
    if(mh.ht > runif(1)){
      opt.ht=opt.ht.star
      sum.vec=sum.vec.star
      accept.opt.ht=accept.opt.ht+1
      
    }
    }
    
    ###    
    ### Sample p
    ###
    
    p.star=as.vector(rdirichlet(1,sum.vec+1))
    
    mh1.p=sum(dcat(x[,2],p.star,log=TRUE)+ddirichlet(p.star,sum.vec+1,log=TRUE))
    mh2.p=sum(dcat(x[,2],p,log=TRUE)+ddirichlet(p,sum.vec+1,log=TRUE))
    
    mh.p=exp(mh1.p-mh2.p)
    
    if(mh.p > runif(1)){
      p = p.star
      accept.p = accept.p+1
      
    }
    
    ###    
    ### Sample x
    ###
    
    x.star = cbind(1,rcat(length(y.seed[,1]),rep(1,length(y.seed[,1]))))
    
    lambda.star = notExp(x.star%*%beta)%*%W
    
 	mh1.x = dnbinom(y.seed,mu=lambda.star,size=N,log=TRUE) + dcat(x.star[,2],p,log=TRUE)
	mh2.x = dnbinom(y.seed,mu=notExp(x%*%beta)%*%W,size=N,log=TRUE) + dcat(x[,2],p,log=TRUE)
	mh.x = exp(mh1.x - mh2.x)
	
	mh.index = which(mh.x > runif(length(mh.x[,1])*length(mh.x[1,])),arr.ind=TRUE)
    x[unique(mh.index[,1]),] = x.star[unique(mh.index[,1]),]

    ###
    ### Sample Betas    
    ###
    
	beta.star = rnorm(2,beta,beta.tune)

	mh1.beta = sum(dnbinom(y.seed,mu=notExp(x%*%beta.star)%*%W,size=N,log=TRUE))+dmvnorm(beta.star,beta.mns,rep(beta.var,2)*diag(2),log=TRUE)
	mh2.beta = sum(dnbinom(y.seed,mu=notExp(x%*%beta)%*%W,size=N,log=TRUE))+dmvnorm(beta,beta.mns,rep(beta.var,2)*diag(2),log=TRUE) #solve(t(x[x[,2]>0,])%*%x[x[,2]>0,])
	mh.beta = exp(mh1.beta-mh2.beta)
   
    if(mh.beta > runif(1)){
      beta = beta.star
      accept.beta = accept.beta+1
      
    }

    ###    
    ### Sample N
    ###
    
    N.star=rnorm(1,N,N.tune)
    
    if (N.star > 0){
    
    mh1.N = sum(dnbinom(y.seed,mu=notExp(x%*%beta)%*%W,size=N.star,log=TRUE) + dlnorm(N.star,N.shape,N.rate,log=TRUE))
    mh2.N = sum(dnbinom(y.seed,mu=notExp(x%*%beta)%*%W,size=N,log=TRUE) + dlnorm(N,N.shape,N.rate,log=TRUE))
    mh.N = exp(mh1.N - mh2.N)
    
    if(mh.N > runif(1)){
      N=N.star
      accept.N=accept.N+1
      
    }
    }
    
    ###
    ### Obtain Predictions 
    ###
    
    mn.nbinom=notExp(x%*%beta)
    ypred=rnbinom(length(mn.nbinom),mu=mn.nbinom,size=N)
    
    mse.y[k]=mean((rowSums(y.seed)-mn.nbinom)^2)
    mse.ypred[k]=mean((ypred-mn.nbinom)^2)
    
    ###
    ### Save Values
    ###
        
	opt.ht.save[k] = opt.ht
    N.save[k] = N
	beta.save[k,] = beta
    x.save[k]=mean(x[,2])
    
  msediffsave[k]=mse.ypred[k]-mse.y[k]
  #ypredsave[,k]=ypred  
        
  }
  
	n.burn=round(.2*n.mcmc)  
	p.value=sum(mse.ypred[n.burn:n.mcmc] > mse.y[n.burn:n.mcmc]) / (n.mcmc - n.burn)
	
    ###
    ### Create Plots
    ###

	layout(matrix(1:2,2,1))
	plot(opt.ht.save[n.burn:n.mcmc],type="l",main="Threshold Height")
	plot(beta.save[n.burn:n.mcmc,1],type="l",main=expression(paste(beta[0])))
	plot(beta.save[n.burn:n.mcmc,2],type="l",main=expression(paste(beta[1]))) 
	plot(N.save[n.burn:n.mcmc],type="l",main="N")
	plot(x.save[n.burn:n.mcmc],type="l",main="x")
	

    
    list(x.save=x.save,accept.x=accept.x,p.value=p.value,N.save=N.save,opt.ht.save=opt.ht.save,beta.save=beta.save,accept.beta=accept.beta,accept.N=accept.N,accept.opt.ht=accept.opt.ht)
}