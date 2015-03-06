change.pt.pois.regMCMC.lambdasolo <- function(n.mcmc,X,y.seed,opt.ht.strt,hshape1,hshape2,opt.ht.tune,lambda.shape,lambda.tune,lambda.rate,lambda.strt){
	
	###
	### Subroutines and Libraries 
	###
	
	library(mvtnorm)
	library(LaplacesDemon)
	runifdisc<-function(n, min=0, max=1) sample(min:max, n, replace=T)


	
	###
	### Storage
	###

	opt.ht.save=rep(0,n.mcmc)
	lambda.save=rep(0,n.mcmc)
	
	###
	### Initial Conditions
	###
	
	opt.ht = opt.ht.strt	
	lambda = lambda.strt
	W = matrix(c(.5,.25,.25),1,3)
    h = as.vector(by(X[,3]>opt.ht,X[,1],sum))/sum(as.vector(by(X[,3]>opt.ht,X[,1],sum)))
	
	accept.opt.ht = 0
	accept.lambda = 0
	
	###
	### Begin MCMC
	###
		
	for(k in 1:n.mcmc){
    if(k%%1000==0) cat(k," ");flush.console()
    
    ###    
    ### Sample Opt.ht
    ###
    
    opt.ht.star = rnorm(1,opt.ht,opt.ht.tune)
    
    if(opt.ht.star <= max(X[,3]) & opt.ht.star>0){
    	
      h.star= as.vector(by(X[,3]>opt.ht.star,X[,1],sum))/sum(as.vector(by(X[,3]>opt.ht.star,X[,1],sum)))
    
 	  mh1.ht = sum(dmultinom(x,h.star,size=NULL,log=TRUE)) + dgamma(opt.ht.star,hshape1,hshape2,log=TRUE)
	  mh2.ht = sum(dmultinom(x,h,size=NULL,log=TRUE)) + dgamma(opt.ht,hshape1,hshape2,log=TRUE)
	  mh.ht = exp(mh1.ht - mh2.ht)
    
    if(mh.ht > runif(1)){
      opt.ht=opt.ht.star
      accept.opt.ht=accept.opt.ht+1
      h=h.star
    }
    }

    ###    
    ### Sample X
    ###
    
    x.star = runifdisc(20,0,190)
      
 	mh1.x = sum(dpois(y.seed[x.star > 0],lambda,log=TRUE)) + dmultinom(x.star,h,size=NULL,log=TRUE)
	mh2.x = sum(dpois(y.seed[x > 0],lambda,log=TRUE)) + dmultinom(x,h,size=NULL,log=TRUE)
	mh.x = exp(mh1.x - mh2.x)
	
	if(mh.x > runif(1)){
      x=x.star

    }
    
    #mh.index = which(mh.x > runif(length(mh.x)),arr.ind=TRUE)
    #x[unique(mh.index)] = x.star[unique(mh.index)]

   
    ###
    ### Sample Betas    
    ###
    
	  beta.star = rnorm(2,beta,beta.tune)

	  lambda.star=exp(beta.star%*%cbind(1,log(x))
	  
	  mh1.beta = sum(dpois(y.seed[x[,2]>0,],lambda.star,log=TRUE))+dmvnorm(beta.star,beta.mns,diag(rep(beta.var,2)),log=TRUE)
	  mh2.beta = sum(dpois(y.seed[x[,2]>0,],exp(x[x[,2]>0,]%*%beta)%*%W,log=TRUE))+dmvnorm(beta,beta.mns,diag(rep(beta.var,2)),log=TRUE) #rep(beta.var,2)*diag(2)#solve(t(x[x[,2]>0,])%*%x[x[,2]>0,])
	  mh.beta = exp(mh1.beta-mh2.beta)

   if(mh.beta > runif(1)){
      beta = beta.star
      accept.beta = accept.beta+1
      
    }
	beta[2]=1

    ###
    ### Save Values
    ###
        
	opt.ht.save[k] = opt.ht
	lambda.save[k] = lambda    
        
  }
    
    ###
    ### Create Plots
    ###

	n.burn=round(.2*n.mcmc)
	layout(matrix(1:3,3,1))
	plot(opt.ht.save[n.burn:n.mcmc],type="l",main="Threshold Height")
	plot(lambda.save[n.burn:n.mcmc],type="l",main=expression(paste(lambda)))

	

    
    list(opt.ht.save=opt.ht.save,lambda.save=lambda.save,accept.lambda=accept.lambda,accept.opt.ht=accept.opt.ht)
}