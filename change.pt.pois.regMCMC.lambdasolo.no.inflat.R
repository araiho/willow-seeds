change.pt.pois.regMCMC.lambdasolo <- function(n.mcmc,X,y.seed,opt.ht.strt,hshape1,hshape2,opt.ht.tune,lambda.shape,lambda.tune,lambda.rate,lambda.strt){
	
	###
	### Subroutines and Libraries 
	###
	
	library(mvtnorm)
	
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
    x = as.vector(by(X[,3]>opt.ht,X[,1],sum))
	
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
        
    x.star = as.vector(by(X[,3] > opt.ht.star,X[,1],sum))
    
 	  mh1.ht = sum(dpois(y.seed[x.star > 0],lambda,log=TRUE)) + dgamma(opt.ht.star,hshape1,hshape2,log=TRUE)
	  mh2.ht = sum(dpois(y.seed[x > 0],lambda,log=TRUE)) + dgamma(opt.ht,hshape1,hshape2,log=TRUE)
	  mh.ht = exp(mh1.ht - mh2.ht)
    
    if(mh.ht > runif(1)){
      opt.ht=opt.ht.star
      accept.opt.ht=accept.opt.ht+1
      x=x.star
    }
    }
   
    ###
    ### Sample Betas    
    ###
    
    if (TRUE){
	  lambda.star = rnorm(1,lambda,lambda.tune)

	if(lambda.star > 0){

	  mh1.lambda = sum(dpois(y.seed[x > 0],lambda.star,log=TRUE))+dgamma(lambda.star,lambda.shape,lambda.rate,log=TRUE)
	  mh2.lambda = sum(dpois(y.seed[x > 0],lambda,log=TRUE))+dgamma(lambda,lambda.shape,lambda.rate,log=TRUE)
	  mh.lambda = exp(mh1.lambda-mh2.lambda)

   if(mh.lambda > runif(1)){
      lambda = lambda.star
      accept.lambda = accept.lambda+1
      
    }
    }
    }
	

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