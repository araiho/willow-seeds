change.pt.pois.regMCMC4 <- function(n.mcmc,X,y.seed,opt.ht.strt,hshape1,hshape2,opt.ht.tune,beta,beta.mns,beta.tune,beta.var){
	
	###
	### Subroutines and Libraries 
	###
	
	library(mvtnorm)
	
	###
	### Storage
	###

	opt.ht.save=rep(0,n.mcmc)
	beta.save=matrix(0,n.mcmc,2)
	
	###
	### Initial Conditions
	###
	
	opt.ht = opt.ht.strt	
	W = matrix(c(.5,.25,.25),1,3)
    x = cbind(1,log(as.vector(by(X[,3]>opt.ht,X[,1],sum))))
   
    beta = beta
	
	accept.opt.ht = 0
	accept.beta = 0
	
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
        
    x.star = cbind(1,log(as.vector(by(X[,3] > opt.ht.star,X[,1],sum))))
    
 	  mh1.ht = sum(dpois(y.seed[x.star[,2] > 0,],exp(x.star[x.star[,2]>0,]%*%beta)%*%W,log=TRUE)) + dnorm(opt.ht.star,hshape1,hshape2,log=TRUE)
	  mh2.ht = sum(dpois(y.seed[x[,2] > 0,],exp(x[x[,2]>0,]%*%beta)%*%W,log=TRUE)) + dnorm(opt.ht,hshape1,hshape2,log=TRUE)
	  mh.ht = exp(mh1.ht - mh2.ht)
    
    if(mh.ht > runif(1)){
      opt.ht=opt.ht.star
      accept.opt.ht=accept.opt.ht+1
      
    }
    }
   
    ###
    ### Sample Betas    
    ###
    
    if (FALSE){
	  beta.star = rnorm(2,beta,beta.tune)

	  mh1.beta = sum(dpois(y.seed[x[,2]>0,],exp(x[x[,2]>0,]%*%beta.star)%*%W,log=TRUE))+dmvnorm(beta.star,beta.mns,diag(rep(beta.var,2)),log=TRUE)
	  mh2.beta = sum(dpois(y.seed[x[,2]>0,],exp(x[x[,2]>0,]%*%beta)%*%W,log=TRUE))+dmvnorm(beta,beta.mns,diag(rep(beta.var,2)),log=TRUE) #rep(beta.var,2)*diag(2)#solve(t(x[x[,2]>0,])%*%x[x[,2]>0,])
	  mh.beta = exp(mh1.beta-mh2.beta)

   if(mh.beta > runif(1)){
      beta = beta.star
      accept.beta = accept.beta+1
      
    }
    }
		beta=c(-1,1)

    ###
    ### Save Values
    ###
        
	opt.ht.save[k] = opt.ht
	beta.save[k,] = beta    
        
  }
    
    ###
    ### Create Plots
    ###

	n.burn=round(.2*n.mcmc)
	layout(matrix(1:3,3,1))
	plot(opt.ht.save[n.burn:n.mcmc],type="l",main="Threshold Height")
	plot(beta.save[n.burn:n.mcmc,1],type="l",main=expression(paste(beta[0])))
	plot(beta.save[n.burn:n.mcmc,2],type="l",main=expression(paste(beta[1]))) 
	

    
    list(opt.ht.save=opt.ht.save,beta.save=beta.save,accept.beta=accept.beta,accept.opt.ht=accept.opt.ht)
}