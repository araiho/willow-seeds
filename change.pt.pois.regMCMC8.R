change.pt.pois.regMCMC8 <- function(n.mcmc,X,y.seed,opt.ht.strt,hshape1,hshape2,opt.ht.tune,beta.tune,beta.mns.strt,beta.mns.tune,beta.var.strt,beta.var.tune){
	
	###
	### Subroutines and Libraries 
	###
	
	library(mvtnorm)
	
	###
	### Storage
	###

	opt.ht.save=rep(0,n.mcmc)
	beta.mns.save=matrix(0,n.mcmc,2)
	beta.var.save=rep(0,n.mcmc)
	
	###
	### Initial Conditions
	###
	
	opt.ht = opt.ht.strt	
	W = matrix(c(.60,.30,.30),1,3)
    x = cbind(1,log(as.vector(by(X[,3]>opt.ht,X[,1],sum))))
   
    beta = matrix(1,length(y.seed[,1]),2)
    beta.var=beta.var.strt
    beta.var.mns=0
    beta.var.var=1
    
    beta.mns=beta.mns.strt
    beta.mns.mns=c(0,0)
    beta.mns.var=1
	
	accept.opt.ht = 0
	accept.beta = 0
	accept.beta.mns = 0
	accept.beta.var = 0
	
	###
	### Begin MCMC
	###
		
	for(k in 1:n.mcmc){
    if(k%%1000==0) cat(k," ");flush.console()
    
    ###    
    ### Sample Opt.ht
    ###
    
    opt.ht.star = rnorm(1,opt.ht,opt.ht.tune)
    
    if(opt.ht.star <= max(X[,3]) & opt.ht.star > 0){
       
    x.star = cbind(1,log(as.vector(by(X[,3] > opt.ht.star,X[,1],sum))))
    
 	  mh1.ht = sum(dpois(y.seed[x.star[,2] > 0,],exp(beta[x.star[,2]>0,1]+x.star[x.star[,2]>0,2]*beta[x.star[,2]>0,2])%*%W,log=TRUE)) + dgamma(opt.ht.star,hshape1,hshape2,log=TRUE)
 	  
	  mh2.ht = sum(dpois(y.seed[x[,2] > 0,],exp(beta[x[,2]>0,1]+x[x[,2]>0,2]*beta[x[,2]>0,2])%*%W,log=TRUE)) + dgamma(opt.ht,hshape1,hshape2,log=TRUE)
	  mh.ht = exp(mh1.ht - mh2.ht)
    
    if(mh.ht > runif(1)){
      opt.ht=opt.ht.star
      accept.opt.ht=accept.opt.ht+1
      
    }
    }
   
    ###
    ### Sample Betas    
    ###
    
	  beta.star = cbind(rnorm(21,beta[,1],beta.tune),rnorm(21,beta[,2],beta.tune))
	  
	  if(min(exp(x[x[,2]>0,1]*beta.star[x[,2]>0,1]+x[x[,2]>0,2]*beta.star[x[,2]>0,2])) > 0){

	  mh1.beta = dpois(y.seed[x[,2]>0,],exp(beta.star[x[,2]>0,1]+x[x[,2]>0,2]*beta.star[x[,2]>0,2])%*%W,log=TRUE)+dmvnorm(beta.star[x[,2]>0,],beta.mns,rep(beta.var,2)*diag(2),log=TRUE)
	  mh2.beta = dpois(y.seed[x[,2]>0,],exp(beta[x[,2]>0,1]+x[x[,2]>0,2]*beta[x[,2]>0,2])%*%W,log=TRUE)+dmvnorm(beta[x[,2]>0,],beta.mns,rep(beta.var,2)*diag(2),log=TRUE)
	  mh.beta = exp(mh1.beta-mh2.beta)
	  
	  mh.beta.index = which(mh.beta > runif(length(mh.beta[,1])*length(mh.beta[,1])),arr.ind=TRUE)
	  
	 # mh.index = which(mh.ht > runif(length(mh.ht[,1])*length(mh.ht[1,])),arr.ind=TRUE)
	  
	  beta[mh.beta.index,] = beta.star[mh.beta.index,]
	 
    }
    ###
    ### Sample Beta Means 
    ###
    
	  beta.mns.star = rnorm(2,beta.mns,beta.mns.tune)

	  mh1.beta.mns = sum(dmvnorm(beta[x[,2]>0,],beta.mns.star,rep(beta.var,2)*diag(2),log=TRUE))+dmvnorm(beta.mns.star,beta.mns.mns,rep(beta.mns.var,2)*diag(2),log=TRUE)
	  mh2.beta.mns = sum(dmvnorm(beta[x[,2]>0,],beta.mns,rep(beta.var,2)*diag(2),log=TRUE))+dmvnorm(beta.mns,beta.mns.mns,rep(beta.mns.var,2)*diag(2),log=TRUE)
	  mh.beta.mns = exp(mh1.beta.mns-mh2.beta.mns)

   if(mh.beta.mns > runif(1)){
      beta.mns = beta.mns.star
      accept.beta.mns = accept.beta.mns+1
      
    }
    
    ###
    ### Sample Beta Var 
    ###
    
	  beta.var.star = rnorm(1,beta.var,beta.var.tune)
	  
	  if(beta.var.star > 0){

	  mh1.beta.var = sum(dmvnorm(beta[x[,2]>0,],beta.mns,rep(beta.var.star,2)*diag(2),log=TRUE))+dunif(beta.var.star,beta.var.mns,beta.var.var,log=TRUE)
	  mh2.beta.var = sum(dmvnorm(beta[x[,2]>0,],beta.mns,rep(beta.var,2)*diag(2),log=TRUE))+dunif(beta.var,beta.var.mns,beta.var.var,log=TRUE)
	  mh.beta.var = exp(mh1.beta.var-mh2.beta.var)

   if(mh.beta.var > runif(1)){
      beta.var = beta.var.star
      accept.beta.var = accept.beta.var+1
      
    }
    }

    ###
    ### Save Values
    ###
        
	opt.ht.save[k] = opt.ht
	beta.mns.save[k,] = beta.mns  
	beta.var.save[k] = beta.var  
        
  }
    
    ###
    ### Create Plots
    ###

	n.burn=round(.2*n.mcmc)
	layout(matrix(1:4,4,1))
	plot(opt.ht.save[n.burn:n.mcmc],type="l",main="Threshold Height")
	plot(beta.mns.save[n.burn:n.mcmc,1],type="l",main=expression(paste(beta[0])))
	plot(beta.mns.save[n.burn:n.mcmc,2],type="l",main=expression(paste(beta[1]))) 
	plot(beta.var.save[n.burn:n.mcmc],type="l",main="Beta Var")
	

    
    list(opt.ht.save=opt.ht.save,accept.beta.mns=accept.beta.mns,accept.beta.var=accept.beta.var,accept.opt.ht=accept.opt.ht)
}