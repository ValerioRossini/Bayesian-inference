df=read.csv("dugongs.txt",sep="\t")

names(df)
x = df$Age
y = df$Length

L = function(alpha, beta, gamma, tau_square, x, y){
  n = length(x)
  constant = 1/(2*pi*tau_square^2)^(n/2)
  s = exp(-0.5*(1/tau_square^2)*(sum((y-alpha+beta*(gamma^x)))^2))
  out = constant*s*(alpha>1 & alpha<Inf)*(beta>1 & beta<Inf)*(gamma>0 & gamma<1)*(tau_square>0 & tau_square<Inf)
  return(out)

}


L2 = function(pars, x, y){
  n = length(x)
  constant = 1/(2*pi*pars[4]^2)^(n/2)
  s = exp(-0.5*(1/pars[4]^2)*(sum((y-pars[1]+pars[2]*(pars[3]^x)))^2))
  out = constant*s*(pars[1]>1 & pars[1]<Inf)*(pars[2]>1 & pars[2]<Inf)*(pars[3]>0 & pars[3]<1)*(pars[4]>0 & pars[4]<Inf)
  return(-out)
  
}

#Likelihood = function(alpha) return(dnorm(x,0,sigma_alpha))

L(2,2,0.5,2,x,y)
#sigma_alpha=10000
sigma_alpha=rexp(1,3)
prior_alpha = function(alpha) return(dnorm(alpha,0,sigma_alpha))

#sigma_beta=10000
sigma_beta=rexp(1,3)
prior_beta = function(beta) return(dnorm(beta,0,sigma_beta))

prior_gamma = function(gamma) return (dunif(gamma))

library(MCMCpack)
a = rexp(1,30)
b = rexp(1,30)
prior_tau_square = function(tau_square,a,b) return(dinvgamma(tau_square, a, b))

joint_prior = function(alpha, beta, gamma, tau_square){
  out = prior_alpha(alpha)*prior_beta(beta)*prior_gamma(gamma)*prior_tau_square(tau_square,a,b)
  return(out)
}

posterior = function (alpha, beta, gamma, tau_square, x, y){
  out = joint_prior(alpha, beta, gamma, tau_square)*L(alpha, beta, gamma, tau_square,x,y)
  return(out)
}


#?optim
optim(c(2,2,0.5,2), L2, x=x, y=y)

joint_prior2 = function(pars){
  out = prior_alpha(pars[1])*prior_beta(pars[2])*prior_gamma(pars[3])*prior_tau_square(pars[4],a,b)
  return(out)
}

posterior2 = function (pars, x, y){
  out = joint_prior2(pars)*L2(pars,x,y)
  return(-out)
}

optim(c(2,2,0.5,4), posterior2, x=x, y=y)






