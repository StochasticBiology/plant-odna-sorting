### computes fits of heteroplasmy measurement sets to kimura distribution, finds maximum likelihood parameters, hypothesis tests for bottleneck size, also computes wonnapinij-style uncertainty estimates and compiles results

# set to T if you don't have kimura package
install.kimura = F

# install kimura package if required
if(install.kimura) {
  install.packages("devtools")
  library("devtools")
  devtools::install_github("lbozhilova/kimura")
}

library(kimura)

# a transformation function to cast any real number onto the interval ]0,1[
transfun = function(x) {
  return(1/(1+exp(-x)))
}

# calculate negative log likelihood for a given set of heteroplasmy measurements h and parameters theta = {logit(p), logit(b)} (we write h0 for p)
# we can do this enforcing a particular h0 value (passed as an argument) or treating h0 as a fit parameter (default)
# the logit transform is used to ensure h0 and b remain in the [0,1] interval regardless of what real-valued argument the numerical optimiser attempts

kimura_neg_loglik = function(theta, h, h0=F) {
  # get kimura parameters from argument
  b = transfun(theta[1])
  # if we haven't provided a specific h0, retrieve it as a parameter (only difference is the transformation step)
  if(h0 == F) {  h0 = transfun(theta[2]) }
  return(-sum(log(dkimura(h, h0, b))))
}

# joint negative log likelihood function for several families' heteroplasmy measurements
# theta = { b, h0.1, h0.2, ... } [use h values if use.h0s=F, otherwise initial heteroplasmies are enforced via h0s]
joint_neg_log_lik = function(theta, hlist, use.h0s=F, h0s=-1) {
  llik = 0
  for(i in 1:length(hlist)) {
    if(use.h0s == F) {
      llik = llik + kimura_neg_loglik(c(theta[1], theta[i+1]), hlist[[i]], h0=F)
    } else {
      llik = llik + kimura_neg_loglik(theta[1], hlist[[i]], h0=h0s[i])
    }
  }
  return(llik)
} 


# compute maximum likelihood parameters and confidence intervals for a given heteroplasmy set
# we can do this while imposing a specific h0 as an argument or allowing a search over h0 values
maxlik = function(h, conf.level = 0.95, h0 = F) {
  # if we have enforced a particular h0
  if(h0 != F) {
    # find best transformed parameter b
    best = optim(c(0.5), kimura_neg_loglik, h=h, h0=h0, method="Brent", lower=-30, upper=30, hessian=T)
    # add back-transformed parameters to return structure
    # h0 is fixed here; n is straightforward function of b
    best$h0.hat = h0
    best$b.hat = transfun(best$par[1])
    best$n.hat = 1/(1-best$b.hat)

    # the hessian from the optimisation process is the fisher information matrix.
    # its inverse (via "solve") can be used to construct confidence intervals on our parameter(s)
    fisher.matrix = best$hessian
    crit = qnorm((1 + conf.level)/2)
    inv.fisher.matrix <- solve(fisher.matrix)
    b.ci = best$par[1] + c(-1, 1) * crit * sqrt(inv.fisher.matrix[1, 1])
    best$h0.ci = c(h0, h0)
    best$b.ci = transfun(b.ci)

    best$n.hat = 1/(1-best$b.hat)
    best$n.ci = 1/(1-best$b.ci)
  } else {
    # find best transformed parameters h0 and b
    best = optim(c(0.5, 0.5), kimura_neg_loglik, h=h, h0=F, hessian=T)
    # add back-transformed parameters to return structure
    best$b.hat = transfun(best$par[1])
    best$h0.hat = transfun(best$par[2])
    best$n.hat = 1/(1-best$b.hat)

    # the hessian from the optimisation process is the fisher information matrix.
    # its inverse (via "solve") can be used to construct confidence intervals on our parameter(s)
    fisher.matrix = best$hessian
    crit = qnorm((1 + conf.level)/2)
    inv.fisher.matrix <- solve(fisher.matrix)
    b.ci = best$par[1] + c(-1, 1) * crit * sqrt(inv.fisher.matrix[1, 1])
    h0.ci = best$par[2] + c(-1, 1) * crit * sqrt(inv.fisher.matrix[2, 2])
    best$b.ci = transfun(b.ci)
    best$h0.ci = transfun(h0.ci)

    best$n.hat = 1/(1-best$b.hat)
    best$n.ci = 1/(1-best$b.ci)
  }
  
  return(best)
}

# compute bootstrap estimates for parameters and confidence intervals for a given heteroplasmy set
# we can do this while imposing a specific h0 as an argument or allowing a search over h0 values
maxlikboot = function(h, nboot=200, conf.level = 0.95, h0=F) {
  # if we have enforced a particular h0
  if(h0 != F) {
    boot.b = c()
    # loop over bootstrap resamples
    for(boot in 1:nboot) {
      # construct bootstrap sample
      hboot = sample(h, replace=T)
      # find best transformed parameter b
      boot.best = optim(c(0.5), kimura_neg_loglik, h=hboot, h0=h0, method="Brent", lower=-30, upper=30)
      # record back-transformed parameter for this resample
      boot.b = c(boot.b, transfun(boot.best$par[1]))
    }
    # do the optimisation for the non-resampled set
    best = optim(c(0.5), kimura_neg_loglik, h=h, h0=h0)
    best$b.hat = transfun(best$par[1])
    best$n.hat = 1/(1-best$b.hat)
    best$h0.hat = h0

    # get stats and confidence intervals from the bootstrap distribution
    best$b.bhat = mean(boot.b)
    best$h0.ci = c(h0,h0)
    conf.1 = (1-conf.level)/2
    conf.2 = 1-conf.1
    best$b.ci = c(quantile(boot.b, conf.1), quantile(boot.b, conf.2))
  
    best$n.bhat = 1/(1-best$b.bhat)
    best$n.ci = 1/(1-best$b.ci)
  } else {
    boot.h0 = boot.b = c()
    # loop over bootstrap resamples
    for(boot in 1:nboot) {
      hboot = sample(h, replace=T)
      # find best transformed parameters h0 and b
      boot.best = optim(c(0.5, 0.5), kimura_neg_loglik, h=hboot, h0=F)
      # record back-transformed parameters for this resample
      boot.h0 = c(boot.h0, transfun(boot.best$par[2]))
      boot.b = c(boot.b, transfun(boot.best$par[1]))
    }
    # do the optimisation for the non-resampled set
    best = optim(c(0.5, 0.5), kimura_neg_loglik, h=h)
    best$h0.hat = transfun(best$par[2])
    best$b.hat = transfun(best$par[1])
    best$n.hat = 1/(1-best$b.hat)
  
    # get stats and confidence intervals from the bootstrap distribution
    best$h0.bhat = mean(boot.h0)
    best$b.bhat = mean(boot.b)
    best$h0.ci = c(quantile(boot.h0, 1-conf.level), quantile(boot.h0, conf.level))
    best$b.ci = c(quantile(boot.b, 1-conf.level), quantile(boot.b, conf.level))
  
    best$n.bhat = 1/(1-best$b.bhat)
    best$n.ci = 1/(1-best$b.ci)
  }
  return(best)
}

# compute standard error on the variance using the formula from Wonnapinij
# can be done using a population or a sample picture, and enforcing a particular h0 or using the sample mean
se.var = function(h, h0=F, usepopn=F) {
  # if we're not enforcing a particular h0
  if(h0 == F) {
    hbar = mean(h)
  } else {
    hbar = h0
  }
  n = length(h)
  # construct population or sample variance
  if(usepopn == F) { s2 = 1/(n-1) * sum((h-hbar)**2) } else { s2 = 1/n * sum((h-hbar)**2) }
  # moments and intermediate values
  mu2 = sum((h-hbar)**2)/n
  mu4 = sum((h-hbar)**4)/n
  D4 = (n-1)/n**3 * ((n**2-3*n+3)*mu4 + 3*(2*n-3)*mu2**2)
  return( sqrt(1/n*(D4 - (n-3)/(n-1)*s2**2)) )
}

# calculate various statistics for a heteroplasmy set h
# can enforce an initial h0 or leave as a free parameter
# can use population or sample statistics
hstats = function(h, h0=F, usepopn=F) {
  # initialise list of results
  statres = list()
  # basic stats
  n = length(h)
  # if we're not enforcing a particular h0, compute summary stats
  if(h0 == F) {
    hbar = mean(h)
    s2 = var(h)
    if(usepopn == T) { s2 = (n-1)*s2/n }
    sehbar = sqrt(s2/n)
  } else {
    # otherwise, compute stats assuming the given h0
    hbar = h0
    if(usepopn == F) { s2 = 1/(n-1) * sum((h-hbar)**2) } else { s2 = 1/n * sum((h-hbar)**2) }
    sehbar = 0
  }
  # standard error on the variance and confidence intervals (Wonnapinij) - - cis rely on normality assumption of course 
  sev = se.var(h, usepopn)
  statres$moment.h0.hat = hbar
  statres$moment.h0.ci = c( hbar-1.96*sehbar, hbar+1.96*sehbar )

# deal with the case where we have zero variance (identical samples). In practise this only happens in homoplasmic cases, where b and n are undefined (no information about bottleneck size) 
  if(s2 == 0) {
    statres$moment.n.hat = NA
    statres$moment.n.ci = c(NA, NA)
  } else {
    statres$moment.n.hat = hbar*(1-hbar)/s2
    statres$moment.n.ci = c( (hbar*(1-hbar)/(s2+1.96*sev)), (hbar*(1-hbar)/(s2-1.96*sev)))
  }
  
  # for Kimura fit, decide whether to use bootstrapping or Fisher information for uncertainty
  # looks like homoplasmic samples reward bootstrapping (numerical convergence issues), and others prefer Fisher
  if(all(h == 0 | h == 1)) {
    print("Homoplasmic samples -- bootstrapping")
    fit = maxlikboot(h, h0=h0)
  } else {
    print("Heteroplasmic samples -- using Fisher information")
    fit = maxlik(h, h0=h0)
  }
  # estimates and confidence intervals from Kimura fit
  statres$fit.h0.hat = fit$h0.hat
  statres$fit.h0.ci = fit$h0.ci
  statres$fit.n.hat = fit$n.hat
  statres$fit.n.ci = fit$n.ci

  if(is.na(statres$moment.n.hat)) {
    statres$preferred = "fit"
  } else {
    # which approach is "better"? if confidence intervals on n go below 1, discard that approach; otherwise if one set of confidence intervals fits inside the other, favour that one
    if(statres$moment.n.ci[1] < 0 | statres$moment.n.ci[2] < 0 | (statres$fit.n.ci[1] > statres$moment.n.ci[1] & statres$fit.n.ci[2] < statres$moment.n.ci[2])) {
      statres$preferred = "fit"
    } else if(statres$fit.n.ci[1] < statres$moment.n.ci[1] & statres$fit.n.ci[2] > statres$moment.n.ci[2]) {
      statres$preferred = "moment"
    } else {
      statres$preferred = "neither"
    }
  }
  
  return(statres)
}

# perform likelihood ratio test exploring difference in bottleneck size between two heteroplasmy samples
kimura_lrt = function(h1, h2, use.h0s=F, h1.h0set=0, h2.h0set=0) {
  
  comp.df = data.frame()

  if(use.h0s) {
    # single-param optimisation for PT, MT, and combined h sets
    h1.best = optim(-3, joint_neg_log_lik, hlist=h1, use.h0s=T, h0s=h1.h0set, method="Brent", lower=-30, upper=30, hessian=T)
    h2.best = optim(-3, joint_neg_log_lik, hlist=h2, use.h0s=T, h0s=h2.h0set, method="Brent", lower=-30, upper=30, hessian=T)
    both.best = optim(-3, joint_neg_log_lik, hlist=c(h1,h2), use.h0s=T, h0s=c(h1.h0set,h2.h0set), method="Brent", lower=-30, upper=30, hessian=T)
  } else {
    # multi-param optimisation for PT, MT, and combined h sets
    h1.best = optim(c(-3, rep(0.5, length(h1))), joint_neg_log_lik, hlist=h1, use.h0s=F, hessian=T)
    h2.best = optim(c(-3, rep(0.5, length(h2))), joint_neg_log_lik, hlist=h2, use.h0s=F, hessian=T)
    both.best = optim(c(-3, rep(0.5, length(c(h1,h2)))), joint_neg_log_lik, hlist=c(h1, h2), use.h0s=F, hessian=T)
  }

  # use Fisher information to get confidence intervals on bottleneck size estimates
  conf.level = 0.95
  crit = qnorm((1 + conf.level)/2)
  h1.ci = h1.best$par[1] + c(-1, 1) * crit * sqrt(solve(h1.best$hessian)[1, 1])
  h2.ci = h2.best$par[1] + c(-1, 1) * crit * sqrt(solve(h2.best$hessian)[1, 1])
  both.ci = both.best$par[1] + c(-1, 1) * crit * sqrt(solve(both.best$hessian)[1, 1])
    
  h1.n.hat = 1/(1-transfun(h1.best$par[1]))
  h1.n.ci = 1/(1-transfun(h1.ci))
  h2.n.hat = 1/(1-transfun(h2.best$par[1]))
  h2.n.ci = 1/(1-transfun(h2.ci))
  both.n.hat = 1/(1-transfun(both.best$par[1]))
  both.n.ci = 1/(1-transfun(both.ci))

  # get log-likelihoods for the separate (MT =/= PT) and both (MT = PT) models (returned values are negative log likelihoods, so take negatives here) 
  sep.llik = -h2.best$value-h1.best$value
  both.llik = -both.best$value

  # likelihood ratio test and p-value from chi-squared distribution (one dof -- one parameter difference)
  lrt = -2*(both.llik - sep.llik)
  pval = pchisq(lrt, 1, lower.tail=F)
  comp.df = rbind(comp.df, data.frame(h1.n.hat, h1.n.ci[1], h1.n.ci[2], h2.n.hat, h2.n.ci[1], h2.n.ci[2], both.n.hat, both.n.ci[1], both.n.ci[2], pval))
  return(comp.df)
}

