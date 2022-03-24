# functions for analysis
source("kimura-functions.R")

# read in data from Excel files
source("data-tables.R")

res.text = c()		# this will store outputs as text
use.popn = T		# to match Amanda's calculations

# condensed rounding function for output
r = function(x) {
  return(round(x, digits=2))
}

## table 2 analysis -- between generations in MSH1

# get and store estimates for MT samples
t2.mt.results = t2.pt.results = list()
for(i in 1:length(t2.mt)) {
  t2.mt.results[[i]] = hstats(t2.mt[[i]], usepopn=use.popn)
  res.text = c(res.text, paste(c("Table 2 mt family ", i, " moment estimate ", r(t2.mt.results[[i]]$moment.n.hat), " (", r(t2.mt.results[[i]]$moment.n.ci[1]), "-", r(t2.mt.results[[i]]$moment.n.ci[2]), "), fit estimate ", r(t2.mt.results[[i]]$fit.n.hat), " (", r(t2.mt.results[[i]]$fit.n.ci[1]), "-", r(t2.mt.results[[i]]$fit.n.ci[2]), ")"), collapse=""))
}
# get and store estimates for PT samples
for(i in 1:length(t2.pt)) {
  t2.pt.results[[i]] = hstats(t2.pt[[i]], usepopn=use.popn)
  res.text = c(res.text, paste(c("Table 2 pt family ", i, " moment estimate ", r(t2.pt.results[[i]]$moment.n.hat), " (", r(t2.pt.results[[i]]$moment.n.ci[1]), "-", r(t2.pt.results[[i]]$moment.n.ci[2]), "), fit estimate ", r(t2.pt.results[[i]]$fit.n.hat), " (", r(t2.pt.results[[i]]$fit.n.ci[1]), "-", r(t2.pt.results[[i]]$fit.n.ci[2]), ")"), collapse=""))
}
# perform and store LRT between organelles
t2.test = kimura_lrt(t2.mt, t2.pt)
res.text = c(res.text, paste(c("Table 2 LRT: MT ", r(t2.test$h1.n.hat), " (", r(t2.test$h1.n.ci.1.), "-", r(t2.test$h1.n.ci.2.), ") vs PT ", r(t2.test$h2.n.hat), " (", r(t2.test$h2.n.ci.1.), "-", r(t2.test$h2.n.ci.2.), ") p = ", t2.test$pval, ""), collapse=""))

## table 3 analysis -- within individuals in MSH1

# get and store estimates for MT samples
t3.mt.results = t3.pt.results = list()
for(i in 1:length(t3.mt)) {
  t3.mt.results[[i]] = hstats(t3.mt[[i]], usepopn=use.popn)
  res.text = c(res.text, paste(c("Table 3 mt individual ", i, " moment estimate ", r(t3.mt.results[[i]]$moment.n.hat), " (", r(t3.mt.results[[i]]$moment.n.ci[1]), "-", r(t3.mt.results[[i]]$moment.n.ci[2]), "), fit estimate ", r(t3.mt.results[[i]]$fit.n.hat), " (", r(t3.mt.results[[i]]$fit.n.ci[1]), "-", r(t3.mt.results[[i]]$fit.n.ci[2]), ")"), collapse=""))
}
# get and store estimates for PT samples
for(i in 1:length(t3.pt)) {
  t3.pt.results[[i]] = hstats(t3.pt[[i]], usepopn=use.popn)
  res.text = c(res.text, paste(c("Table 3 pt individual ", i, " moment estimate ", r(t3.pt.results[[i]]$moment.n.hat), " (", r(t3.pt.results[[i]]$moment.n.ci[1]), "-", r(t3.pt.results[[i]]$moment.n.ci[2]), "), fit estimate ", r(t3.pt.results[[i]]$fit.n.hat), " (", r(t3.pt.results[[i]]$fit.n.ci[1]), "-", r(t3.pt.results[[i]]$fit.n.ci[2]), ")"), collapse=""))
}
# perform and store LRTs:
t3.test.1 = kimura_lrt(t3.mt, t3.pt)	# MT vs PT within individuals
t3.test.2 = kimura_lrt(t3.mt, t2.mt)	# MT within individuals vs MT between generations
t3.test.3 = kimura_lrt(t3.pt, t2.pt)	# PT within individuals vs PT between generations
res.text = c(res.text, paste(c("Table 3 LRT: Tissue MT ", r(t3.test.1$h1.n.hat), " (", r(t3.test.1$h1.n.ci.1.), "-", r(t3.test.1$h1.n.ci.2.), ") vs Tissue PT ", r(t3.test.1$h2.n.hat), " (", r(t3.test.1$h2.n.ci.1.), "-", r(t3.test.1$h2.n.ci.2.), ") p = ", t3.test.1$pval, ""), collapse=""))
res.text = c(res.text, paste(c("Table 3 LRT: Tissue MT ", r(t3.test.2$h1.n.hat), " (", r(t3.test.2$h1.n.ci.1.), "-", r(t3.test.2$h1.n.ci.2.), ") vs Gen MT ", r(t3.test.2$h2.n.hat), " (", r(t3.test.2$h2.n.ci.1.), "-", r(t3.test.2$h2.n.ci.2.), ") p = ", t3.test.2$pval, ""), collapse=""))
res.text = c(res.text, paste(c("Table 3 LRT: Tissue PT ", r(t3.test.3$h1.n.hat), " (", r(t3.test.3$h1.n.ci.1.), "-", r(t3.test.3$h1.n.ci.2.), ") vs Gen PT ", r(t3.test.3$h2.n.hat), " (", r(t3.test.3$h2.n.ci.1.), "-", r(t3.test.3$h2.n.ci.2.), ") p = ", t3.test.3$pval, ""), collapse=""))

## table 4 analysis -- between generations in WT background

# get and store estimates for MT samples
t4.mt.results = list()
for(i in 1:length(t4.mt)) {
  t4.mt.results[[i]] = hstats(t4.mt[[i]], usepopn=use.popn)
  res.text = c(res.text, paste(c("Table 4 mt family ", i, " moment estimate ", r(t4.mt.results[[i]]$moment.n.hat), " (", r(t4.mt.results[[i]]$moment.n.ci[1]), "-", r(t4.mt.results[[i]]$moment.n.ci[2]), "), fit estimate ", r(t4.mt.results[[i]]$fit.n.hat), " (", r(t4.mt.results[[i]]$fit.n.ci[1]), "-", r(t4.mt.results[[i]]$fit.n.ci[2]), ")"), collapse=""))
}
# perform and store LRT between WT and MSH1 backgrounds
t4.test = kimura_lrt(t2.mt, t4.mt)
res.text = c(res.text, paste(c("Table 4 LRT: MSH1 ", r(t4.test$h1.n.hat), " (", r(t4.test$h1.n.ci.1.), "-", r(t4.test$h1.n.ci.2.), ") vs WT ", r(t4.test$h2.n.hat), " (", r(t4.test$h2.n.ci.1.), "-", r(t4.test$h2.n.ci.2.), ") p = ", t4.test$pval, ""), collapse=""))

# sanity checks -- should fail to reject null, and get reasonable parameter estimates
kimura_lrt(t2.mt[1], t2.mt[2])
kimura_lrt(t3.mt[2], t3.mt[3]) # 1 is an outlier
maxlik(rkimura(10, 0.5, 0.5))
maxlik(rkimura(100, 0.5, 0.5))
maxlik(rkimura(10, 0.2, 0.5))
maxlik(rkimura(20, 0.7, 0.9))

#### redo with outlier/awkward case removal
t2.pt.prime = t2.pt[-c(1,2)]
t3.mt.prime = t3.mt[-1]
t4.mt.prime = t4.mt[-1]
t2.test.prime = kimura_lrt(t2.mt, t2.pt.prime)
t3.test.1.prime = kimura_lrt(t3.mt.prime, t3.pt)	
t3.test.2.prime = kimura_lrt(t3.mt.prime, t2.mt)
t4.test.prime = kimura_lrt(t2.mt, t4.mt.prime)	

res.text = c(res.text, paste(c("Table 2* LRT: MT ", r(t2.test.prime$h1.n.hat), " (", r(t2.test.prime$h1.n.ci.1.), "-", r(t2.test.prime$h1.n.ci.2.), ") vs PT ", r(t2.test.prime$h2.n.hat), " (", r(t2.test.prime$h2.n.ci.1.), "-", r(t2.test.prime$h2.n.ci.2.), ") p = ", t2.test.prime$pval, ""), collapse=""))
res.text = c(res.text, paste(c("Table 3* LRT: Tissue MT ", r(t3.test.1.prime$h1.n.hat), " (", r(t3.test.1.prime$h1.n.ci.1.), "-", r(t3.test.1.prime$h1.n.ci.2.), ") vs Tissue PT ", r(t3.test.1.prime$h2.n.hat), " (", r(t3.test.1.prime$h2.n.ci.1.), "-", r(t3.test.1.prime$h2.n.ci.2.), ") p = ", t3.test.1.prime$pval, ""), collapse=""))
res.text = c(res.text, paste(c("Table 3* LRT: Tissue MT ", r(t3.test.2.prime$h1.n.hat), " (", r(t3.test.2.prime$h1.n.ci.1.), "-", r(t3.test.2.prime$h1.n.ci.2.), ") vs Gen MT ", r(t3.test.2.prime$h2.n.hat), " (", r(t3.test.2.prime$h2.n.ci.1.), "-", r(t3.test.2.prime$h2.n.ci.2.), ") p = ", t3.test.2.prime$pval, ""), collapse=""))
res.text = c(res.text, paste(c("Table 4* LRT: MSH1 ", r(t4.test.prime$h1.n.hat), " (", r(t4.test.prime$h1.n.ci.1.), "-", r(t4.test.prime$h1.n.ci.2.), ") vs WT MT ", r(t4.test.prime$h2.n.hat), " (", r(t4.test.prime$h2.n.ci.1.), "-", r(t4.test.prime$h2.n.ci.2.), ") p = ", t4.test.prime$pval, ""), collapse=""))

# output results
res.text
