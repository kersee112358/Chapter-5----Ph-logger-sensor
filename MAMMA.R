

####################################################################################
# A Compilation of R Functions Used for Data Analyses Presented in the Boook
# Measuring Agreement: Models, Methods, and Applications 
# (By P. K. Choudhary and H. N. Nagaraja)
# 11/13/2017
# Note: This file is called by R code for each data analysis as source("MAMMA.R")
####################################################################################

require(nlme)
require(numDeriv)
require(mvtnorm)
require(lattice)
require(latticeExtra)
require(multcomp)
require(Matrix)

####################################################################################
# chapter 12 #
####################################################################################

################################
# 1. Bangdiwala chart function #
################################

# Chaitra H. Nagaraja
# 11/23/2015
# updated 6/23/2017

# constructed using "The agreement chart" by Shrikant I Bangdiwala
# http://www.stat.ncsu.edu/information/library/mimeo.archive/ISMS_1988_1859.pdf

# Also used agreementplot() function in vcd R package

bangdiwala.chart <- function(x, main="", 
                             xlab=names(dimnames(x))[2], 
                             ylab= names(dimnames(x))[1], 
                             cat.label=colnames(x),
                             las=TRUE,
                             add.legend=TRUE,	
                             location.legend="bottomright",	
                             add.B.hat=TRUE	,
                             location.B.hat="topleft",
                             add.x.margin=TRUE,
                             add.y.margin=TRUE,
                             ordinal.data=TRUE,						  
                             col.line="grey"
                             ,
                             lwd.line=2,
                             lty.line=1){
  
  # x: data as a table object in R
  
  # note: columns of x will turn into x-axis
  #       rows of x will turn into y-axis
  #       to reverse, use t() function (i.e., transpose)
  
  # function arguments for labels:
  # main: title of plot
  # xlab, ylab: title of x-axis and y-axis respectively
  #			  (default: extracted from x)
  # cat.label: category labels
  #			  (default: extracted from x)	
  # las: y-axis category labels rotated
  # add.legend: TRUE/FALSE, add a legend to the plot
  # location.legend: either "bottomright" or "topleft" 
  #		      (only relevant if add.legend=TRUE)
  # add.B.hat: TRUE/FALSE, add Bangdiwala agreement measure to plot
  #			 (otherwise it is only given as the function's output)
  #            (agreement measure without weights)
  # location.B.hat: either "bottomright" or "topleft"
  #			  (only relevant if add.B.hat=TRUE)
  # add.x.margin/add.y.margin: TRUE/FALSE, add marginal values for 
  #            columns (x) or rows (y) to plot
  
  # function arguments for rectangles on chart:
  # ordinal.data: TRUE/FALSE, if true, partial agreement squares are added
  
  # function arguments for 45 degree line:
  # col.line: color of line (default: "firebrick")
  # lwd.line: width of line (default: 2)
  # lty.line: type of line (solid, dotted, etc.; default: solid)
  
  #####
  # 0. preliminary checks and calculations
  
  # ensure data is of the correct type:
  if(any(!is.numeric(x))) return("ERROR: some cells have non-numeric counts")
  if(any(is.na(x))) return("ERROR: some cells have missing counts")
  if(any(x < 0))  return("ERROR: some cells have negative counts")
  
  # ensure data is in the form of a table:
  if(!is.table(x)) x <- as.table(x)
  
  n <- sum(x)
  k <- length(cat.label)
  
  #####
  # 1. draw n x n square
  
  plot(1:n, 1:n, col="white", main=main, xlab=xlab, ylab=ylab,  
       xlim=c(0,n), ylim=c(0,n), xaxs="i", yaxs="i", xaxt="n", yaxt="n")
  
  #####
  # 2. add axis labels 
  x.margin <- c(0, cumsum(as.vector(apply(x, 2, sum))))
  y.margin <- c(0, cumsum(as.vector(apply(x, 1, sum))))
  # axis labels should be at the center of each rectangle in step 3
  axis(side=1, at=(x.margin[1:k] + x.margin[2:(k+1)])/2, labels=cat.label)
  axis(side=2, at=(y.margin[1:k] + y.margin[2:(k+1)])/2, labels=cat.label, las=las)
  
  #####	
  # 3. k rectangles with width = column marginals, height = row marginals
  for(i in 1:k){
    rect(xleft=x.margin[i], ybottom=y.margin[i], xright=x.margin[i+1], ytop=y.margin[i+1])
  } # end for loop
  rm(i)	
  
  #####
  # 4. partial agreement squares Using Bangdiwala partial agreement definition (1988, p 11)
  #    -- add only if your rating categories are ordinal
  
  if(ordinal.data){
    
    # darker gray --> closer to full agreement 
    gray.vector <- gray.colors(n=k-2, start=0.3, end=0.6)
    
    for(j in rev(1:(k-2))){ # loop through colors
      
      # within +/- j of correct rating (make larger rectangles first)
      for(i in 1:k){ # loop through rating values
        
        # summation vector
        temp <- sort((i-j):(i+j), decreasing=FALSE)
        temp <- temp[ (temp > 0) & (temp <= k) ]
        
        partial.agreement.x <- sum(as.vector(x[temp,i]))
        partial.agreement.y <- sum(as.vector(x[i,temp]))
        
        if(temp[1]==1){
          add.x <- 0
          add.y <- 0
        }else{
          add.x <- sum(as.vector(x[1:(temp[1]-1),i]))
          add.y <- sum(as.vector(x[i,1:(temp[1]-1)]))
        } # end if/else
        
        rect(xleft=x.margin[i]+add.x, ybottom=y.margin[i]+add.y, 
             xright=x.margin[i]+add.x+partial.agreement.x, 
             ytop=y.margin[i]+add.y+partial.agreement.y, 
             col=gray.vector[j])
        
        rm(temp, partial.agreement.x, partial.agreement.y, add.x, add.y)
      } # end inner for loop
      
    } # end outer for loop
  } # end if for ordinal.data
  
  #####
  # 5. full agreement squares (in black)
  full.agreement <- as.vector(diag(x))
  for(i in 1:k){
    if(i==1){
      add.x <- 0
      add.y <- 0
    }else{
      add.x <- sum(as.vector(x[1:(i-1),i]))
      add.y <- sum(as.vector(x[i,1:(i-1)]))
    } # end if/else
    
    rect(xleft=x.margin[i]+add.x, ybottom=y.margin[i]+add.y, 
         xright=x.margin[i]+add.x+full.agreement[i], ytop=y.margin[i]+add.y+full.agreement[i], 
         col="black")
    rm(add.x, add.y)
  } # end for loop
  rm(i)
  
  ####
  # 6. draw 45 degree diagonal line 
  abline(a=0, b=1, col=col.line, lwd=lwd.line, lty=lty.line)
  
  ####
  # 7. add marginal counts to plot?
  if(add.x.margin){
    temp <- as.vector(apply(x, 2, sum))
    # axis labels should be at the center of each rectangle in step 3
    axis(side=3, at=(x.margin[1:k] + x.margin[2:(k+1)])/2, labels=temp, 
         cex.axis=0.8, line=-0.8, tick=FALSE)	
    rm(temp)	
  } # end if
  
  if(add.y.margin){		
    temp <- as.vector(apply(x, 1, sum))
    # axis labels should be at the center of each rectangle in step 3		
    axis(side=4, at=(y.margin[1:k] + y.margin[2:(k+1)])/2, labels=temp, las=las, 
         cex.axis=0.8, line=-0.8, tick=FALSE)
    rm(temp)
  } # end if
  
  ####
  # 8. add legend?
  
  if(add.legend){
    
    legend(location.legend, 
           legend=c("full", paste("partial (+/- ",1:(k-1), ")", sep=""), 
                    # PKC:                    expression(paste(45*degree, " line", sep=""))),
                    paste(45, " degree line", sep="")),
           col=c("black", gray.vector, "black", col.line),
           pch=c(rep(15, times=k-1),22, NA),
           pt.cex=2,
           lty=c(rep(NA, times=k), lty.line),
           lwd=c(rep(NA, times=k), lwd.line), bty="n")
    
  } # end if statement
  
  ####	
  # 9. compute Bangdiwala agreement measure 	
  B.hat <- sum(full.agreement^2)/sum((as.vector(apply(x, 2, sum))*as.vector(apply(x, 1, sum))))
  # note: don't have to get rid of the first element of x.margin or y.margin because it is 0
  
  # add B.hat to plot? 
  if(add.B.hat){
    legend(location.B.hat, legend=bquote(hat(B)==.(round(B.hat, digits=3))), 
           fill=NA, bty="n", border="white")
  } # end if
  
  # output Bangdiwala agreement measure
  return(B.hat)
  
} # end function


####################################################################################
# chapter 11 #
####################################################################################

#########################################################################
# compute true values of agreement measures	in case of unlinked repeated
# measurements data					       #
#########################################################################

true.measures <- function(param, prob) {
  q.chisq <- qchisq(prob, 1)
  
  ltdi <- ltdi.unlinked.fun(param, prob)
  ltdi.1 <- (1/2) * (log(2 * q.chisq) + param["lsigmasq.e1"])
  ltdi.2 <- (1/2) * (log(2 * q.chisq) + param["lsigmasq.e2"])	
  
  zccc <- zccc.unlinked.fun(param)
  zccc.1 <- zccc.1.unlinked.fun(param)
  zccc.2 <- zccc.2.unlinked.fun(param)
  
  results <- c(exp(c(ltdi, ltdi.1, ltdi.2)), tanh(c(zccc, zccc.1, zccc.2)))
  names(results) <- c("tdi", "tdi.1", "tdi.2", "ccc", "ccc.1", "ccc.2")
  return(results)
}

###########################################################################
# simulate data and compute confidence bounds in case of unlinked repeated
# measurements data						       #
###########################################################################

simulate.bounds <- function() {  
  # simulate random effects
  bsim <- rnorm(n, mean = mu.b, sd = sqrt(exp(lsigmasq.b)))
  b1sim <- rnorm(n, sd = sqrt(exp(lsigmasq.bI)))
  b2sim <- rnorm(n, sd = sqrt(exp(lsigmasq.bI)))
  
  # simulate errors
  e1sim <- rnorm(n*m, sd = sqrt(exp(lsigmasq.e1)))
  e2sim <- rnorm(n*m, sd = sqrt(exp(lsigmasq.e2)))
  
  # get simulated data
  y1sim <- rep(bsim, each = m) + rep(b1sim, each = m) + e1sim
  y2sim <- alpha + rep(bsim, each = m) + rep(b2sim, each = m) + e2sim
  ysim <- c(y1sim, y2sim)
  
  # create a new dataframe
  subject <- as.factor(rep(rep(1:n, each = m), 2))
  method <- as.factor(c(rep("method1", n*m), rep("method2", n*m)))
  rep <- rep(rep(1:m, n), 2)
  
  sim.dat  <-  data.frame(subject, method, rep, ysim)
  colnames(sim.dat) <- c("subject","method","rep","meas")
  
  # fit model and get confidence bounds
  fit.results <- lme.unlinked.repmeas.fit(sim.dat)
  conf.results <- conf.measures.unlinked(sim.dat, fit.results$param.hat, conf.level, prob)
  
  results <- c(conf.results$tdi, conf.results$rep.tdi[1, ],
               conf.results$rep.tdi[2, ], conf.results$ccc, 
               conf.results$rep.ccc[1, ], conf.results$rep.ccc[2, ])
  names(results) <- c("ltdi", "se.ltdi", "ucb.ltdi", "tdi", "ucb.tdi", 
                      "ltdi.1", "se.ltdi.1", "ucb.ltdi.1", "tdi.1", "ucb.tdi.1", 
                      "ltdi.2", "se.ltdi.2", "ucb.ltdi.2", "tdi.2", "ucb.tdi.2",
                      "zccc", "se.zccc", "lcb.zccc", "ccc", "lcb.ccc", 
                      "zccc.1", "se.zccc.1", "lcb.zccc.1", "ccc.1", "lcb.ccc.1",
                      "zccc.2", "se.zccc.2", "lcb.zccc.2", "ccc.2", "lcb.ccc.2")
  
  return(results)
}


#########################################################################
# compute true values of agreement measures	in case of bivariate normal
# data					       #
#########################################################################

true.measures.bvn <- function(param, prob) {
  ltdi <- ltdi.fun(param, prob, "bvn")
  zccc <- zccc.fun(param, "bvn")
  results <- c(exp(ltdi), tanh(zccc))
  names(results) <- c("tdi", "ccc")
  return(results)
}


###########################################################################
# simulate data and compute confidence bounds in case of bivariate
# normal data						       #
###########################################################################

simulate.bounds.bvn <- function() {  
  require("mvtnorm")
  # simulate data
  cov.mat <- matrix(c(sigmasq.1, rho * sqrt(sigmasq.1 * sigmasq.2), rho * 
                        sqrt(sigmasq.1 * sigmasq.2), sigmasq.2), byrow = T, ncol = 2)
  
  sim.dat <- rmvnorm(n, mean = c(mu.1, mu.2), sigma = cov.mat)
  colnames(sim.dat) <- c("method1", "method2")
  
  # fit model and get confidence bounds
  fit.results <- bvn.fit(sim.dat)
  conf.results <- conf.measures(sim.dat, fit.results$param.hat, "bvn", "no", 
                                conf.level, prob)
  
  results <- c(conf.results$tdi, conf.results$ccc)
  names(results) <- c("ltdi", "se.ltdi", "ucb.ltdi", "tdi", "ucb.tdi", 
                      "zccc", "se.zccc", "lcb.zccc", "ccc", "lcb.ccc")
  return(results)
}


####################################################################################
# chapter 10 #
####################################################################################


######################
# Function for computing estimates and confidence bounds and intervals 
# for performing nonparametric inference on measures for J = 2 or 3 methods
# assuming balanced design and all pairwise comparisons
# Arguments: 
#	dat = a data frame with four columns --- subject (a factor), 
#			  method (a factor with 3 levels ), 
#			rep (replication number, unused) and meas
#	alpha = level of confidence
#	prob = probability cutoff for TDI
# 	comparison.type = type of multiple comparisons (only "all.pairwise" comparisons
#				implemented at present)
# Value:
#		A bunch of lists (need to specify)

npm.inference <- function(dat, alpha = 0.05, prob = 0.9, comparison.type = "all.pairwise") {
  # drop unused levels
  dat <- droplevels(dat)
  # some initial checks 
  if (comparison.type != "all.pairwise") {
    stop("invalid multiple comparisons type")
  }
  
  J <- nlevels(dat$method)
  method.names <- levels(dat$method)
  
  if (!(J %in% c(2, 3))) {
    stop("number of methods is not either 2 or 3")
  }
  
  nsub <- nlevels(dat$subject)
  nrep <- table(dat$subject, dat$method)
  
  if (J == 2) {
    if (mean((nrep[, 1] == nrep[, 2])) != 1) {
      stop("design is not balanced")
    }
    nrep <- unique(nrep)
    wt <- 1/(nsub * (nrep^2))
    wt.paired <- rep(wt, nsub * (nrep^2))
    results <- conf.measures.2methods.npm(dat, nsub, nrep, wt.paired, method.names, 
                                          alpha, prob)
  }
  
  if (J == 3) {
    if (mean((nrep[, 1] == nrep[, 2]) & (nrep[, 2] == nrep[, 3])) != 1) {
      stop("design is not balanced")
    }
    nrep <- unique(nrep)
    wt <- 1/(nsub * (nrep^2))
    wt.paired <- rep(wt, nsub * (nrep^2))
    results <- conf.measures.3methods.npm(dat, nsub, nrep, wt.paired, method.names, 
                                          alpha, prob)
  }
  return(results)
}

################
#
# Function for computing nonparametric estimates for measures and their 
# appropriate confidence intervals for J = 2 methods 
# Arguments: 
#	dat = a data frame with four columns --- subject (a factor), 
#			  method (a factor with 3 levels ), 
#			rep (replication number, unused) and meas
#	nsub = number of subjects
#	nrep = common number of replications for balanced design
# 	wt.paired = vector of length same as nrow(dat) containing weights for paired
#				measurements
#	method.names = a vector of length 2 containing names of the two methods
#	alpha = level of confidence
#	prob = probabilty cutoff for TDI
# Value:
#		A bunch of lists (need to specify)
#

conf.measures.2methods.npm <- function(dat, nsub, nrep, wt.paired, method.names, 
                                       alpha, prob) {
  
  stopifnot(nlevels(dat$method) == 2)
  
  # split the data by subject, form pairs for each method pair, and get estimates
  dat.split <- split(dat, dat$subject)
  
  # methods 1 and 2
  dat12 <- data.frame(do.call("rbind", lapply(dat.split, make.pairs.unlinked, 
                                              nrep = nrep, method.pair = method.names[c(1, 2)])))
  est12 <- plugin.estimates(dat12, wt.paired, prob)
  
  # moment estimates
  moments.hat <- c(est12$est["mean1"], est12$est["mean2"], est12$est["sd1"], 
                   est12$est["sd2"], est12$est["corr"])
  names(moments.hat) <- c(paste0("mean", method.names), paste0("sd", method.names), 
                          paste0("corr", method.names[1], method.names[2]))
  
  # estimate and lower bound for ccc
  ccc.hat <- est12$est["ccc"]
  ccc.var <- diag.Omega(est12$infl$ccc, dat12$subject, nsub, nrep)/nsub
  ci.zccc <- ci.simul.paramfun(ccc.hat, ccc.var, atanh, diag(1), "greater", 
                               1 - alpha)
  ci.zccc <- list(ci.fun = ci.zccc$ci.fun, zresult = ci.zccc$result, result = tanh(ci.zccc$result[, 
                                                                                                  c("estimate", "lcl", "ucl"), drop = FALSE]))
  
  # estimate and upper bound for cp at estimated TDI
  cp.hat <- est12$est["cp"]
  cp.var <- diag.Omega(est12$infl$cp, dat12$subject, nsub, nrep)/nsub
  ci.cp <- ci.simul.paramfun(cp.hat, cp.var, identity, diag(1), "less", 1 - 
                               alpha)
  
  # estimate and upper bound for TDI
  tdi.hat <- est12$est["tdi"]
  cp.cval <- attr(ci.cp$ci.fun$confint, "calpha")
  new.prob <- min(prob + cp.cval * sqrt(cp.var), 1)
  ucb.tdi <- tdi.fun.npm(dat12, wt.paired, new.prob)
  ci.tdi <- list(tdi.hat = tdi.hat, ucb.tdi = ucb.tdi)
  
  # estimate and interval for mean difference
  meandiff.hat <- est12$est["meandiff"]
  meandiff.var <- diag.Omega(est12$infl$meandiff, dat12$subject, nsub, nrep)/nsub
  ci.meandiff <- ci.simul.paramfun(meandiff.hat, meandiff.var, identity, 
                                   diag(1), "two.sided", 1 - alpha)
  
  # estimate and interval for variance ratio
  varratio.hat <- est12$est["varratio"]
  varratio.var <- diag.Omega(est12$infl$varratio, dat12$subject, nsub, nrep)/nsub
  ci.lvarratio <- ci.simul.paramfun(varratio.hat, varratio.var, log, diag(1), 
                                    "two.sided", 1 - alpha)
  ci.lvarratio <- list(ci.fun = ci.lvarratio$ci.fun, lresult = ci.lvarratio$result, 
                       result = exp(ci.lvarratio$result[, c("estimate", "lcl", "ucl"), drop = FALSE]))
  
  results <- list(moments = moments.hat, mean.diff = ci.meandiff, varratio = ci.lvarratio, 
                  cp = ci.cp, tdi = ci.tdi, ccc = ci.zccc)
  return(results)
  
}


################
#
# Function for computing nonparametric estimates for measures and their 
# appropriate simultaneous confidence intervals for J = 3 methods 
# Arguments: 
#	dat = a data frame with four columns --- subject (a factor), 
#			  method (a factor with 3 levels ), 
#			rep (replication number, unused) and meas
#	nsub = number of subjects
#	nrep = common number of replications for balanced design
# 	wt.paired = vector of length same as nrow(dat) containing weights for paired
#				measurements
#	method.names = a vector of length 3 containing names of the three methods
#	alpha = level of confidence
#	prob = probabilty cutoff for TDI
# Value:
#		A bunch of lists (need to specify)
#

conf.measures.3methods.npm <- function(dat, nsub, nrep, wt.paired, method.names, 
                                       alpha, prob) {
  
  stopifnot(nlevels(dat$method) == 3)
  
  # split the data by subject, form pairs for each method pair, and get estimates
  dat.split <- split(dat, dat$subject)
  
  # methods 1 and 2
  dat12 <- data.frame(do.call("rbind", lapply(dat.split, make.pairs.unlinked, 
                                              nrep = nrep, method.pair = method.names[c(1, 2)])))
  est12 <- plugin.estimates(dat12, wt.paired, prob)
  
  # methods 1 and 3
  dat13 <- data.frame(do.call("rbind", lapply(dat.split, make.pairs.unlinked, 
                                              nrep = nrep, method.pair = method.names[c(1, 3)])))
  est13 <- plugin.estimates(dat13, wt.paired, prob)
  
  # methods 2 and 3
  dat23 <- data.frame(do.call("rbind", lapply(dat.split, make.pairs.unlinked, 
                                              nrep = nrep, method.pair = method.names[c(2, 3)])))
  est23 <- plugin.estimates(dat23, wt.paired, prob)
  
  # moment estimates
  moments.hat <- c(est12$est["mean1"], est23$est["mean1"], est23$est["mean2"], 
                   est12$est["sd1"], est23$est["sd1"], est23$est["sd2"], est12$est["corr"], 
                   est13$est["corr"], est23$est["corr"])
  names(moments.hat) <- c(paste0("mean", method.names), paste0("sd", method.names), 
                          paste0("corr", method.names[1], method.names[2]), paste0("corr", method.names[1], 
                                                                                   method.names[3]), paste0("corr", method.names[2], method.names[3]))
  
  # estimates and simultaneous lower bounds for ccc
  contrast.matrix <- diag(3)
  rownames(contrast.matrix) <- c(paste0(method.names[1], method.names[2]), 
                                 paste0(method.names[1], method.names[3]), paste0(method.names[2], method.names[3]))
  colnames(contrast.matrix) <- rownames(contrast.matrix)
  
  ccc.hat <- c(est12$est["ccc"], est13$est["ccc"], est23$est["ccc"])
  ccc.covmat <- Omega.matrix(est12$infl$ccc, est13$infl$ccc, est23$infl$ccc, 
                             dat12$subject, dat13$subject, dat23$subject, nsub, nrep)/nsub
  ci.zccc <- ci.simul.paramfun(ccc.hat, ccc.covmat, atanh, contrast.matrix, 
                               "greater", 1 - alpha)
  ci.zccc <- list(ci.fun = ci.zccc$ci.fun, zresult = ci.zccc$result, result = tanh(ci.zccc$result[, 
                                                                                                  c("estimate", "lcl", "ucl"), drop = FALSE]))
  
  # estimates and simultaneous upper bounds for cp at estimated TDI
  cp.hat <- c(est12$est["cp"], est13$est["cp"], est23$est["cp"])
  cp.covmat <- Omega.matrix(est12$infl$cp, est13$infl$cp, est23$infl$cp, 
                            dat12$subject, dat13$subject, dat23$subject, nsub, nrep)/nsub
  ci.cp <- ci.simul.paramfun(cp.hat, cp.covmat, identity, contrast.matrix, 
                             "less", 1 - alpha)
  
  # estimates and simultaneous upper bounds for TDI
  tdi.hat <- c(est12$est["tdi"], est13$est["tdi"], est23$est["tdi"])
  names(tdi.hat) <- rownames(contrast.matrix)
  cp.cval <- attr(ci.cp$ci.fun$confint, "calpha")
  new.prob <- pmin(prob + cp.cval * sqrt(diag(cp.covmat)), rep(1, 3))
  ucb.tdi <- c(tdi.fun.npm(dat12, wt.paired, new.prob[1]), tdi.fun.npm(dat13, 
                                                                       wt.paired, new.prob[2]), tdi.fun.npm(dat23, wt.paired, new.prob[3]))
  names(ucb.tdi) <- names(tdi.hat)
  ci.tdi <- list(tdi.hat = tdi.hat, ucb.tdi = ucb.tdi)
  
  # estimates and simultaneous intervals for mean difference
  meandiff.hat <- c(est12$est["meandiff"], est13$est["meandiff"], est23$est["meandiff"])
  meandiff.covmat <- Omega.matrix(est12$infl$meandiff, est13$infl$meandiff, 
                                  est23$infl$meandiff, dat12$subject, dat13$subject, dat23$subject, nsub, nrep)/nsub
  ci.meandiff <- ci.simul.paramfun(meandiff.hat, meandiff.covmat, identity, 
                                   contrast.matrix, "two.sided", 1 - alpha)
  
  # estimates and simultaneous intervals for variance ratio
  varratio.hat <- c(est12$est["varratio"], est13$est["varratio"], est23$est["varratio"])
  varratio.covmat <- Omega.matrix(est12$infl$varratio, est13$infl$varratio, 
                                  est23$infl$varratio, dat12$subject, dat13$subject, dat23$subject, nsub, nrep)/nsub
  ci.lvarratio <- ci.simul.paramfun(varratio.hat, varratio.covmat, log, contrast.matrix, 
                                    "two.sided", 1 - alpha)
  ci.lvarratio <- list(ci.fun = ci.lvarratio$ci.fun, lresult = ci.lvarratio$result, 
                       result = exp(ci.lvarratio$result[, c("estimate", "lcl", "ucl"), drop = FALSE]))
  
  results <- list(moments = moments.hat, mean.diff = ci.meandiff, varratio = ci.lvarratio, 
                  cp = ci.cp, tdi = ci.tdi, ccc = ci.zccc)
  
  return(results)
}


################
#
# Function for computing the Omega matrix 
# Arguments: 
#		infl12, infl13, infl23 = vectors of estimated influence
#		subid12, subid13, subid23 = vectors of corresponding subject ids
#		nrep = common number of replications for a balanced design
# Value:
#		the 3x3 Omega matrix
#

Omega.matrix <- function(infl12, infl13, infl23, subid12, subid13, subid23, 
                         nsub, nrep) {
  diag12 <- diag.Omega(infl12, subid12, nsub, nrep)
  diag13 <- diag.Omega(infl13, subid13, nsub, nrep)
  diag23 <- diag.Omega(infl23, subid23, nsub, nrep)
  offdiag <- offdiag.Omega(infl12, infl13, infl23, subid12, subid13, subid23, 
                           nsub, nrep)
  Omega <- matrix(NA, nrow = 3, ncol = 3)
  Omega[, 1] <- c(diag12, offdiag["12.13"], offdiag["12.23"])
  Omega[, 2] <- c(offdiag["12.13"], diag13, offdiag["13.23"])
  Omega[, 3] <- c(offdiag["12.23"], offdiag["13.23"], diag23)
  rownames(Omega) <- c("12.13", "12.23", "13.23")
  colnames(Omega) <- c("12.13", "12.23", "13.23")
  return(Omega)
}

################
#
# Function for computing off diagonal elements of the Omega matrix 
# Arguments: 
#		infl12, infl13, infl23 = vectors of estimated influence
#		subid12, subid13, subid23 = vectors of corresponding subject ids
#		nrep = common number of replications for a balanced design
# Value:
#		a vector containing the three off diagonal elements
#
offdiag.Omega <- function(infl12, infl13, infl23, subid12, subid13, subid23, 
                          nsub, nrep) {
  infl12 <- split(infl12, subid12)
  infl13 <- split(infl13, subid13)
  infl23 <- split(infl23, subid23)
  
  ss <- rowSums(mapply(FUN = "inner.sums.offdiag.Omega", x12 = infl12, x13 = infl13, 
                       x23 = infl23, nrep = nrep))
  
  # 12 vs 13
  cov1 <- ss["sum12.13"]/(nsub * (nrep^3))
  if (nrep == 1) {
    cov12.13 <- cov1/nrep
  }
  if (nrep > 1) {
    cov2 <- (ss["asum12.13"] - ss["sum12.13"])/(nsub * (nrep^3) * (nrep - 
                                                                     1))
    cov12.13 <- (cov1 + (nrep - 1) * cov2)/nrep
  }
  
  # 12 vs 23
  cov1 <- ss["sum12.23"]/(nsub * (nrep^3))
  if (nrep == 1) {
    cov12.23 <- cov1/nrep
  }
  if (nrep > 1) {
    cov2 <- (ss["asum12.23"] - ss["sum12.23"])/(nsub * (nrep^3) * (nrep - 
                                                                     1))
    cov12.23 <- (cov1 + (nrep - 1) * cov2)/nrep
  }
  
  # 13 vs 23
  cov1 <- ss["sum13.23"]/(nsub * (nrep^3))
  if (nrep == 1) {
    cov13.23 <- cov1/nrep
  }
  if (nrep > 1) {
    cov2 <- (ss["asum13.23"] - ss["sum13.23"])/(nsub * (nrep^3) * (nrep - 
                                                                     1))
    cov13.23 <- (cov1 + (nrep - 1) * cov2)/nrep
  }
  
  result <- c(cov12.13, cov12.23, cov13.23)
  names(result) <- c("12.13", "12.23", "13.23")
  return(result)
}

################
#
# Function for computing a diagonal element of the Omega matrix 
# Arguments: 
#		infl = vector of estimated influence
#		subid = vector of subject ids corresponding to the influences
#		nrep = common number of replications for a balanced design
# Value:
#		a scalar giving the desired diagonal element
#

diag.Omega <- function(infl, subid, nsub, nrep) {
  sums <- rowSums(sapply(split(infl, subid), FUN = "inner.sums.diag.Omega", 
                         nrep = nrep))
  t1 <- sums[1]/(nsub * (nrep^2))
  if (nrep == 1) {
    avar <- t1/(nrep^2)
  }
  if (nrep > 1) {
    t2 <- sums[2]/(nsub * (nrep^2) * (nrep - 1))
    t3 <- sums[3]/(nsub * (nrep^2) * (nrep - 1))
    t4 <- sums[4]/(nsub * (nrep^2) * ((nrep - 1)^2))
    avar <- (t1 + (nrep - 1) * t2 + (nrep - 1) * t3 + ((nrep - 1)^2) * 
               t4)/(nrep^2)
  }
  return(unname(avar))
}
################
#
# Function for computing the inner sums needed for computing off diagonal elements of 
# the Omega matrix with J = 3
# Arguments: 
#		x12, x13, x23 = vectors of estimated influences for a subject
#		nrep = common number of replications for a balanced design
# Value:
#		a vector containing six appropriate inner sums	
#

inner.sums.offdiag.Omega <- function(x12, x13, x23, nrep) {
  x12 <- matrix(x12, nrow = nrep, ncol = nrep, byrow = T)
  x13 <- matrix(x13, nrow = nrep, ncol = nrep, byrow = T)
  x23 <- matrix(x23, nrow = nrep, ncol = nrep, byrow = T)
  sum12.13 <- sum(rowSums(x12) * rowSums(x13))
  sum12.23 <- sum(colSums(x12) * rowSums(x23))
  sum13.23 <- sum(colSums(x13) * colSums(x23))
  asum12.13 <- sum(x12) * sum(x13)
  asum12.23 <- sum(x12) * sum(x23)
  asum13.23 <- sum(x13) * sum(x23)
  return(c(sum12.13 = sum12.13, sum12.23 = sum12.23, sum13.23 = sum13.23, 
           asum12.13 = asum12.13, asum12.23 = asum12.23, asum13.23 = asum13.23))
}


################
#
# Function for computing the four inner sums needed for computing diagonal elements of 
# the Omega matrix
# Arguments: 
#		x = vector of estimated influences for a subject
#		nrep = common number of replications for a balanced design
# Value:
#		a vector containing four inner sums	
#

inner.sums.diag.Omega <- function(x, nrep) {
  x <- matrix(x, nrow = nrep, ncol = nrep, byrow = T)
  a1 <- sum(x^2)
  a2 <- sum(rowSums(x)^2)
  a3 <- sum(colSums(x)^2)
  a4 <- sum(x)^2
  return(c(sum1 = a1, sum2 = a2 - a1, sum3 = a3 - a1, sum4 = a4 - a3 - a2 + 
             a1))
}

######################
#
# Function for computing plug-in estimates of measures and their influence functions
# Arguments: 
#		paired.dat = a data frame containing all possible pairs
#		wt.paired = a vector of same length as paired.dat containing the weights
#					for each observation pair
#		prob = probability cutoff for tdi 
# Value:
#		a list of two items --- one a vector of various estimates, and a list 
#		of estimated influence vectors computed using paired.dat 		
#

plugin.estimates <- function(paired.dat, wt.paired, prob) {
  
  x1 <- paired.dat$x1
  x2 <- paired.dat$x2
  absd <- abs(x2 - x1)
  
  mean1 <- drop(x1 %*% wt.paired)
  mean2 <- drop(x2 %*% wt.paired)
  var1 <- drop(x1^2 %*% wt.paired) - mean1^2
  var2 <- drop(x2^2 %*% wt.paired) - mean2^2
  cov12 <- drop((x1 * x2) %*% wt.paired) - (mean1 * mean2)
  
  x1dev <- x1 - mean1
  x2dev <- x2 - mean2
  
  # tdi at prob
  tdi <- tdi.fun.npm(paired.dat, wt.paired, prob)
  
  # cp at tdi
  cp <- cp.fun.npm(paired.dat, wt.paired, tdi)
  cp.infl <- 1 * (absd <= tdi) - cp
  
  # ccc
  denom <- var1 + var2 + (mean2 - mean1)^2
  ccc <- 2 * cov12/denom
  numer <- 2 * x1dev * x2dev - ccc * (x1^2 + x2^2 - 2 * x1 * mean2 - 2 * 
                                        x2 * mean1 + 2 * mean1 * mean2)
  ccc.infl <- numer/denom
  
  # mean difference
  meandiff <- mean2 - mean1
  meandiff.infl <- (x2 - x1) - meandiff
  
  # variance ratio
  varratio <- var2/var1
  varratio.infl <- (var1 * (x2dev^2 - var2) - var2 * (x1dev^2 - var1))/var1^2
  
  est <- c(mean1 = mean1, mean2 = mean2, var1 = var1, var2 = var2, cov12 = cov12, 
           corr = cov12/sqrt(var1 * var2), sd1 = sqrt(var1), sd2 = sqrt(var2), 
           tdi = tdi, ccc = ccc, cp = cp, meandiff = meandiff, varratio = varratio)
  
  infl <- list(cp = cp.infl, ccc = ccc.infl, meandiff = meandiff.infl, varratio = varratio.infl)
  
  result <- list(est = est, infl = infl)
  return(result)
}

######################
#
# Function for computing the TDI agreement measure
# Arguments: 
#		paired.dat = a data frame containing all possible pairs
#		wt.paired = a vector of same length as paired.dat containing the weights
#					for each observation pair
#		prob = the probability
# Value:
#		the value of TDI (prob)		
#
tdi.fun.npm <- function(paired.dat, wt.paired, prob) {
  x1 <- paired.dat$x1
  x2 <- paired.dat$x2
  absd <- abs(x2 - x1)
  absd.pmf <- aggregate(wt.paired, by = list(absd = absd), FUN = sum)
  colnames(absd.pmf) <- c("value", "pmf")
  absd.cdf <- cumsum(absd.pmf$pmf)
  tdi <- min(absd.pmf[absd.cdf >= prob, "value"])
  return(tdi)
}

######################
#
# Function for computing the coverage probability (CP) agreement measure
# Arguments: 
#		paired.dat = a data frame containing all possible pairs
#		wt.paired = a vector of same length as paired.dat containing the weights
#					for each observation pair
#		delta = the acceptable margin
# Value:
#		the value of CP measure at delta		
#
cp.fun.npm <- function(paired.dat, wt.paired, delta) {
  x1 <- paired.dat$x1
  x2 <- paired.dat$x2
  absd <- abs(x2 - x1)
  cp <- sum(wt.paired * (absd <= delta))
  return(cp)
}



######################
#
# Function for forming all possible pairs of unlinked repeated measurements
# from one subject for a specified pair of methods
# Arguments: 
#		dat.one = a data frame containing data on one subject
#		nrep = common number of replications
#		method.pair = a character vector of length two specifying the two methods
# Value:
#		a data frame with three columns and nrep^2 rows

make.pairs.unlinked <- function(dat.one, nrep, method.pair) {
  subid <- dat.one$subject[1]
  x1 <- dat.one[dat.one$method == method.pair[1], "meas"]
  x2 <- dat.one[dat.one$method == method.pair[2], "meas"]
  paired.dat <- data.frame(rep(subid, nrep^2), rep(x1, each = nrep), rep(x2, 
                                                                         nrep))
  #	names(paired.dat) <- c("subject", method.pair[1], method.pair[2])
  names(paired.dat) <- c("subject", "x1", "x2")
  return(paired.dat)
}

####################################################################################
# chapter 9 #
####################################################################################

####################################
#
# Function for computing confidence intervals and individual confidence bands
# on a grid
# Arguments: 
#		param.hat = vector of estimated transformed model parameters
#		hess = Hessian matrix of parameter estimates
#		Xgrid = a matrix ngrid x number of fixed effects parameters. ngrid is the 
#			number of points on the grid on which the simultaneous confidence band 
#			is to be computed
#		X.colnames = a character vector containing names of the columns of the 
#			fixed effects design matrix
#		method.levels = a character vector containing names of the two methods
#		interaction = a logical indicating whether subject x method interaction
#			is included in the model
#		prob = probability cutoff for TDI
#		seed = seed needed for multcomp because it uses a randomized algorithm 
#			to compute critical points
#
# Note: For this function, it does not matter whether err.cor = T or F

conf.measures.longitudinal <- function(param.hat, hess, interaction, Xgrid = X.g, 
                                       method.levels = method.names, X.colnames = colnames(X), conf.level = 0.95, 
                                       prob = 0.9, seed = 12345) {
  require(numDeriv)
  require(multcomp)
  set.seed(seed)
  z.twosided <- qnorm((1 + conf.level)/2)
  
  # invert the Hessian matrix
  inv.hess <- solve(hess)
  rownames(inv.hess) <- names(param.hat)
  colnames(inv.hess) <- names(param.hat)
  
  # param - individual intervals for transformed model parameters
  se.param.hat <- sqrt(diag(inv.hess))
  ci.param <- cbind(param.hat - z.twosided * se.param.hat, param.hat + z.twosided * 
                      se.param.hat)
  ci.param <- cbind(param.hat, se.param.hat, ci.param)
  rownames(ci.param) <- names(param.hat)
  colnames(ci.param) <- c("estimate", "se", "lcl", "ucl")
  
  # lambda - individual interval for precision ratio
  nm <- paste0("lsigmasq.e", method.names)
  llambda.hat <- param.hat[nm[1]] - param.hat[nm[2]]
  names(llambda.hat) <- NULL
  se.llambda.hat <- sqrt(inv.hess[nm[1], nm[1]] + inv.hess[nm[2], nm[2]] - 
                           2 * inv.hess[nm[1], nm[2]])
  ci.llambda <- llambda.hat + c(-1, 1) * z.twosided * se.llambda.hat
  ci.llambda <- list(lresult = c(est = llambda.hat, se = se.llambda.hat, 
                                 lcl = ci.llambda[1], ucl = ci.llambda[2]), result = exp(c(est = llambda.hat, 
                                                                                           lcl = ci.llambda[1], ucl = ci.llambda[2])))
  # individual confidence bands on grid
  cmat <- diag(nrow(Xgrid))
  
  # mean function for method 1 -individual two-sided band on the grid
  ci.mean.1 <- ci.ind.paramfun(param.hat, inv.hess, mean.longitudinal.fun, 
                               cmat, "two.sided", conf.level, mean.flag = "1", method.levels = method.levels, 
                               Xgrid = Xgrid, X.colnames = X.colnames)
  
  # mean function for method 2 - individual two-sided band on the grid
  ci.mean.2 <- ci.ind.paramfun(param.hat, inv.hess, mean.longitudinal.fun, 
                               cmat, "two.sided", conf.level, mean.flag = "2", method.levels = method.levels, 
                               Xgrid = Xgrid, X.colnames = X.colnames)
  
  # mean function for difference (method 2 - method 1) - individual two-sided band on the grid
  ci.mean.D <- ci.ind.paramfun(param.hat, inv.hess, mean.longitudinal.fun, 
                               cmat, "two.sided", conf.level, mean.flag = "D", method.levels = method.levels, 
                               Xgrid = Xgrid, X.colnames = X.colnames)
  
  # z(ccc) - individual lower band on the grid
  ci.zccc <- ci.ind.paramfun(param.hat, inv.hess, zccc.longitudinal.fun, 
                             cmat, "greater", conf.level, Xgrid = Xgrid, X.colnames = X.colnames,
                             method.levels = method.levels, interaction = interaction)
  ci.zccc <- list(ci.fun = ci.zccc$ci.fun, zresult = ci.zccc$result, result = tanh(ci.zccc$result[, 
                                                                                                  c("estimate", "lcl", "ucl"), drop = FALSE]))
  
  # log(tdi) - individual upper band on the grid
  ci.ltdi <- ci.ind.paramfun(param.hat, inv.hess, ltdi.longitudinal.fun, 
                             cmat, "less", conf.level, prob = prob, Xgrid = Xgrid, X.colnames = X.colnames,
                             method.levels = method.levels, interaction = interaction)
  ci.ltdi <- list(ci.fun = ci.ltdi$ci.fun, lresult = ci.ltdi$result, result = exp(ci.ltdi$result[, 
                                                                                                 c("estimate", "lcl", "ucl"), drop = FALSE]))
  
  results <- list(param = ci.param, mean.1 = ci.mean.1, mean.2 = ci.mean.2, 
                  mean.D = ci.mean.D, precision.ratio = ci.llambda, tdi = ci.ltdi, ccc = ci.zccc)
  return(results)
}
##################
#
# Function for computing log (TDI) as a function of model parameters
# Arguments: 
#		param = transformed model parameters
#		Xgrid =
#		X.colnames =
#		method.levels = 
#		interaction = 
#		prob = 
# (these agruments are same as those for conf.measures.longitudinal)
#		ulim = upper limit for the interval neeed for the root finding algorithm
# Note: For this function, it does not matter whether err.cor = T or F

ltdi.longitudinal.fun <- function(param, prob, Xgrid, X.colnames, method.levels, 
                                  interaction) {
  # get the moments
  sigmasq.b <- exp(param["lsigmasq.b"])
  sigmasq.bstar <- exp(param["lsigmasq.bstar"])
  sigmasq.e <- exp(param[paste0("lsigmasq.e", method.levels)])
  cov.12 <- sigmasq.b + sigmasq.bstar
  if (interaction == F) {
    var.1 <- sigmasq.b + sigmasq.bstar + sigmasq.e[1]
    var.2 <- sigmasq.b + sigmasq.bstar + sigmasq.e[2]
  }
  if (interaction == T) {
    psisq <- exp(param["lpsisq"])
    var.1 <- sigmasq.b + psisq + sigmasq.bstar + sigmasq.e[1]
    var.2 <- sigmasq.b + psisq + sigmasq.bstar + sigmasq.e[2]
  }
  mean.D <- mean.longitudinal.fun(param, Xgrid = Xgrid, X.colnames = X.colnames, 
                                  method.levels = method.levels, mean.flag = "D")
  sd.D <- sqrt(var.1 + var.2 - 2 * cov.12)
  result <- mapply(ltdi, dmean = mean.D, MoreArgs = list(dsd = sd.D, prob = prob))
  names(result) <- NULL
  return(result)
}

##################
#
# Function for computing Fisher's z-transformation of CCC as a function of model parameters
# Arguments: 
#		param = transformed model parameters
#		Xgrid =
#		X.colnames =
#		method.levels = 
#		interaction = 
# (these agruments are same as those for conf.measures.longitudinal)
# Note: For this function, it does not matter whether err.cor = T or F
#
zccc.longitudinal.fun <- function(param, Xgrid, X.colnames, method.levels, 
                                  interaction) {
  # get the moments
  sigmasq.b <- exp(param["lsigmasq.b"])
  sigmasq.bstar <- exp(param["lsigmasq.bstar"])
  sigmasq.e <- exp(param[paste0("lsigmasq.e", method.levels)])
  cov.12 <- sigmasq.b + sigmasq.bstar
  mean.D <- mean.longitudinal.fun(param, Xgrid = Xgrid, X.colnames = X.colnames, 
                                  method.levels = method.levels, mean.flag = "D")
  if (interaction == F) {
    var.1 <- sigmasq.b + sigmasq.bstar + sigmasq.e[1]
    var.2 <- sigmasq.b + sigmasq.bstar + sigmasq.e[2]
  }
  if (interaction == T) {
    psisq <- exp(param["lpsisq"])
    var.1 <- sigmasq.b + psisq + sigmasq.bstar + sigmasq.e[1]
    var.2 <- sigmasq.b + psisq + sigmasq.bstar + sigmasq.e[2]
  }
  result <- mapply(zccc, dmean = mean.D, MoreArgs = list(var1 = var.1, var2 = var.2, 
                                                         cov12 = cov.12))
  names(result) <- NULL
  return(result)
}

##################
#
# Function for computing mean of first method as a function of model parameters
# Arguments: 
#		param = transformed model parameters
#		Xgrid =
#		X.colnames =
#		method.levels = 
# (these agruments are same as those for conf.measures.longitudinal)
#		mean.flag = possible values "1" (for method 1), "2" (for method 2) 
#		or "D" (for method 2 - method 1)
# Note: For this function, it does not matter whether err.cor = T or F
#
mean.longitudinal.fun <- function(param, Xgrid, X.colnames, method.levels, 
                                  mean.flag) {
  beta <- param[X.colnames]
  if (mean.flag %in% c("1", "2", "D") == F) {
    stop("not a valid option for mean.flag in mean.longitudinal.fun")
  }
  if (mean.flag == "1") {
    result <- Xgrid %*% beta[c(paste0("method", method.levels[1]), paste0("method", 
                                                                          method.levels[1], ":time"), paste0("method", method.levels[1], 
                                                                                                             ":I(time^2)"), paste0("method", method.levels[1], ":I(time^3)"))]
  }
  if (mean.flag == "2") {
    result <- Xgrid %*% beta[c(paste0("method", method.levels[2]), paste0("method", 
                                                                          method.levels[2], ":time"), paste0("method", method.levels[2], 
                                                                                                             ":I(time^2)"), paste0("method", method.levels[2], ":I(time^3)"))]
  }
  if (mean.flag == "D") {
    result1 <- Xgrid %*% beta[c(paste0("method", method.levels[1]), paste0("method", 
                                                                           method.levels[1], ":time"), paste0("method", method.levels[1], 
                                                                                                              ":I(time^2)"), paste0("method", method.levels[1], ":I(time^3)"))]
    result2 <- Xgrid %*% beta[c(paste0("method", method.levels[2]), paste0("method", 
                                                                           method.levels[2], ":time"), paste0("method", method.levels[2], 
                                                                                                              ":I(time^2)"), paste0("method", method.levels[2], ":I(time^3)"))]
    result <- result2 - result1
  }
  return(as.numeric(result))
}


############################

inv.logit <- function(x) {
  1/(1 + exp(-x))
}

###################
#
# Function for log-likelihood function for the fitted model (with 
#	or without subject x methd interaction)
# to bivariate longitudinal data
# Arguments: 
#		param = vector of ML estimates of transformed parameters
#		dat = a data frame with five columns --- subject (a factor), 
#			  method (a factor with two levels), 
#			rep (a factor representing measurement occasion number),
#			time (the time of measurement) and meas (the observation)
#		Xmat = fixed-effects design matrix
#		Zmat = design matrix for random effects
# 		method.levels = a character vector containing levels of the method factor
#		rep.levels = a character vector containing levels of the rep factor
#		err.cor = logical indicating whether a continuous-AR(1) correlation
#			structure is assumed for within-subject errors
# Value: 
#		log-likelihood for the fitted model
# Note: this function is very slow, but it works!

loglik.fun <- function(param, dat, Xmat = X, Zmat = Z.noint, method.levels = method.names, 
                       rep.levels = rep.names, err.cor = T, interaction = F) {
  require(mvtnorm)
  require(Matrix)
  require(mvtnorm)
  require(multcomp)
  
  subid <- unique(dat$subject)
  nsub <- length(subid)
  nobs <- nrow(dat)
  
  # extract parameters
  beta <- param[colnames(Xmat)]
  sigmasq.b <- exp(param["lsigmasq.b"])
  sigmasq.bstar <- exp(param["lsigmasq.bstar"])
  sigmasq.e <- exp(param[paste0("lsigmasq.e", method.levels)])
  
  if (err.cor == T) {
    phi <- 1/(1 + exp(-param["logitphi"]))
  }
  
  # covariance matrices of random vectors
  if (interaction == F) {
    var.re.diag <- rep(c(sigmasq.b, rep(sigmasq.bstar, length(rep.levels))), 
                       nsub)
  }
  
  if (interaction == T) {
    psisq <- exp(param["lpsisq"])
    var.re.diag <- rep(c(sigmasq.b, rep(psisq, 2), rep(sigmasq.bstar, length(rep.levels))), 
                       nsub)
  }
  
  var.re <- diag(var.re.diag, nrow = length(var.re.diag))
  
  # covariance matrix of errors
  if (err.cor == F) {
    var.e.diag <- Xmat[, paste0("method", method.levels[1])] * sigmasq.e[1] + 
      Xmat[, paste0("method", method.levels[2])] * sigmasq.e[2]
    var.e <- diag(var.e.diag, nrow = length(var.e.diag))
  }
  
  if (err.cor == T) {
    idx <- with(dat, interaction(subject, method, drop = T, lex.order = T))
    var.e <- lapply(split(dat, idx, drop = T), function(df) {
      vv1.diag <- (df$method == method.levels[1]) * sigmasq.e[1] + (df$method == 
                                                                      method.levels[2]) * sigmasq.e[2]
      vv1 <- diag(vv1.diag, nrow = length(vv1.diag))
      vv2 <- phi^abs(outer(df$time, df$time, "-"))
      return(vv1 %*% vv2)
    })
    var.e <- as.matrix(bdiag(var.e))
  }
  
  # components in the covariance matrix of y
  var.re <- Zmat %*% var.re %*% t(Zmat)
  
  mean.y <- as.numeric(Xmat %*% beta)
  var.y <- var.re + var.e
  
  loglik <- dmvnorm(dat$meas, mean = mean.y, sigma = var.y, log = T)
  return(loglik)
}


################################################################
#
# Function for extracting estimates of parameters of model and 
# computing estimates of various parameter functions
#		fm = fitted model
#		interaction = a logical indicating whether subject x method interaction
#			is included in the model
# 		method.levels = a character vector containing levels of the method factor
#		rep.levels = a character vector containing levels of the rep factor
#		Xgrid = a matrix ngrid x number of fixed effects parameters. ngrid is the 
#			number of points on the grid on which the simultaneous confidence band 
#			is to be computed
#		X.colnames = a character vector containing names of the columns of the 
#			fixed effects design matrix
#		prob = probability cutoff for TDI
#		err.cor = logical indicating whether a continuous-AR(1) correlation
#			structure is assumed for within-subject errors

get.param.estimates <- function(fm, interaction, method.levels, rep.levels, 
                                Xgrid, X.colnames, prob, err.cor = T) {
  # fixed effects
  beta.hat <- fixef(fm)
  names(beta.hat) <- X.colnames
  
  # variances of random effects
  sigma.hat <- fm$sigma
  mat1 <- diag(as.matrix(fm$modelStruct$reStruct[[1]]) * sigma.hat^2)
  sigmasq.b.hat <- mat1["(Intercept)"]
  sigmasq.bstar.hat <- mat1[paste0("rep", rep.levels[1])]
  
  # error variances
  sigma.ratio.hat <- coef(fm$modelStruct$varStruct, unconstrained = FALSE, 
                          allCoef = TRUE)[method.levels]
  sigmasq.e.hat <- (sigma.hat * sigma.ratio.hat)^2
  names(sigmasq.e.hat) <- method.levels
  
  if (interaction == F) {
    var.1.hat <- sigmasq.b.hat + sigmasq.bstar.hat + sigmasq.e.hat[1]
    var.2.hat <- sigmasq.b.hat + sigmasq.bstar.hat + sigmasq.e.hat[2]
    # transformed parameter estimates
    if (err.cor == F) {
      param.hat <- c(beta.hat, log(sigmasq.b.hat), log(sigmasq.bstar.hat), 
                     log(sigmasq.e.hat))
      names(param.hat) <- c(X.colnames, "lsigmasq.b", "lsigmasq.bstar", 
                            paste0("lsigmasq.e", method.levels))
    }
    if (err.cor == T) {
      phi.hat <- coef(fm$modelStruct$corStruct, unconstrained = FALSE)
      param.hat <- c(beta.hat, log(sigmasq.b.hat), log(sigmasq.bstar.hat), 
                     log(sigmasq.e.hat), log(phi.hat/(1 - phi.hat)))
      names(param.hat) <- c(X.colnames, "lsigmasq.b", "lsigmasq.bstar", 
                            paste0("lsigmasq.e", method.levels), "logitphi")
    }
  }
  
  if (interaction == T) {
    psisq.hat <- mat1[paste0("method", method.levels[1])]
    var.1.hat <- sigmasq.b.hat + psisq.hat + sigmasq.bstar.hat + sigmasq.e.hat[1]
    var.2.hat <- sigmasq.b.hat + psisq.hat + sigmasq.bstar.hat + sigmasq.e.hat[2]
    # transformed parameter estimates
    if (err.cor == F) {
      param.hat <- c(beta.hat, log(sigmasq.b.hat), log(psisq.hat), log(sigmasq.bstar.hat), 
                     log(sigmasq.e.hat))
      names(param.hat) <- c(X.colnames, "lsigmasq.b", "lpsisq", "lsigmasq.bstar", 
                            paste0("lsigmasq.e", method.levels))
    }
    if (err.cor == T) {
      phi.hat <- coef(fm$modelStruct$corStruct, unconstrained = FALSE)
      param.hat <- c(beta.hat, log(sigmasq.b.hat), log(psisq.hat), log(sigmasq.bstar.hat), 
                     log(sigmasq.e.hat), log(phi.hat/(1 - phi.hat)))
      names(param.hat) <- c(X.colnames, "lsigmasq.b", "lpsisq", "lsigmasq.bstar", 
                            paste0("lsigmasq.e", method.levels), "logitphi")
    }
  }
  
  cov.12.hat <- sigmasq.b.hat + sigmasq.bstar.hat
  cor.12.hat <- cov.12.hat/sqrt(var.1.hat * var.2.hat)
  names(var.1.hat) <- NULL
  names(var.2.hat) <- NULL
  names(cov.12.hat) <- NULL
  names(cor.12.hat) <- NULL
  
  # log(lambda)
  llambda.hat <- log(sigmasq.e.hat[1]/sigmasq.e.hat[2])
  names(llambda.hat) <- NULL
  
  # parameter values on grid
  
  # means
  mu.1.hat <- as.numeric(Xgrid %*% beta.hat[c(paste0("method", method.levels[1]), 
                                              paste0("method", method.levels[1], ":time"), paste0("method", method.levels[1], 
                                                                                                  ":I(time^2)"), paste0("method", method.levels[1], ":I(time^3)"))])
  
  mu.2.hat <- as.numeric(Xgrid %*% beta.hat[c(paste0("method", method.levels[2]), 
                                              paste0("method", method.levels[2], ":time"), paste0("method", method.levels[2], 
                                                                                                  ":I(time^2)"), paste0("method", method.levels[2], ":I(time^3)"))])
  
  # mean and SD of difference
  mu.D.hat <- mu.2.hat - mu.1.hat
  sd.D.hat <- sqrt(var.1.hat + var.2.hat - 2 * cov.12.hat)
  
  # LOA 
  loa <- matrix(0, nrow = nrow(Xgrid), ncol = 2)
  colnames(loa) <- c("lower", "upper")
  loa[, 1] <- mu.D.hat - 1.96 * sd.D.hat
  loa[, 2] <- mu.D.hat + 1.96 * sd.D.hat
  
  # CCC
  zccc.hat <- mapply(zccc, dmean = mu.D.hat, MoreArgs = list(var1 = var.1.hat, 
                                                             var2 = var.2.hat, cov12 = cov.12.hat))
  
  # TDI
  ltdi.hat <- mapply(ltdi, dmean = mu.D.hat, MoreArgs = list(dsd = sd.D.hat, 
                                                             prob = prob))
  
  return(list(param = param.hat, bvn = list(mean.1 = mu.1.hat, mean.2 = mu.2.hat, 
                                            sd.1 = sqrt(var.1.hat), sd.2 = sqrt(var.2.hat), cor = cor.12.hat), 
              mean.D = mu.D.hat, sd.D = sd.D.hat, llambda = llambda.hat, loa = loa, 
              zccc = zccc.hat, ltdi = ltdi.hat))
}


###############################################################
# Function for estimating semvariograms of errors of the two methods

variogram.resid <- function(dat, res, robust = F) {
  dat.res <- data.frame(dat, res = res)
  method.names <- levels(dat$method)
  dat.1 <- subset(dat.res, method == method.names[1])
  result.1 <- new.variogram.lme(dat.1, robust = robust)
  dat.2 <- subset(dat.res, method == method.names[2])
  result.2 <- new.variogram.lme(dat.2, robust = robust)
  result <- list(result.1, result.2)
  names(result) <- method.names
  return(result)
}

# The new.varigraom.lme function is a just a modification of Variogram.lme function,  
# obtained as getAnywhere("Variogram.lme"), which accepts a data frame with residuals
# and time instead of extracting them from the fitted model

# dat <- fat
# res <- resid(fit, type = "n")
# dat.res <- data.frame(dat, res = res)
# df <- dat.res
# new.variogram.lme (df) and Variogram(fit, form = ~time, resType = "n") provide
# identical answers

new.variogram.lme <- function(df, maxDist, length.out = 50, collapse = c("quantiles", 
                                                                         "fixed", "none"), nint = 20, breaks, robust = F) {
  grps <- df$subject
  covar <- split(df$time, grps)
  covar <- covar[sapply(covar, function(el) nrow(as.matrix(el))) > 1]
  distance <- lapply(covar, function(el) dist(as.matrix(el)))
  res <- split(df$res, grps)
  res <- res[sapply(res, length) > 1]
  levGrps <- levels(grps)
  val <- vector("list", length(levGrps))
  names(val) <- levGrps
  for (i in levGrps) {
    val[[i]] <- Variogram(res[[i]], distance[[i]])
  }
  val <- do.call("rbind", val)
  if (!missing(maxDist)) {
    val <- val[val$dist <= maxDist, ]
  }
  collapse <- match.arg(collapse)
  if (collapse != "none") {
    dst <- val$dist
    udist <- sort(unique(dst))
    ludist <- length(udist)
    if (!missing(breaks)) {
      if (min(breaks) > udist[1L]) {
        breaks <- c(udist[1L], breaks)
      }
      if (max(breaks) < udist[2L]) {
        breaks <- c(breaks, udist[2L])
      }
      if (!missing(nint) && nint != (length(breaks) - 1L)) {
        stop("'nint' is not consistent with 'breaks'")
      }
      nint <- length(breaks) - 1
    }
    if (nint < ludist) {
      if (missing(breaks)) {
        if (collapse == "quantiles") {
          breaks <- unique(quantile(dst, seq(0, 1, 1/nint)))
        } else {
          breaks <- seq(udist[1L], udist[length(udist)], length = nint + 
                          1L)
        }
      }
      cutDist <- cut(dst, breaks)
    } else {
      cutDist <- dst
    }
    val <- lapply(split(val, cutDist), function(el, robust) {
      nh <- nrow(el)
      vrg <- el$variog
      
      if (robust) {
        vrg <- ((mean(vrg^0.25))^4)/(0.457 + 0.494/nh)
      } else {
        vrg <- mean(vrg)
      }
      dst <- median(el$dist)
      data.frame(variog = vrg, dist = dst)
    }, robust = robust)
    
    val <- do.call("rbind", val)
    val$n.pairs <- as.vector(table(na.omit(cutDist)))
    val <- na.omit(val)
  }
  row.names(val) <- 1:nrow(val)
  return(val)
}


###############################################################

#
# Function for ACF of residuals 
# Arguments: 
#		dat = a data frame with five columns --- subject (a factor), 
#			  method (a factor with two levels), 
#			rep (a factor representing measurement occasion number),
#			time (the time of measurement) and meas (the observation);
#			time and meas are not used
#		res = residuals from the mixed-model fit to dat 
#		lag.max = maximum number of lags
# 		cor.ind = logical indicating autocorrelation or autocovariance is desired
# 		alpha = level of test of null hypothesis of independence
# Value:
#		A list of length 2 contanining acf, etc.

acf.resid <- function(dat, res, lag.max, cor.ind = T, alpha = 0.05) {
  
  #	function for computing sum of products for a given subject
  getSum <- function(df, lag, rep.max) {
    ss <- 0
    nn <- 0
    if (lag <= (rep.max - 1)) {
      res <- df[, "res"]
      for (k in 1:rep.max) {
        if ((k %in% df$rep) && (k + lag) %in% df$rep) {
          idx1 <- match(k, df$rep)
          idx2 <- match((k + lag), df$rep)
          ss <- ss + res[idx1] * res[idx2]
          nn <- nn + 1
        }
      }
    }
    return(c(nn, ss))
  }
  
  # function for computing autocovariance at a fixed lag
  acov.fun <- function(dat, lag, rep.max) {
    ss <- Reduce("+", lapply(split(dat, dat$subject, drop = T), getSum, 
                             lag = lag, rep.max = rep.max))
    return(c(ss[1], ss[2]/ss[1]))
  }
  
  # function for computing acf at all lags from 0 to lag.max
  acf.lags <- function(dat, lag.max, rep.max, cor.ind, alpha) {
    lags <- 0:lag.max
    covs <- data.frame(lags, t(sapply(lags, acov.fun, dat = dat, rep.max = rep.max)))
    colnames(covs) <- c("lag", "n.used", "acf")
    cval <- qnorm(1 - alpha/2)/sqrt(covs[, "n.used"])
    cval[1] <- NA
    if (cor.ind == T) {
      acf.val <- covs[, "acf"]/covs[1, "acf"]
      covs <- data.frame(lag = covs[, "lag"], n.used = covs[, "n.used"], 
                         acf = acf.val, cval = cval, accept.null = (abs(acf.val) < cval))
    }
    return(covs)
  }
  # do the computation
  
  dat.res <- data.frame(dat, res = res)
  method.names <- levels(dat$method)
  rep.max <- nlevels(dat$rep)
  dat.1 <- subset(dat.res, method == method.names[1])
  dat.1 <- with(dat.1, dat.1[order(subject, rep), ])
  result.1 <- acf.lags(dat.1, lag.max, rep.max, cor.ind, alpha)
  dat.2 <- subset(dat.res, method == method.names[2])
  dat.2 <- with(dat.2, dat.2[order(subject, rep), ])
  result.2 <- acf.lags(dat.2, lag.max, rep.max, cor.ind, alpha)
  result <- list(result.1, result.2)
  names(result) <- method.names
  return(result)
}

####################################################################################
# chapter 8 #
####################################################################################
#
# Function for fitting a heteroscedastic mixed model to unlinked repeated 
# measurements data without method x subject interaction
# Arguments: 
#		dat = a data frame with four columns --- subject (a factor), 
#			  method (a factor with two levels method1 and method2), 
#			rep (replication number, unused) and meas
#		vtilde = a vector of length nrow(dat) containing the subject-specific 
#				value of the variance covariate 
#		varfunc = a character with one of two possible values---
#				"power" or "exponential"--- indicating the common variance function
#				model for the two methods
# Value: The fitted object from the lme fit
#

lme.noint.hetero.fit <- function(dat, vtilde, varfunc) {
  require(nlme)
  
  if (!(varfunc %in% c("power", "exponential"))) {
    stop(paste(varfunc, "is not a valid option in fit.hetero.lme"))
  }
  
  # fit model
  dat <- data.frame(dat, vtilde)
  dat <- groupedData(meas ~ method | subject, data = dat, order.groups = FALSE)
  
  if (varfunc == "power") {
    fm.power <- lme(meas ~ method - 1, data = dat, random = ~1 | subject, method = "ML", 
                    weights = varComb(varIdent(form = ~1 | method), varPower(form = ~vtilde | 
                                                                               method)), control = lmeControl(maxIter = 1000, msMaxIter = 5000, niterEM = 2000, 
                                                                                                              msMaxEval = 1500))
    return(fm.power)
  }
  
  if (varfunc == "exponential") {
    fm.exp <- lme(meas ~ method - 1, data = dat, random = ~1 | subject, method = "ML", 
                  weights = varComb(varIdent(form = ~1 | method), varExp(form = ~vtilde | method)), 
                  control = lmeControl(maxIter = 1000, msMaxIter = 5000, niterEM = 2000, msMaxEval = 1500))
    return(fm.exp)
  }
  
}

#########################################################################

#
# Function for computing estimates, their SEs and confidence intervals
#		 for measures of similarity, agreement and repeatability
# Arguments: 
#		param.hat = vector of ML estimates of transformed parameters
#		vcov.param.hat = covariance matrix (inverse Hessian matrix) for param.hat
#		method.names = a character vector of length two containing names of methods 
#		vtilde = a vector of length nrow(dat) containing the subject-specific 
#				value of the variance covariate 
#		varfunc = a character with one of two possible values---
#				"power" or "exponential"--- indicating the common variance function
#				model for the two methods
#		vgrid = a grid of values on which to evaluate the fitted error variance 
#					functions and the paramerters that involve them
#		mean.covariates = logical indicate whether mean.covariates are included in the model
#		contr.mat = if mean.covariates are included, then contr.mat provides the matrix of 
#			contrast coefficients which when premultiplied to fixed(fit) provides the
#			vector of difference in means of two methods at each covariate setting
#		conf.level = level of confidence (sometimes individual, sometimes simultaneous)
#		prob = probability cutoff used in TDI's definition
#		seed = seed to fix because multcomp uses a randomized algorithm to compute critical point
# Value: 
#		lists containing a bunch of stuff, including
#		estimate, se and confidence intervals for measures of similarity
#		estimate, se and ucb of log(TDI), and estimate and ucb of TDI
#		estimate, se and lcb of z(CCC), and estimate and lcb of CCC
#		estimates, se and ucb of log(TDI.1) and log(TDI.2), and estimates and ucb for 
#				TDI.1 and TDI.2
#		estimates, se and lcb for z(CCC.1) and z(CCC.2), and estimates and lcb for 
#				CCC.1 and CCC.2
#	Note: some of the intervals and bounds are individual whereas others are simultaneous - check the code


conf.measures.noint.hetero.reg <- function(param.hat, vcov.param.hat, method.names, varfunc, 
                                           vgrid, mean.covariates, contr.mat, conf.level = 0.95, prob = 0.9, seed = 54321) {
  require(multcomp)
  set.seed(seed)
  q.chisq <- qchisq(prob, 1)
  inv.hess <- vcov.param.hat
  
  # evaluation of similarity
  # confidence intervals for mean difference
  if (mean.covariates == F) {
    contr.mat <- diag(1, nrow = 1)
  }
  ci.mean <- ci.ind.paramfun(param.hat, inv.hess, mean.diff.fun.reg, diag(1, nrow = nrow(contr.mat)), 
                             "two.sided", conf.level, method.names = method.names, mean.covariates = mean.covariates, 
                             contr.mat = contr.mat)
  
  # confidence intervals for precision ratio
  contrast.matrix <- diag(1, nrow = length(vgrid))
  ci.lprec <- ci.simul.paramfun(param.hat, inv.hess, lprec.ratio.fun.reg, contrast.matrix, 
                                "two.sided", conf.level, varfunc = varfunc, vgrid = vgrid, method.names = method.names)
  ci.lprec <- list(ci.fun = ci.lprec$ci.fun, lresult = ci.lprec$result, result = exp(ci.lprec$result[, 
                                                                                                     c("estimate", "lcl", "ucl"), drop = FALSE]))
  
  # evaluation of agreement 
  # confidence intervals for tdi for the given covariate.level
  ci.ltdi <- vector("list", nrow(contr.mat))
  for (i in 1:nrow(contr.mat)) {
    rr <- ci.simul.paramfun(param.hat, inv.hess, ltdi.fun.reg, contrast.matrix, "less", 
                            conf.level, varfunc = varfunc, vgrid = vgrid, method.names = method.names, 
                            mean.covariates = mean.covariates, contr.mat = contr.mat, covariate.level = i, 
                            prob = prob)
    ci.ltdi[[i]] <- list(ci.fun = rr$ci.fun, lresult = rr$result, result = exp(rr$result[, 
                                                                                         c("estimate", "lcl", "ucl"), drop = FALSE]))
  }
  
  # confidence intervals for ccc
  ci.zccc <- vector("list", nrow(contr.mat))
  for (i in 1:nrow(contr.mat)) {
    rr <- ci.simul.paramfun(param.hat, inv.hess, zccc.fun.reg, contrast.matrix, "greater", 
                            conf.level, varfunc = varfunc, vgrid = vgrid, method.names = method.names, 
                            mean.covariates = mean.covariates, contr.mat = contr.mat, covariate.level = i)
    ci.zccc[[i]] <- list(ci.fun = rr$ci.fun, zresult = rr$result, result = tanh(rr$result[, 
                                                                                          c("estimate", "lcl", "ucl"), drop = FALSE]))
  }
  
  # evaluation of repetability
  # confidence intervals for tdi for method 1
  rr <- ci.simul.paramfun(param.hat, inv.hess, ltdi.1.fun.reg, contrast.matrix, "less", 
                          conf.level, varfunc = varfunc, vgrid = vgrid, method.names = method.names, q.chisq = q.chisq)
  ci.ltdi.1 <- list(ci.fun = rr$ci.fun, lresult = rr$result, result = exp(rr$result[, 
                                                                                    c("estimate", "lcl", "ucl"), drop = FALSE]))
  
  # confidence intervals for tdi for method 2
  rr <- ci.simul.paramfun(param.hat, inv.hess, ltdi.2.fun.reg, contrast.matrix, "less", 
                          conf.level, varfunc = varfunc, vgrid = vgrid, method.names = method.names, q.chisq = q.chisq)
  ci.ltdi.2 <- list(ci.fun = rr$ci.fun, lresult = rr$result, result = exp(rr$result[, 
                                                                                    c("estimate", "lcl", "ucl"), drop = FALSE]))
  
  ci.rep.ltdi <- list(tdi.1 = ci.ltdi.1, tdi.2 = ci.ltdi.2)
  
  # confidence intervals for ccc for method 1
  rr <- ci.simul.paramfun(param.hat, inv.hess, zccc.1.fun.reg, contrast.matrix, "greater", 
                          conf.level, varfunc = varfunc, vgrid = vgrid, method.names = method.names)
  ci.zccc.1 <- list(ci.fun = rr$ci.fun, zresult = rr$result, result = tanh(rr$result[, 
                                                                                     c("estimate", "lcl", "ucl"), drop = FALSE]))
  
  # confidence intervals for ccc for method 2
  rr <- ci.simul.paramfun(param.hat, inv.hess, zccc.2.fun.reg, contrast.matrix, "greater", 
                          conf.level, varfunc = varfunc, vgrid = vgrid, method.names = method.names)
  ci.zccc.2 <- list(ci.fun = rr$ci.fun, zresult = rr$result, result = tanh(rr$result[, 
                                                                                     c("estimate", "lcl", "ucl"), drop = FALSE]))
  
  ci.rep.zccc <- list(ccc.1 = ci.zccc.1, ccc.2 = ci.zccc.2)
  
  
  return(list(mean.diff = ci.mean, precision.ratio = ci.lprec, tdi = ci.ltdi, ccc = ci.zccc, 
              rep.tdi = ci.rep.ltdi, rep.ccc = ci.rep.zccc))
  
}

#########################################################################

#
# Function for computing z(CCC.1) as a function of model parameters on vgrid
# Arguments: 
#		param = transformed model parameters
#		varfunc = a character with one of two possible values---
#				"power" or "exponential"--- indicating the common variance function
#				model for the two methods
#		vgrid = a grid of values on which to evaluate the fitted error variance 
#					functions and the paramerters that involve them
#		method.names = a character vector of length two containing names of methods 
#	Value:
#		z(CCC.1) on vgrid
#


zccc.1.fun.reg <- function(param, varfunc, vgrid, method.names) {
  # extract parameters
  sigmasq.b <- exp(param["lsigmasq.b"])
  sigmasq.e <- exp(param[paste0("lsigmasq.e", method.names[1])])
  delta <- param[paste0("delta", method.names[1])]
  
  # get error variances on vgrid
  if (varfunc == "power") {
    err.var <- sigmasq.e * (abs(vgrid)^(2 * delta))
  }
  if (varfunc == "exponential") {
    err.var <- sigmasq.e * exp(2 * delta * vgrid)
  }
  result <- atanh(sigmasq.b/(sigmasq.b + err.var))
  return(result)
}
#########################################################################
#
# Function for computing z(CCC.2) as a function of model parameters on vgrid
# Arguments: 
#		param = transformed model parameters
#		varfunc = a character with one of two possible values---
#				"power" or "exponential"--- indicating the common variance function
#				model for the two methods
#		vgrid = a grid of values on which to evaluate the fitted error variance 
#					functions and the paramerters that involve them
#		method.names = a character vector of length two containing names of methods 
#	Value:
#		z(CCC.2) on vgrid
#

zccc.2.fun.reg <- function(param, varfunc, vgrid, method.names) {
  # extract parameters
  sigmasq.b <- exp(param["lsigmasq.b"])
  sigmasq.e <- exp(param[paste0("lsigmasq.e", method.names[2])])
  delta <- param[paste0("delta", method.names[2])]
  
  # get error variances on vgrid
  if (varfunc == "power") {
    err.var <- sigmasq.e * (abs(vgrid)^(2 * delta))
  }
  if (varfunc == "exponential") {
    err.var <- sigmasq.e * exp(2 * delta * vgrid)
  }
  result <- atanh(sigmasq.b/(sigmasq.b + err.var))
  return(result)
}

#########################################################################

#
# Function for computing z(CCC) on vgrid as a function of mean covariates
# Arguments: 
#		param = vector of transformed model parameters
#		varfunc = a character with one of two possible values---
#				"power" or "exponential"--- indicating the common variance function
#				model for the two methods
#		vgrid = a grid of values on which to evaluate the fitted error variance 
#					functions and the paramerters that involve them
#		method.names = a character vector of length two containing names of methods 
#		mean.covariates = logical indicate whether mean.covariates are included in the model
#		contr.mat = if mean.covariates are included, then contr.mat provides the matrix of 
#			contrast coefficients which when premultiplied to fixed(fit) provides the
#			vector of difference in means of two methods at each covariate setting
#		covariate.level = an integer indicating the given covariate level
# Value: 
#		z(CCC) on vgrid for the given covariate level


zccc.fun.reg <- function(param, varfunc, vgrid, method.names, mean.covariates, contr.mat, 
                         covariate.level) {
  # extract parameters
  sigmasq.b <- exp(param["lsigmasq.b"])
  names(sigmasq.b) <- "sigmasq.b"
  sigmasq.e <- exp(param[paste0("lsigmasq.e", method.names)])
  names(sigmasq.e) <- method.names
  delta <- param[paste0("delta", method.names)]
  names(delta) <- method.names
  
  # get error variances on vgrid
  err.var <- matrix(NA, nrow = length(vgrid), ncol = 2)
  colnames(err.var) <- method.names
  if (varfunc == "power") {
    err.var[, method.names[1]] <- sigmasq.e[method.names[1]] * (abs(vgrid)^(2 * delta[method.names[1]]))
    err.var[, method.names[2]] <- sigmasq.e[method.names[2]] * (abs(vgrid)^(2 * delta[method.names[2]]))
  }
  
  if (varfunc == "exponential") {
    err.var[, method.names[1]] <- sigmasq.e[method.names[1]] * exp(2 * delta[method.names[1]] * 
                                                                     vgrid)
    err.var[, method.names[2]] <- sigmasq.e[method.names[2]] * exp(2 * delta[method.names[2]] * 
                                                                     vgrid)
  }
  
  # get mean of D 
  Dmean <- mean.diff.fun.reg(param, method.names, mean.covariates, contr.mat) # may be a vector
  names(Dmean) <- NULL
  
  # get var and cov of Y
  Yvar <- sigmasq.b + err.var
  colnames(Yvar) <- method.names
  Ycov <- sigmasq.b
  
  # get z(CCC) on vgrid for Dmean at the given covariate level and Dsd
  result <- mapply(zccc, Dmean[covariate.level], Yvar[, method.names[1]], Yvar[, method.names[2]], 
                   Ycov)
  return(result)
}

#########################################################################
#
# Function for computing log(TDI.1) as a function of model parameters on vgrid
# Arguments: 
#		param = transformed model parameters
#		varfunc = a character with one of two possible values---
#				"power" or "exponential"--- indicating the common variance function
#				model for the two methods
#		vgrid = a grid of values on which to evaluate the fitted error variance 
#					functions and the paramerters that involve them
#		method.names = a character vector of length two containing names of methods 
#		q.chisq = chi-square percentile used in the definition of TDI (instead of prob)
#	Value:
#		log(TDI.1) on vgrid
#

ltdi.1.fun.reg <- function(param, varfunc, vgrid, method.names, q.chisq) {
  # extract parameters
  sigmasq.e <- exp(param[paste0("lsigmasq.e", method.names[1])])
  delta <- param[paste0("delta", method.names[1])]
  
  # get error variances on vgrid
  if (varfunc == "power") {
    err.var <- sigmasq.e * (abs(vgrid)^(2 * delta))
  }
  if (varfunc == "exponential") {
    err.var <- sigmasq.e * exp(2 * delta * vgrid)
  }
  result <- log(sqrt(2 * err.var * q.chisq))
  return(result)
}

#########################################################################
#
# Function for computing log(TDI.2) as a function of model parameters on vgrid
# Arguments: 
#		param = transformed model parameters
#		varfunc = a character with one of two possible values---
#				"power" or "exponential"--- indicating the common variance function
#				model for the two methods
#		vgrid = a grid of values on which to evaluate the fitted error variance 
#					functions and the paramerters that involve them
#		method.names = a character vector of length two containing names of methods 
#		q.chisq = chi-square percentile used in the definition of TDI  (instead of prob)
#	Value:
#		log(TDI.2) on vgrid
#

ltdi.2.fun.reg <- function(param, varfunc, vgrid, method.names, q.chisq) {
  # extract parameters
  sigmasq.e <- exp(param[paste0("lsigmasq.e", method.names[2])])
  delta <- param[paste0("delta", method.names[2])]
  
  # get error variances on vgrid
  if (varfunc == "power") {
    err.var <- sigmasq.e * (abs(vgrid)^(2 * delta))
  }
  if (varfunc == "exponential") {
    err.var <- sigmasq.e * exp(2 * delta * vgrid)
  }
  result <- log(sqrt(2 * err.var * q.chisq))
  return(result)
}


#########################################################################

#
# Function for computing log(TDI) on vgrid as a function of mean covariates
# Arguments: 
#		param = vector of transformed model parameters
#		varfunc = a character with one of two possible values---
#				"power" or "exponential"--- indicating the common variance function
#				model for the two methods
#		vgrid = a grid of values on which to evaluate the fitted error variance 
#					functions and the paramerters that involve them
#		method.names = a character vector of length two containing names of methods 
#		mean.covariates = logical indicate whether mean.covariates are included in the model
#		contr.mat = if mean.covariates are included, then contr.mat provides the matrix of 
#			contrast coefficients which when premultiplied to fixed(fit) provides the
#			vector of difference in means of two methods at each covariate setting
#		covariate.level = an integer indicating the given covariate level
#		prob = probability cutoff used in TDI's definition
#		ulim = upper limit of the interval that is expected to cotain the root that determines TDI
# Value: 
#		log(TDI) on vgrid for the given covariate level


ltdi.fun.reg <- function(param, varfunc, vgrid, method.names, mean.covariates, contr.mat, 
                         covariate.level, prob, ulim = 500) {
  # extract parameters
  sigmasq.b <- exp(param["lsigmasq.b"])
  names(sigmasq.b) <- "sigmasq.b"
  sigmasq.e <- exp(param[paste0("lsigmasq.e", method.names)])
  names(sigmasq.e) <- method.names
  delta <- param[paste0("delta", method.names)]
  names(delta) <- method.names
  
  # get error variances on vgrid
  err.var <- matrix(NA, nrow = length(vgrid), ncol = 2)
  colnames(err.var) <- method.names
  if (varfunc == "power") {
    err.var[, method.names[1]] <- sigmasq.e[method.names[1]] * (abs(vgrid)^(2 * delta[method.names[1]]))
    err.var[, method.names[2]] <- sigmasq.e[method.names[2]] * (abs(vgrid)^(2 * delta[method.names[2]]))
  }
  
  if (varfunc == "exponential") {
    err.var[, method.names[1]] <- sigmasq.e[method.names[1]] * exp(2 * delta[method.names[1]] * 
                                                                     vgrid)
    err.var[, method.names[2]] <- sigmasq.e[method.names[2]] * exp(2 * delta[method.names[2]] * 
                                                                     vgrid)
  }
  
  # get mean of D 
  Dmean <- mean.diff.fun.reg(param, method.names, mean.covariates, contr.mat) # may be a vector
  names(Dmean) <- NULL
  
  # get var and cov of Y
  Yvar <- sigmasq.b + err.var
  colnames(Yvar) <- method.names
  Ycov <- sigmasq.b
  
  # get SD of D on vgrid
  Dsd <- sqrt(Yvar[, method.names[1]] + Yvar[, method.names[2]] - 2 * Ycov)
  names(Dsd) <- NULL
  
  # get log(TDI) on vgrid for Dmean at the given covariate level and Dsd
  result <- mapply(ltdi, dmean = Dmean[covariate.level], dsd = Dsd, MoreArgs = list(prob = prob, 
                                                                                    ulim = ulim))
  
  return(result)
}


#########################################################################

#
# Function for computing log(precision ratio) as a function of model parameters on vgrid
# Arguments: 
#		param = transformed model parameters
#		varfunc = a character with one of two possible values---
#				"power" or "exponential"--- indicating the common variance function
#				model for the two methods
#		vgrid = a grid of values on which to evaluate the fitted error variance 
#					functions and the paramerters that involve them
#		method.names = a character vector of length two containing names of methods 
#	Value:
#		log(precision ration) on vgrid
#

lprec.ratio.fun.reg <- function(param, varfunc, vgrid, method.names) {
  
  sigmasq.e <- exp(param[paste0("lsigmasq.e", method.names)])
  names(sigmasq.e) <- method.names
  delta <- param[paste0("delta", method.names)]
  names(delta) <- method.names
  
  
  err.var <- matrix(NA, nrow = length(vgrid), ncol = 2)
  colnames(err.var) <- method.names
  if (varfunc == "power") {
    err.var[, method.names[1]] <- sigmasq.e[method.names[1]] * (abs(vgrid)^(2 * delta[method.names[1]]))
    err.var[, method.names[2]] <- sigmasq.e[method.names[2]] * (abs(vgrid)^(2 * delta[method.names[2]]))
  }
  
  if (varfunc == "exponential") {
    err.var[, method.names[1]] <- sigmasq.e[method.names[1]] * exp(2 * delta[method.names[1]] * 
                                                                     vgrid)
    err.var[, method.names[2]] <- sigmasq.e[method.names[2]] * exp(2 * delta[method.names[2]] * 
                                                                     vgrid)
  }
  
  result <- log(err.var[, method.names[1]]/err.var[, method.names[2]])
  
  return(result)
}

#########################################################################

#
# Function for computing mean difference as a function of model parameters
# Arguments: 
#		param = vector of transformed model parameters
#		method.names = a character vector of length two containing names of methods 
#		mean.covariates = logical indicate whether mean.covariates are included in the model
#		contr.mat = if mean.covariates are included, then contr.mat provides the matrix of 
#			contrast coefficients which when premultiplied to fixed(fit) provides the
#			vector of difference in means of two methods at each covariate setting
#	Value:
#		a vector containing values of mean difference function at given covariates
#

mean.diff.fun.reg <- function(param, method.names, mean.covariates = FALSE, contr.mat = NULL) {
  if (mean.covariates == FALSE) {
    beta <- param[paste0("method", method.names)]
    mean.diff <- beta[2] - beta[1]
    names(mean.diff) <- paste0(method.names[2], "-", method.names[1])
  }
  if (mean.covariates == TRUE) {
    beta <- param[1:ncol(contr.mat)]
    mean.diff <- as.numeric(contr.mat %*% beta)
  }
  return(mean.diff)
}

#########################################################################

#
# Function for extracting parameter estimates from the fitted object and computing
# their covariance matrix (inverse Hessian) and confidence intervals of parameters
# Arguments: 
#		dat = the data frame
#		vtilde = a vector of length nrow(dat) containing the subject-specific 
#				value of the variance covariate 
#		varfunc = a character with one of two possible values---
#				"power" or "exponential"--- indicating the common variance function
#				model for the two methods
#		vgrid = a grid of values on which to evaluate the fitted error variance 
#					functions and the paramerters that involve them
#		fm = fitted lme object
#		conf.level = level of confidence (sometimes individual, sometimes simultaneous)
# Value: 
#	A bunch of stuff including estimates of transformed model parameter vector, 
# 	inverse Hessian, etc.


estimates.noint.hetero.reg <- function(dat, vtilde, varfunc, fm, vgrid, conf.level = 0.95) {
  require(numDeriv)
  require(nlme)
  z.twosided <- qnorm((1 + conf.level)/2)
  
  # get method names
  method.names <- levels(droplevels(dat$method))
  
  # estimate of beta
  beta.hat <- fixef(fm)
  
  # extract random effect variances
  sigma.hat <- fm$sigma
  sigmasq.b.hat <- as.matrix(fm$modelStruct$reStruct[[1]])["(Intercept)", "(Intercept)"] * 
    sigma.hat^2
  
  # extract the remaining estimates
  sigma.ratio.hat <- coef(fm$modelStruct$varStruct$A, unconstrained = FALSE, allCoef = TRUE)[method.names]
  sigmasq.e.hat <- (sigma.hat * sigma.ratio.hat)^2
  
  # get the heteroscedasticity parameters
  delta.hat <- coef(fm$modelStruct$varStruct$B)[method.names]
  
  param.hat <- c(beta.hat, log(c(sigmasq.b.hat, sigmasq.e.hat)), delta.hat)
  names(param.hat) <- c(names(beta.hat), "lsigmasq.b", paste0("lsigmasq.e", method.names), 
                        paste0("delta", method.names))
  
  # get Hessian and invert it
  hess <- -hessian(loglik.noint.hetero.reg, param.hat, dat = dat, vtilde = vtilde, 
                   varfunc = varfunc, fm = fm)
  inv.hess <- solve(hess)
  rownames(inv.hess) <- names(param.hat)
  colnames(inv.hess) <- names(param.hat)
  
  # SEs of param.hat and confidence intevals
  se.param.hat <- sqrt(diag(inv.hess))
  ci.param <- cbind(param.hat - z.twosided * se.param.hat, param.hat + z.twosided * 
                      se.param.hat)
  results.param <- cbind(param.hat, se.param.hat, ci.param)
  rownames(results.param) <- names(param.hat)
  colnames(results.param) <- c("estimate", "se", "lcl", "ucl")
  
  # estimated error variances
  err.var.hat <- matrix(NA, nrow = length(vgrid), ncol = 2)
  colnames(err.var.hat) <- method.names
  if (varfunc == "power") {
    err.var.hat[, method.names[1]] <- sigmasq.e.hat[method.names[1]] * (abs(vgrid)^(2 * 
                                                                                      delta.hat[method.names[1]]))
    err.var.hat[, method.names[2]] <- sigmasq.e.hat[method.names[2]] * (abs(vgrid)^(2 * 
                                                                                      delta.hat[method.names[2]]))
  }
  
  if (varfunc == "exponential") {
    err.var.hat[, method.names[1]] <- sigmasq.e.hat[method.names[1]] * exp(2 * delta.hat[method.names[1]] * 
                                                                             vgrid)
    err.var.hat[, method.names[2]] <- sigmasq.e.hat[method.names[2]] * exp(2 * delta.hat[method.names[2]] * 
                                                                             vgrid)
  }
  
  # fitted bivariate distribution
  var.hat <- sigmasq.b.hat + err.var.hat
  colnames(var.hat) <- method.names
  corr.hat <- sigmasq.b.hat/sqrt(var.hat[, 1] * var.hat[, 2])
  bvn.hat <- list(beta = beta.hat, var = var.hat, corr = corr.hat)
  
  # fitted variance for the difference
  var.D.hat <- err.var.hat[, 1] + err.var.hat[, 2]
  
  return(list(method.names = method.names, param.hat = param.hat, inv.hess = inv.hess, 
              results.param = results.param, vgrid = vgrid, err.var.hat = err.var.hat, bvn.hat = bvn.hat, 
              var.D.hat = var.D.hat))
}

##################################################################################

#
# Function for computing loglikelihood function as a function of model parameters
# for unlinked repeated measures data without method x subject interaction
# Arguments: 
#		param = transformed model parameters
#		dat = a data frame whose columns include --- subject (a factor), 
#			  method (a factor with two levels), meas, mean covariates
#		vtilde = a vector of length nrow(dat) containing the subject-specific 
#				values of the variance covariate 
#		varfunc = a character with one of two possible values---
#				"power" or "exponential"--- indicating the common variance function
#				model for the two methods
#		fm = fitted model object
# Value: The loglikelihood function


loglik.noint.hetero.reg <- function(param, dat, vtilde, varfunc, fm) {
  require(mvtnorm)
  require(nlme)
  
  # function for computing loglik for one subject
  loglik.one <- function(dd, beta, G, varfunc, fm) {
    
    # design matrices
    X <- model.matrix(fm, data = dd)
    Z <- model.matrix(fm$modelStruct$reStruct, data = dd)
    
    # variance covariate
    vv <- unique(dd$vtilde)
    if (length(vv) != 1) {
      stop("variance covariate in loglik.noint.hetero.reg is not unique to a subject")
    }
    
    # error variance matrix		
    if (varfunc == "power") {
      err.var.1 <- sigmasq.e[method.names[1]] * (abs(vv)^(2 * delta[method.names[1]]))
      err.var.2 <- sigmasq.e[method.names[2]] * (abs(vv)^(2 * delta[method.names[2]]))
    }
    if (varfunc == "exponential") {
      err.var.1 <- sigmasq.e[method.names[1]] * exp(2 * delta[method.names[1]] * 
                                                      vv)
      err.var.2 <- sigmasq.e[method.names[2]] * exp(2 * delta[method.names[2]] * 
                                                      vv)
    }
    
    R <- diag(err.var.1 * (dd$method == method.names[1]) + err.var.2 * (dd$method == 
                                                                          method.names[2]), nrow = nrow(X))
    
    # loglikelihood
    return(dmvnorm(dd$meas, mean = (X %*% beta), sigma = (Z %*% G %*% t(Z) + R), 
                   log = T))
  }
  
  if (!(varfunc %in% c("power", "exponential"))) {
    stop(paste(varfunc, "is not a valid option in loglik.hetero.repmeas.data"))
  }
  
  # get method names
  method.names <- levels(droplevels(dat$method))
  
  # parameters
  beta <- param[names(fixef(fm))]
  sigmasq.b <- exp(param["lsigmasq.b"])
  sigmasq.e <- exp(param[paste0("lsigmasq.e", method.names)])
  names(sigmasq.e) <- method.names
  delta <- param[paste0("delta", method.names)]
  names(delta) <- method.names
  
  G <- sigmasq.b
  
  # split data and variance covariate 
  dat <- data.frame(dat, vtilde = vtilde)
  datlist <- split(dat, dat$subject)
  
  # compute loglik for each subject
  loglik <- mapply(loglik.one, datlist, MoreArgs = list(beta = beta, G = G, varfunc = varfunc, 
                                                        fm = fm))
  
  return(sum(as.numeric(loglik)))
}


####################################################################################
# chapter 7 #
####################################################################################

# begin some code from an older version of R package multcomp
# (needed for computing simultaneous intervals in Chapter 7)

chkcorr <- 
  function (x) 
  {
    if (!is.matrix(x)) 
      return(FALSE)
    rownames(x) <- colnames(x) <- NULL
    storage.mode(x) <- "numeric"
    ONE <- 1 + sqrt(.Machine$double.eps)
    ret <- (min(x) < -ONE || max(x) > ONE) || !isTRUE(all.equal(diag(x), 
                                                                rep(1, nrow(x))))
    !ret
  }

checkmvArgs <- 
  function (lower, upper, mean, corr, sigma) 
  {
    UNI <- FALSE
    if (!is.numeric(lower) || !is.vector(lower)) 
      stop(sQuote("lower"), " is not a numeric vector")
    if (!is.numeric(upper) || !is.vector(upper)) 
      stop(sQuote("upper"), " is not a numeric vector")
    if (!is.numeric(mean) || !is.vector(mean)) 
      stop(sQuote("mean"), " is not a numeric vector")
    if (is.null(lower) || any(is.na(lower))) 
      stop(sQuote("lower"), " not specified or contains NA")
    if (is.null(upper) || any(is.na(upper))) 
      stop(sQuote("upper"), " not specified or contains NA")
    rec <- cbind(lower, upper, mean)
    lower <- rec[, "lower"]
    upper <- rec[, "upper"]
    if (!all(lower <= upper)) 
      stop("at least one element of ", sQuote("lower"), " is larger than ", 
           sQuote("upper"))
    mean <- rec[, "mean"]
    if (any(is.na(mean))) 
      stop("mean contains NA")
    if (is.null(corr) && is.null(sigma)) {
      corr <- diag(length(lower))
    }
    if (!is.null(corr) && !is.null(sigma)) {
      sigma <- NULL
      warning("both ", sQuote("corr"), " and ", sQuote("sigma"), 
              " specified: ignoring ", sQuote("sigma"))
    }
    if (!is.null(corr)) {
      if (!is.numeric(corr)) 
        stop(sQuote("corr"), " is not numeric")
      if (!is.matrix(corr)) {
        if (length(corr) == 1) 
          UNI <- TRUE
        if (length(corr) != length(lower)) 
          stop(sQuote("diag(corr)"), " and ", sQuote("lower"), 
               " are of different length")
      }
      else {
        if (length(corr) == 1) {
          UNI <- TRUE
          corr <- corr[1, 1]
          if (length(lower) != 1) 
            stop(sQuote("corr"), " and ", sQuote("lower"), 
                 " are of different length")
        }
        else {
          if (length(diag(corr)) != length(lower)) 
            stop(sQuote("diag(corr)"), " and ", sQuote("lower"), 
                 " are of different length")
          if (!chkcorr(corr)) 
            stop(sQuote("corr"), " is not a correlation matrix")
        }
      }
    }
    if (!is.null(sigma)) {
      if (!is.numeric(sigma)) 
        stop(sQuote("sigma"), " is not numeric")
      if (!is.matrix(sigma)) {
        if (length(sigma) == 1) 
          UNI <- TRUE
        if (length(sigma) != length(lower)) 
          stop(sQuote("diag(sigma)"), " and ", sQuote("lower"), 
               " are of different length")
      }
      else {
        if (length(sigma) == 1) {
          UNI <- TRUE
          sigma <- sigma[1, 1]
          if (length(lower) != 1) 
            stop(sQuote("sigma"), " and ", sQuote("lower"), 
                 " are of different length")
        }
        else {
          if (length(diag(sigma)) != length(lower)) 
            stop(sQuote("diag(sigma)"), " and ", sQuote("lower"), 
                 " are of different length")
          if (!isTRUE(all.equal(sigma, t(sigma))) || any(diag(sigma) < 
                                                         0)) 
            stop(sQuote("sigma"), " is not a covariance matrix")
        }
      }
    }
    list(lower = lower, upper = upper, mean = mean, corr = corr, 
         sigma = sigma, uni = UNI)
  }



dots2GenzBretz <- 
  function (...) 
  {
    addargs <- list(...)
    fm1 <- sapply(names(addargs), function(x) length(grep(x, 
                                                          names(formals(GenzBretz)))) == 1)
    fm2 <- sapply(names(addargs), function(x) length(grep(x, 
                                                          names(formals(uniroot)))) == 1)
    algorithm <- NULL
    uniroot <- NULL
    if (any(fm1)) 
      algorithm <- do.call("GenzBretz", addargs[fm1])
    if (any(fm2)) 
      uniroot <- addargs[fm2]
    list(algorithm = algorithm, uniroot = uniroot)
  }



chkdots <- 
  function (...) 
  {
    lst <- list(...)
    if (length(lst) > 0) {
      warning("Argument(s) ", sQuote(names(lst)), " passed to ", 
              sQuote("..."), " are ignored", call. = TRUE)
    }
  }


confint.glht <- 
  function (object, parm, level = 0.95, calpha = adjusted_calpha(), 
            ...) 
  {
    chkdots(...)
    type <- attr(calpha, "type")
    if (is.function(calpha)) 
      calpha <- calpha(object, level)
    if (!is.numeric(calpha) || length(calpha) != 1) 
      stop(sQuote("calpha"), " is not a scalar")
    error <- attr(calpha, "error")
    attributes(calpha) <- NULL
    betahat <- coef(object)
    ses <- sqrt(diag(vcov(object)))
    switch(object$alternative, two.sided = {
      LowerCL <- betahat - calpha * ses
      UpperCL <- betahat + calpha * ses
    }, less = {
      LowerCL <- rep(-Inf, length(ses))
      UpperCL <- betahat + calpha * ses
    }, greater = {
      LowerCL <- betahat + calpha * ses
      UpperCL <- rep(Inf, length(ses))
    })
    ci <- cbind(LowerCL, UpperCL)
    colnames(ci) <- c("lower", "upper")
    object$confint <- cbind(betahat, ci)
    colnames(object$confint) <- c("Estimate", "lwr", "upr")
    attr(object$confint, "conf.level") <- level
    attr(object$confint, "calpha") <- calpha
    attr(object$confint, "error") <- error
    if (is.null(type)) 
      type <- "univariate"
    attr(object, "type") <- type
    class(object) <- c("confint.glht", "glht")
    return(object)
  }


qmvt <- 
  function (p, 
            tail = c("lower.tail", "upper.tail", 
                     "both.tails"), df = 1, delta = 0, corr = NULL, sigma = NULL, 
            algorithm = GenzBretz(), type = c("Kshirsagar", "shifted"), 
            ...) 
  {
    # added by PKC 
    if (tail == "lower.tail") {interval <- c(0, 4)}
    if (tail == "upper.tail") {interval <- c(-4, 0)}
    if (tail == "both.tails") {interval <- c(0, 4)}
    
    
    if (length(p) != 1 || (p <= 0 || p >= 1)) 
      stop(sQuote("p"), " is not a double between zero and one")
    dots <- dots2GenzBretz(...)
    if (!is.null(dots$algorithm) && !is.null(algorithm)) 
      algorithm <- dots$algorithm
    type <- match.arg(type)
    tail <- match.arg(tail)
    if (tail == "both.tails" && p < 0.5) 
      stop("cannot compute two-sided quantile for p < 0.5")
    dim <- 1
    if (!is.null(corr)) 
      dim <- NROW(corr)
    if (!is.null(sigma)) 
      dim <- NROW(sigma)
    lower <- upper <- rep.int(0, dim)
    args <- checkmvArgs(lower, upper, delta, corr, sigma)
    if (args$uni) {
      if (tail == "both.tails") 
        p <- ifelse(p < 0.5, p/2, 1 - (1 - p)/2)
      if (df == 0 || isInf(df)) {
        q <- qnorm(p, mean = args$mean, lower.tail = (tail != 
                                                        "upper.tail"))
      }
      else {
        q <- qt(p, df = df, ncp = args$mean, lower.tail = (tail != 
                                                             "upper.tail"))
      }
      qroot <- list(quantile = q, f.quantile = 0)
      return(qroot)
    }
    dim <- length(args$mean)
    pfct <- function(q) {
      switch(tail, both.tails = {
        low <- rep(-abs(q), dim)
        upp <- rep(abs(q), dim)
      }, upper.tail = {
        low <- rep(q, dim)
        upp <- rep(Inf, dim)
      }, lower.tail = {
        low <- rep(-Inf, dim)
        upp <- rep(q, dim)
      } )
      pmvt(lower = low, upper = upp, df = df, delta = args$mean, 
           corr = args$corr, sigma = args$sigma, algorithm = algorithm, 
           type = type) - p
    }
    if (is.null(interval) || length(interval) != 2) {
      if (FALSE) {
        cr <- args$corr
        interval <- approx_interval(p = p, tail = tail, corr = cr, 
                                    df = df)
      }
      else {
        interval <- c(if (tail == "both.tails") 0 else -10, 
                      10)
      }
    }
    qroot <- if (is.null(dots$uniroot)) {
      uniroot(pfct, interval = interval, extendInt = "yes")
    }
    else {
      do.call("uniroot", list(pfct, interval = interval, dots$uniroot))
    }
    names(qroot)[1:2] <- c("quantile", "f.quantile")
    qroot
  }

adjusted_calpha <- 
  function (...) 
  {
    ret <- function(object, level) {
      pqglht(object)$qfunction(level, adjusted = TRUE, ...)
    }
    attr(ret, "type") <- "adjusted"
    ret
  }


pqglht <- function (object) 
{
  betahat <- coef(object)
  covm <- vcov(object)
  m <- coef(object, rhs = TRUE)
  df <- object$df
  ses <- sqrt(diag(covm))
  tstat <- (betahat - m)/ses
  cr <- cov2cor(covm)
  dim <- ncol(cr)
  pfunction <- function(type = c("univariate", "adjusted", 
                                 p.adjust.methods), ...) {
    type <- match.arg(type)
    pfct <- function(q) {
      switch(object$alternative, two.sided = {
        low <- rep(-abs(q), dim)
        upp <- rep(abs(q), dim)
      }, less = {
        low <- rep(q, dim)
        upp <- rep(Inf, dim)
      }, greater = {
        low <- rep(-Inf, dim)
        upp <- rep(q, dim)
      })
      pmvt(lower = low, upper = upp, df = df, corr = cr, 
           ...)
    }
    switch(object$alternative, two.sided = {
      if (df > 0) pvals <- 2 * (1 - pt(abs(tstat), df)) else pvals <- 2 * 
          (1 - pnorm(abs(tstat)))
    }, less = {
      if (df > 0) pvals <- pt(tstat, df) else pvals <- pnorm(tstat)
    }, greater = {
      if (df > 0) pvals <- 1 - pt(tstat, df) else pvals <- 1 - 
          pnorm(tstat)
    })
    if (type == "univariate") 
      return(pvals)
    if (type == "adjusted") {
      ret <- numeric(length(tstat))
      error <- 0
      for (i in 1:length(tstat)) {
        tmp <- pfct(tstat[i])
        if (attr(tmp, "msg") != "Normal Completion" && 
            length(grep("^univariate", attr(tmp, "msg"))) == 
            0) 
          warning(attr(tmp, "msg"))
        if (error < attr(tmp, "error")) 
          error <- attr(tmp, "error")
        ret[i] <- tmp
      }
      ret <- 1 - ret
      attr(ret, "error") <- error
      return(ret)
    }
    return(p.adjust(pvals, method = type))
  }
  qfunction <- function(conf.level, adjusted = TRUE, ...) {
    tail <- switch(object$alternative, two.sided = "both.tails", 
                   less = "lower.tail", greater = "upper.tail")
    if (adjusted) {
      calpha <- qmvt(conf.level, df = df, corr = cr, tail = tail, 
                     ...)
    }
    else {
      calpha <- qmvt(conf.level, df = df, corr = matrix(1), 
                     tail = tail, ...)
    }
    ret <- calpha$quantile
    attr(ret, "error") <- calpha$estim.prec
    return(ret)
  }
  RET <- list(pfunction = pfunction, qfunction = qfunction, 
              coefficients = betahat, sigma = ses, tstat = tstat)
  class(RET) <- "pqglht"
  RET
}

# end some code from an older version of R package multcomp

####################################################################################
# chapter 7 (continued)
####################################################################################

####################
#
# Function for computing estimates, their SEs and confidence intervals
#		 for various parameters and functions for unreplicated data
#		from multiple methods
# Arguments: 
#		dat = a data frame with four columns --- subject (a factor), 
#			  method (a factor with two levels method1 and method2), 
#			rep (replication number, unused) and meas
#		param.hat = vector of ML estimates of transformed parameters
# Value: 
#		standard errors of estimates of parameters and their confidence intervals
#		estimate, se and ucb of log(TDI), and estimate and ucb of TDI
#		estimate, se and lcb of z(CCC), and estimate and lcb of CCC
#		estimates, se and ucb of log(TDI.1) and log(TDI.2), and estimates and ucb for 
#				TDI.1 and TDI.2
#		estimates, se and lcb for z(CCC.1) and z(CCC.2), and estimates and lcb for 
#				CCC.1 and CCC.2
#

conf.measures.unreplicated.multi <- function(dat, param.hat, method.names, 
                                             contrast.matrix, conf.level = 0.95, prob = 0.9, seed=12345) {
  require(numDeriv)
  require(multcomp)
  set.seed(seed) # multcomp uses a randomized algorithm to compute critical point
  z.twosided <- qnorm((1 + conf.level)/2)
  
  # get Hessian and invert it
  hess <- -hessian(loglik.unreplicated.data.multi, param.hat, dat = dat, method.names = method.names)
  inv.hess <- solve(hess)
  rownames(inv.hess) <- names(param.hat)
  colnames(inv.hess) <- names(param.hat)
  
  # SEs of param.hat and confidence intevals
  se.param.hat <- sqrt(diag(inv.hess))
  ci.param <- cbind(param.hat - z.twosided * se.param.hat, param.hat + z.twosided * 
                      se.param.hat)
  ci.param <- cbind(param.hat, se.param.hat, ci.param)
  rownames(ci.param) <- names(param.hat)
  colnames(ci.param) <- c("estimate", "se", "lcl", "ucl")
  
  # evaluation of similarity 
  # estimates and simultaneous confidence intervals for similarity measures
  # means
  ci.mean <- ci.simul.paramfun(param.hat, inv.hess, mean.fun, contrast.matrix, 
                               "two.sided", conf.level, method.names = method.names)
  
  # precisions (precision ratio = lambda)
  ci.lprec <- ci.simul.paramfun(param.hat, inv.hess, lprec.fun.multi, contrast.matrix, 
                                "two.sided", conf.level, method.names = method.names, include.sigmasq.bI = FALSE)
  result1 <- exp(ci.lprec$result[, c("estimate", "lcl", "ucl"), drop = FALSE])
  rownames(result1) <- gsub("-", rep = "/", rownames(contrast.matrix))
  ci.lprec <- list(ci.fun = ci.lprec$ci.fun, lresult=ci.lprec$result, result=result1)
  
  # evaluation of agreement 
  # estimates and simultaneous confidence intervals for agreement measures
  contrast.matrix.agr.meas <- diag(nrow(contrast.matrix))
  rownames(contrast.matrix.agr.meas) <- gsub("-", rep = ",", rownames(contrast.matrix))
  
  # tdi
  ci.ltdi <- ci.simul.paramfun(param.hat, inv.hess, ltdi.unreplicated.fun.multi, 
                               contrast.matrix.agr.meas, "less", conf.level, method.names = method.names, 
                               contrast.matrix = contrast.matrix, prob = prob)
  ci.ltdi <- list(ci.fun = ci.ltdi$ci.fun, lresult = ci.ltdi$result, 
                  result = exp(ci.ltdi$result[, c("estimate", "lcl", "ucl"), drop = FALSE]))
  
  # ccc
  ci.zccc <- ci.simul.paramfun(param.hat, inv.hess, zccc.unreplicated.fun.multi, 
                               contrast.matrix.agr.meas, "greater", conf.level, method.names = method.names, 
                               contrast.matrix = contrast.matrix)
  ci.zccc <- list(ci.fun = ci.zccc$ci.fun, zresult = ci.zccc$result, 
                  result = tanh(ci.zccc$result[, c("estimate", "lcl", "ucl"), drop = FALSE]))
  
  results <- list(param = ci.param, mean.diff = ci.mean, precision.ratio = ci.lprec, 
                  tdi = ci.ltdi, ccc = ci.zccc)
  return(results)
}

######################
#
# Function for computing Fisher's z-transformation of CCC as a 
#		function of model parameters for unreplicated data
# Arguments: 
#		param = transformed model parameters
#		method.names = levels of the factor method
# 		contrast.matrix = a matrix when multiplied by the vector Y gives
#			the vector of differences of interest for multiple comparisons
# Value:
#		A vector of values of z(CCC)
#

zccc.unreplicated.fun.multi <- function(param, method.names, contrast.matrix) {
  # get the mean and variance of Y vector
  Ymean <- mean.fun(param, method.names)
  sigmasq.b <- exp(param["lsigmasq.b"])
  sigmasq.e <- exp(param[paste0("lsigmasq.e", method.names)])
  Yvar <- sigmasq.b + diag(sigmasq.e) 
  rownames(Yvar) <- method.names
  colnames(Yvar) <- method.names
  
  # get method pairs of interest
  method.pairs <- matrix(unlist(strsplit(rownames(contrast.matrix), 
                                         "-")), byrow = T, ncol = 2)		
  rownames(method.pairs) <- rownames(contrast.matrix)
  
  result <- mapply(function(i, j) {
    zccc(Ymean[i] - Ymean[j], Yvar[i, i], Yvar[j, j], Yvar[i, j])
  }, i = method.pairs[, 1], j = method.pairs[, 2])
  names(result) <- rownames(contrast.matrix)
  return(result)
}

######################
#
# Function for computing log(TDI) as a function of model parameters for unreplicated data 
# Arguments: 
#		param = transformed model parameters
#		method.names = levels of the factor method
# 		contrast.matrix = a matrix when multiplied by the vector Y gives
#			the vector of differences of interest for multiple comparisons
# 		prob = probability
#		ulim = upper limit of the interval in which to search for the root 
# Value:
#		A vector of values of log(TDI)
#
ltdi.unreplicated.fun.multi <- function(param, method.names, contrast.matrix, 
                                        prob, ulim = 500) {
  # get the mean and variance of Y vector
  Ymean <- mean.fun(param, method.names)
  sigmasq.b <- exp(param["lsigmasq.b"])
  sigmasq.e <- exp(param[paste0("lsigmasq.e", method.names)])
  Yvar <- sigmasq.b + diag(sigmasq.e) 
  rownames(Yvar) <- method.names
  colnames(Yvar) <- method.names
  
  # get the mean and standard deviations of the elements of D vector
  Dmean <- as.numeric(contrast.matrix %*% Ymean)
  Dsd <- sqrt(diag(contrast.matrix %*% Yvar %*% t(contrast.matrix)))
  
  # get log(TDI)
  result <- mapply(ltdi, dmean = Dmean, dsd = Dsd, MoreArgs = list(prob = prob, 
                                                                   ulim = ulim))
  names(result) <- rownames(contrast.matrix)
  return(result)
}


######################
#
# Function for computing loglikelihood function as a function of model parameters
# for unreplicated data
# Arguments: 
#		param = transformed model parameters
#		method.names = levels of the factor method
# 		dat = data in long format, used for model fitting
# Value:
#		loglikelihood function

loglik.unreplicated.data.multi <- function(param, method.names, dat) {
  require(mvtnorm)
  # function for computing loglik for one subject
  loglik.unlinked.one <- function(yy) {
    return(dmvnorm(yy, mean=beta, sigma=(sigmasq.b + diag(sigmasq.e)), log=T))
  }
  
  J <- length(method.names)
  if (J <= 2) {stop("the factor method should have more than two levels")}
  
  if (J > 2) {
    alpha <- param[paste0("alpha", method.names[2:J])]
    mu.b <- param["mu.b"]
    sigmasq.b <- exp(param["lsigmasq.b"])
    sigmasq.e <- exp(param[paste0("lsigmasq.e", method.names)])
    names(sigmasq.e) <- method.names
  }
  
  # parameters
  beta <- c(mu.b, alpha + mu.b)
  names(beta) <- method.names
  
  # get data in wide format
  ylist <- long2wide.unlinked(dat)
  
  # compute loglik for each subject
  loglik <- lapply(ylist, loglik.unlinked.one)
  
  return(sum(as.numeric(loglik)))
}


##################
#
# Function for fitting a mixed model to unreplicated method comparison measurements
# data with multiple (J > 2 ) methods (note strict inequality)
# Arguments: 
#		dat = a data frame with three columns --- subject (a factor), 
#			  	method (a factor), and meas (a numeric)
#		comp.type = a character with two possible values ("all.pairwise" or 
#			"reference") indicating the type of comparisons (default = "reference")
# Value: A list with following components
#		fit = fitted lme object
#		param.hat = ML estimates of transformed parameters
#		var.hat = ML estimates of variance components and variance of difference
#		ydist.hat = parameters of fitted multivariate normal distribution
#		dist.hat = parameters of the fitted distribution of the differences
#		reliability.hat = ML estimates of reliabilities of the two methods
#		method.names = a character vector of names of methods
#		contrast.matrix = the contrast matrix determined by comp.type
#

lme.unreplicated.fit.multi <- function(dat, comp.type = "reference") {
  require(nlme)
  method.names.used <- droplevels(dat$method)
  J <- nlevels(method.names.used)
  method.names <- levels(method.names.used)
  
  if (J <= 2) {stop("the factor method should have more than two levels")}
  
  # fit model
  dat <- groupedData(meas ~ method | subject, data = dat, order.groups = FALSE)
  fit <- lme(meas ~ method - 1, data = dat, weights = varIdent(form = ~1 | 
                                                                 method), method = "ML", random = ~1, 
             control = lmeControl(maxIter = 1000, msMaxIter = 5000, niterEM = 2000, 
                                  msMaxEval = 1500))
  
  # extract random effect variances
  sigma.hat <- fit$sigma
  mat1 <- diag(as.matrix(fit$modelStruct$reStruct[[1]]) * sigma.hat^2)
  sigmasq.b.hat <- mat1["(Intercept)"]
  
  # extract the remaining estimates
  mu.hat <- fixef(fit)[paste0(rep("method", J), method.names)]
  names(mu.hat) <- method.names
  
  sigma.ratio.hat <- coef(fit$modelStruct$varStruct, unconstrained = FALSE, 
                          allCoef = TRUE)[method.names]
  sigmasq.e.hat <- (sigma.hat * sigma.ratio.hat)^2
  
  # cast estimates in the desired form
  mu.b.hat <- mu.hat[method.names[1]]
  alpha.hat <- mu.hat[method.names[2:J]] - mu.hat[method.names[1]] # a (J-1)-vector
  
  param.hat <- c(alpha.hat, mu.b.hat, log(sigmasq.b.hat), log(sigmasq.e.hat))
  names(param.hat) <- c(paste0("alpha", method.names[2:J]), "mu.b", "lsigmasq.b", 
                        paste0("lsigmasq.e", method.names))
  
  # estimated variance parameters
  varpar.hat <- list(sigmasq.b = sigmasq.b.hat, sigmasq.e = sigmasq.e.hat)
  
  # fitted J-variate distribution 
  var.hat <- sigmasq.b.hat + diag(sigmasq.e.hat)
  rownames(var.hat) <- method.names
  colnames(var.hat) <- method.names
  ydist.hat <- list(ymean = mu.hat, ycov = var.hat, ycor = cov2cor(var.hat), 
                    yvar = diag(var.hat))
  
  # fitted distribution of difference vector of interest 
  contrast.matrix <- make.contrast.matrix(method.names, comp.type)
  dmean.hat <- as.numeric(contrast.matrix %*% mu.hat)
  names(dmean.hat) <- rownames(contrast.matrix)
  dvar.hat <- contrast.matrix %*% var.hat %*% t(contrast.matrix)
  ddist.hat <- list(dmean = dmean.hat, dcov = dvar.hat, dvar = diag(dvar.hat))
  
  # reliability
  rlbty.hat <- sigmasq.b.hat/diag(var.hat)
  names(rlbty.hat) <- method.names
  
  return(list(fit = fit, method.names = method.names, param.hat = param.hat, 
              varpar.hat = varpar.hat, ydist.hat = ydist.hat, contrast.matrix = contrast.matrix, 
              ddist.hat = ddist.hat, reliability.hat = rlbty.hat))
}

####################################

#
# Function for computing estimates, their SEs and confidence intervals
#		 for various parameters and functions for unlinked repeated measurements data
#		from multiple (>= 2) methods
# Arguments: 
#		dat = a data frame with four columns --- subject (a factor), 
#			  method (a factor with J (>=2) levels ), 
#			rep (replication number, unused) and meas
#		param.hat = vector of ML estimates of transformed parameters
#		method.names = levels of the factor method
# 		contrast.matrix = a matrix when multiplied by the vector Y gives
#			the vector of differences needed for agreement evaluation in pairs of methods
#		conf.level =level of confidence 
# 		prob = probability for TDI
#		seed = seed for random number generator involved in comuptation for 
#			critical point needed for simultaneous intervals
# Value: 
#		A bunch of lists (need to specify)
#

conf.measures.unlinked.multi <- function(dat, param.hat, method.names, 
                                         contrast.matrix, conf.level = 0.95, prob = 0.9, seed=12345) {
  require(numDeriv)
  require(multcomp)
  set.seed(seed) # multcomp uses a randomized algorithm to compute critical point
  z.twosided <- qnorm((1 + conf.level)/2)
  
  # get Hessian and invert it
  hess <- -hessian(loglik.unlinked.data.multi, param.hat, dat = dat, method.names = method.names)
  inv.hess <- solve(hess)
  rownames(inv.hess) <- names(param.hat)
  colnames(inv.hess) <- names(param.hat)
  
  # SEs of param.hat and confidence intevals
  se.param.hat <- sqrt(diag(inv.hess))
  ci.param <- cbind(param.hat - z.twosided * se.param.hat, param.hat + z.twosided * 
                      se.param.hat)
  ci.param <- cbind(param.hat, se.param.hat, ci.param)
  rownames(ci.param) <- names(param.hat)
  colnames(ci.param) <- c("estimate", "se", "lcl", "ucl")
  
  # evaluation of similarity 
  # estimates and simultaneous confidence intervals for similarity measures
  # means
  ci.mean <- ci.simul.paramfun(param.hat, inv.hess, mean.fun, contrast.matrix, 
                               "two.sided", conf.level, method.names = method.names)
  
  # precisions (precision ratio = lambda)
  # don't include interaction variance
  ci.lprec <- ci.simul.paramfun(param.hat, inv.hess, lprec.fun.multi, contrast.matrix, 
                                "two.sided", conf.level, method.names = method.names, include.sigmasq.bI = FALSE)
  result1 <- exp(ci.lprec$result[,c("estimate", "lcl", "ucl"), drop = FALSE])
  rownames(result1) <- gsub("-", rep = "/", rownames(contrast.matrix))
  ci.lprec <- list(ci.fun=ci.lprec$ci.fun, lresult=ci.lprec$result, result=result1)
  
  # include interaction variance
  ci.lnewprec <- ci.simul.paramfun(param.hat, inv.hess, lprec.fun.multi, contrast.matrix, 
                                   "two.sided", conf.level, method.names = method.names, include.sigmasq.bI = TRUE)
  result2 <- exp(ci.lnewprec$result[, c("estimate", "lcl", "ucl"), drop = FALSE])
  rownames(result2) <- gsub("-", rep = "/", rownames(contrast.matrix))
  ci.lnewprec <- list(ci.fun = ci.lnewprec$ci.fun, lresult = ci.lnewprec$result, 
                      result = result2)
  
  # evaluation of agreement 
  # estimates and simultaneous confidence intervals for agreement measures
  contrast.matrix.agr.meas <- diag(nrow(contrast.matrix))
  rownames(contrast.matrix.agr.meas) <- gsub("-", rep = ",", rownames(contrast.matrix))
  
  # tdi
  ci.ltdi <- ci.simul.paramfun(param.hat, inv.hess, ltdi.fun.multi, contrast.matrix.agr.meas, 
                               "less", conf.level, method.names = method.names, contrast.matrix = contrast.matrix, 
                               prob = prob)
  ci.ltdi <- list(ci.fun = ci.ltdi$ci.fun, lresult = ci.ltdi$result, result = exp(ci.ltdi$result[, 
                                                                                                 c("estimate", "lcl", "ucl"), drop = FALSE]))
  
  # ccc
  ci.zccc <- ci.simul.paramfun(param.hat, inv.hess, zccc.fun.multi, contrast.matrix.agr.meas, 
                               "greater", conf.level, method.names = method.names, contrast.matrix = contrast.matrix)
  ci.zccc <- list(ci.fun = ci.zccc$ci.fun, zresult = ci.zccc$result, result = tanh(ci.zccc$result[, 
                                                                                                  c("estimate", "lcl", "ucl"), drop = FALSE]))
  
  # evaluation of repeatability 
  # estimates and **individual** confidence intervals for repeatability measures
  contrast.matrix.rep.meas <- diag(length(method.names))
  rownames(contrast.matrix.rep.meas) <- method.names
  
  # tdi
  ci.rep.ltdi <- ci.ind.paramfun(param.hat, inv.hess, rep.ltdi.fun, contrast.matrix.rep.meas, 
                                 "less", conf.level, method.names = method.names, prob = prob)
  ci.rep.ltdi <- list(ci.fun = ci.rep.ltdi$ci.fun, lresult = ci.rep.ltdi$result, 
                      result = exp(ci.rep.ltdi$result[, c("estimate", "lcl", "ucl"), drop = FALSE]))
  
  # ccc
  ci.rep.zccc <- ci.ind.paramfun(param.hat, inv.hess, rep.zccc.fun, contrast.matrix.rep.meas, 
                                 "greater", conf.level, method.names = method.names)
  ci.rep.zccc <- list(ci.fun = ci.rep.zccc$ci.fun, zresult = ci.rep.zccc$result, 
                      result = tanh(ci.rep.zccc$result[, c("estimate", "lcl", "ucl"), drop = FALSE]))
  
  results <- list(param = ci.param, mean.diff = ci.mean, precision.ratio = ci.lprec, 
                  new.precision.ratio = ci.lnewprec, tdi = ci.ltdi, ccc = ci.zccc, rep.tdi = ci.rep.ltdi, 
                  rep.ccc = ci.rep.zccc)
  return(results)
}

######################
#
# Function for computing Fisher's z-transformation of repeatability CCC as a 
#		function of model parameters
# Arguments: 
#		param = transformed model parameters
#		method.names = levels of the factor method
# Value:
#		A vector of values of repeatability z(CCC)
#

rep.zccc.fun <- function(param, method.names) {
  J <- length(method.names)
  sigmasq.b <- exp(param["lsigmasq.b"])
  sigmasq.e <- exp(param[paste0("lsigmasq.e", method.names)])
  if (J == 2) {
    sigmasq.bI <- exp(param["lsigmasq.bI"])
  }
  if (J > 2) {
    sigmasq.bI <- exp(param[paste0("lsigmasq.bI", method.names)])
  }
  result <- atanh((sigmasq.b + sigmasq.bI)/(sigmasq.b + sigmasq.bI + sigmasq.e))
  names(result) <- method.names
  return(result)
}

######################
#
# Function for computing repeatability log(TDI) as a function of model parameters
# Arguments: 
#		param = transformed model parameters
#		method.names = levels of the factor method
# 		prob = probability
# Value:
#		A vector of values of repetability log(TDI)
#
rep.ltdi.fun <- function(param, method.names, prob) {
  q.chisq <- qchisq(prob, 1)
  result <- (1/2) * (log(2 * q.chisq) + param[paste0("lsigmasq.e", method.names)])
  names(result) <- method.names
  return(result)
}


######################
#
# Function for computing Fisher's z-transformation of CCC as a 
#		function of model parameters
# Arguments: 
#		param = transformed model parameters
#		method.names = levels of the factor method
# 		contrast.matrix = a matrix when multiplied by the vector Y gives
#			the vector of differences of interest for multiple comparisons
# Value:
#		A vector of values of z(CCC)
#

zccc.fun.multi <- function(param, method.names, contrast.matrix) {
  J <- length(method.names)
  Z <- cbind(rep(1, J), diag(J))
  
  # get the mean and variance of Y vector
  Ymean <- mean.fun(param, method.names)
  sigmasq.b <- exp(param["lsigmasq.b"])
  sigmasq.e <- exp(param[paste0("lsigmasq.e", method.names)])
  if (J == 2) {
    sigmasq.bI <- rep(exp(param["lsigmasq.bI"]), J)
  }
  if (J > 2) {
    sigmasq.bI <- exp(param[paste0("lsigmasq.bI", method.names)])
  }
  Yvar <- Z %*% diag(c(sigmasq.b, sigmasq.bI)) %*% t(Z) + diag(sigmasq.e)
  rownames(Yvar) <- method.names
  colnames(Yvar) <- method.names
  
  # get method pairs of interest
  method.pairs <- matrix(unlist(strsplit(rownames(contrast.matrix), 
                                         "-")), byrow = T, ncol = 2)		
  rownames(method.pairs) <- rownames(contrast.matrix)
  
  result <- mapply(function(i, j) {
    zccc(Ymean[i] - Ymean[j], Yvar[i, i], Yvar[j, j], Yvar[i, j])
  }, i = method.pairs[, 1], j = method.pairs[, 2])
  names(result) <- rownames(contrast.matrix)
  return(result)
}


######################
#
# Function for computing log(TDI) as a function of model parameters
#	for unlinked repeated measurements data
# Arguments: 
#		param = transformed model parameters
#		method.names = levels of the factor method
# 		contrast.matrix = a matrix when multiplied by the vector Y gives
#			the vector of differences of interest for multiple comparisons
# 		prob = probability
#		ulim = upper limit of the interval in which to search for the root 
# Value:
#		A vector of values of log(TDI)
#
ltdi.fun.multi <- function(param, method.names, contrast.matrix, prob, ulim = 500) {
  J <- length(method.names)
  Z <- cbind(rep(1, J), diag(J))
  
  # get the mean and variance of Y vector
  Ymean <- mean.fun(param, method.names)
  sigmasq.b <- exp(param["lsigmasq.b"])
  sigmasq.e <- exp(param[paste0("lsigmasq.e", method.names)])
  if (J == 2) {
    sigmasq.bI <- rep(exp(param["lsigmasq.bI"]), J)
  }
  if (J > 2) {
    sigmasq.bI <- exp(param[paste0("lsigmasq.bI", method.names)])
  }
  Yvar <- Z %*% diag(c(sigmasq.b, sigmasq.bI)) %*% t(Z) + diag(sigmasq.e)
  
  # get the mean and standard deviations of the elements of D vector
  Dmean <- as.numeric(contrast.matrix %*% Ymean)
  Dsd <- sqrt(diag(contrast.matrix %*% Yvar %*% t(contrast.matrix)))
  
  # get log(TDI)
  result <- mapply(ltdi, dmean = Dmean, dsd = Dsd, MoreArgs = list(prob = prob, 
                                                                   ulim = ulim))
  names(result) <- rownames(contrast.matrix)
  return(result)
}


##################
#
# Function for computing log(precisions) of methods as a function of model parameters
#	for unreplicated data as well as unlinked repeated measurements data
# Arguments: 
#		param = transformed model parameters
#		method.names = levels of the factor method
#		include.sigmasq.bI = a logical indicating whether to include 
#			in precision interaction variance 
# Value:
#		A vector of length equal to number of levels of method containing 
#		the log(precisions) [precision = 1/variance]
#
lprec.fun.multi <- function(param, method.names, include.sigmasq.bI = TRUE) {
  J <- length(method.names)
  sigmasq.e <- exp(param[paste0("lsigmasq.e", method.names)])
  
  if (J == 2) {
    if (include.sigmasq.bI) {
      sigmasq.bI <- exp(param["lsigmasq.bI"])
      result <- -log(sigmasq.bI + sigmasq.e)
    } else {
      result <- -log(sigmasq.e)
    }
  }
  
  if (J > 2) {
    if (include.sigmasq.bI) {
      sigmasq.bI <- exp(param[paste0("lsigmasq.bI", method.names)])
      result <- -log(sigmasq.bI + sigmasq.e)
    } else {
      result <- -log(sigmasq.e)
    }
  }
  names(result) <- paste0("lprec", method.names)
  return(result)
}

######################
#
# Function for computing means of methods as a function of model parameters
#	for unreplicated data as well as unlinked repeated measurements data
# Arguments: 
#		param = transformed model parameters
#		method.names = levels of the factor method
# Value:
#		A vector of length equal to number of levels of method containing the means

mean.fun <- function(param, method.names) {
  J <- length(method.names)
  if (J == 2) {
    alpha <- param["alpha"]
  }
  if (J > 2) {
    alpha <- param[paste0("alpha", method.names[2:J])]
  }
  mu.b <- param["mu.b"]
  mu <- c(mu.b, alpha + mu.b)
  names(mu) <- method.names
  return(mu)
}
######################
#
# Function for computing estimate of m-vector valued function of model parameters
# its m x m estimated covariance matrix, and their individual confidence bounds or
# intervals
# Arguments: 
#		param.hat = estimated transformed model parameter vector
#		vcov.param.hat = estimated covariance matrix of param.hat
#		FUN = Name of the m-vector valued parameter function of interest
#		comp.matrix = a matrix that pre-multiplies the vector of parameter
#	 		functions to provide the comparisons of interest
#			(Note: This is may be a diagonal matrix)
# 		alternative =  "two.sided" or "less" (upper bound) or "greater" (lower bound)
# 		level = level of confidence
# Value:
#		A list contanining the estimated parameter function, its covariance matrix,
#			and individual confidence intervals


ci.ind.paramfun <- function(param.hat, vcov.param.hat, FUN, comp.matrix, 
                            alternative, level, ...) {
  require(numDeriv)
  require(multcomp, quietly = T, warn.conflicts = F)
  # estimate and covariance matrix
  fun.hat <- FUN(param.hat, ...)
  fun.grad <- jacobian(FUN, param.hat, ...)
  fun.vcov <- fun.grad %*% vcov.param.hat %*% t(fun.grad)
  
  # univariate confidence intervals or bounds
  
  fun.mcomp <- glht(model = parm(coef = fun.hat, vcov = fun.vcov), linfct = comp.matrix, 
                    alternative = alternative)
  ci.fun <- confint(fun.mcomp, level = level, calpha=univariate_calpha())
  result <- cbind(estimate = coef(fun.mcomp), se = sqrt(diag(vcov(fun.mcomp))), 
                  lcl = ci.fun$confint[, "lwr"], ucl = ci.fun$confint[, "upr"])
  return(list(ci.fun = ci.fun, result = result))
}



######################
#
# Function for computing estimate of m-vector valued function of model parameters
# and its m x m estimated covariance matrix, and their simultaneous confidence bounds
# and intervals
# Arguments: 
#		param.hat = estimated transformed model parameter vector
#		vcov.param.hat = estimated covariance matrix of param.hat
#		FUN = Name of the m-vector valued parameter function of interest
#		comp.matrix = a matrix that pre-multiplies the vector of parameter
#	 		functions to provide the comparisons of interest
#			(Note: This is may be a diagonal matrix)
# 		alternative =  "two.sided" or "less" (upper bound) or "greater" (lower bound)
# 		level = level of confidence
# Value:
#		A list contanining the estimated parameter function, its covariance matrix,
#			and simultaneous confidence intervals


ci.simul.paramfun <- function(param.hat, vcov.param.hat, FUN, comp.matrix, 
                              alternative, level, ...) {
  require(numDeriv)
  require(multcomp, quietly = T, warn.conflicts = F)
  # estimate and covariance matrix
  fun.hat <- FUN(param.hat, ...)
  fun.grad <- jacobian(FUN, param.hat, ...)
  fun.vcov <- fun.grad %*% vcov.param.hat %*% t(fun.grad)
  
  # simultaneous confidence intervals or bounds
  fun.mcomp <- glht(model = parm(coef = fun.hat, vcov = fun.vcov), linfct = comp.matrix, 
                    alternative = alternative)
  #  	ci.fun <- confint(fun.mcomp, level = level)
  # detach the multcomp package to allow using the qmvt function in multcomp_stuff.R
  detach("package:multcomp", unload = T)
  ci.fun <- confint.glht(fun.mcomp, level = level)
  result <- cbind(estimate = coef(fun.mcomp), se = sqrt(diag(vcov(fun.mcomp))), 
                  lcl = ci.fun$confint[, "lwr"], ucl = ci.fun$confint[, "upr"])
  return(list(ci.fun = ci.fun, result = result))
}


######################
#
# Function for computing loglikelihood function as a function of model parameters
# for unlinked repeated measurements data with multiple methods (J >=2)
# Arguments: 
#		param = transformed model parameters
#		method.names = levels of the factor method
# 		dat = data in long format, used for model fitting
# Value:
#		loglikelihood function

loglik.unlinked.data.multi <- function(param, method.names, dat) {
  require(mvtnorm)
  # function for computing loglik for one subject
  loglik.unlinked.one <- function(yy) {
    X <- model.matrix(~-1 + names(yy))
    Z <- cbind(rep(1, length(yy)), X)
    R <- diag(as.numeric(X %*% sigmasq.e))
    return(dmvnorm(yy, mean = (X %*% beta), sigma = (Z %*% G %*% t(Z) + 
                                                       R), log = T))
  }
  
  J <- length(method.names)
  
  if (J == 2) {
    alpha <- param["alpha"]
    mu.b <- param["mu.b"]
    sigmasq.b <- exp(param["lsigmasq.b"])
    sigmasq.bI <- exp(param["lsigmasq.bI"])
    sigmasq.e <- exp(param[paste0("lsigmasq.e", method.names)])
    names(sigmasq.e) <- method.names
    G <- diag(c(sigmasq.b, rep(sigmasq.bI, J)))
  }
  
  if (J > 2) {
    alpha <- param[paste0("alpha", method.names[2:J])]
    mu.b <- param["mu.b"]
    sigmasq.b <- exp(param["lsigmasq.b"])
    sigmasq.bI <- exp(param[paste0("lsigmasq.bI", method.names)])
    sigmasq.e <- exp(param[paste0("lsigmasq.e", method.names)])
    names(sigmasq.e) <- method.names
    G <- diag(c(sigmasq.b, sigmasq.bI))
  }
  
  # parameters
  beta <- c(mu.b, alpha + mu.b)
  names(beta) <- method.names
  
  # get data in wide format
  ylist <- long2wide.unlinked(dat)
  
  # compute loglik for each subject
  loglik <- lapply(ylist, loglik.unlinked.one)
  
  return(sum(as.numeric(loglik)))
}


######################
#
# Function for coverting data in long format, used for model fitting, to wide
# format useful for computing log-likelihood
# Arguments: 
#		dat = data in long format
# Value:
#	data in wide format, one list for each subject
#
long2wide.unlinked <- function(dat) {
  assign.method <- function(obs, method.ind) {
    names(obs) <- method.ind
    return(obs)
  }
  # split the meas and method columns	
  obs <- split(dat$meas, dat$subject)
  method.ind <- split(dat$method, dat$subject)
  return(mapply(assign.method, obs, method.ind, SIMPLIFY = FALSE))
}

##################
#
# Function for fitting a mixed model to unlinked repeated measurements data
#		with multiple (J >=2 ) methods
# Arguments: 
#		dat = a data frame with four columns --- subject (a factor), 
#			  	method (a factor), 
#				rep (replication number, unused) and meas
#		comp.type = a character with two possible values ("all.pairwise" or 
#			"reference") indicating the type of comparisons (default = "reference")
# Value: A list with following components
#		fit = fitted lme object
#		param.hat = ML estimates of transformed parameters
#		var.hat = ML estimates of variance components and variance of difference
#		ydist.hat = parameters of fitted multivariate normal distribution
#		dist.hat = parameters of the fitted distribution of the differences
#		reliability.hat = ML estimates of reliabilities of the two methods
#		method.names = a character vector of names of methods
#		contrast.matrix = the contrast matrix determined by comp.type
#
# Note: (a) This function can be used for J=2 as well.
#		(b) Unlike the version of this function for J=2 case, the levels of the factor
# 			method don't have to be method1, method2, ..., methodJ

lme.unlinked.repmeas.fit.multi <- function(dat, comp.type = "reference") {
  require(nlme)
  method.names.used <- droplevels(dat$method)
  J <- nlevels(method.names.used)
  method.names <- levels(method.names.used)
  
  if (J < 2) {
    stop("the factor method should have at least two levels")
  }
  
  # fit model
  dat <- groupedData(meas ~ method | subject, data = dat, order.groups = FALSE)
  
  if (J == 2) {
    fit <- lme(meas ~ method - 1, data = dat, weights = varIdent(form = ~1 | 
                                                                   method), method = "ML", random = pdBlocked(list(~1, pdIdent(~method - 
                                                                                                                                 1))), control = lmeControl(maxIter = 1000, msMaxIter = 5000, niterEM = 2000, 
                                                                                                                                                            msMaxEval = 1500))
    
    # extract random effect variances
    sigma.hat <- fit$sigma
    mat1 <- diag(as.matrix(fit$modelStruct$reStruct[[1]]) * sigma.hat^2)
    sigmasq.b.hat <- mat1["(Intercept)"]
    sigmasq.bI.hat <- mat1[paste0("method", method.names[1])]
    
    # extract the remaining estimates
    mu.hat <- fixef(fit)[paste0(rep("method", J), method.names)]
    names(mu.hat) <- method.names
    
    sigma.ratio.hat <- coef(fit$modelStruct$varStruct, unconstrained = FALSE, 
                            allCoef = TRUE)[method.names]
    sigmasq.e.hat <- (sigma.hat * sigma.ratio.hat)^2
    
    # cast estimates in the desired form
    mu.b.hat <- mu.hat[method.names[1]]
    alpha.hat <- mu.hat[method.names[2:J]] - mu.hat[method.names[1]] # a (J-1)-vector
    
    param.hat <- c(alpha.hat, mu.b.hat, log(sigmasq.b.hat), log(sigmasq.bI.hat), 
                   log(sigmasq.e.hat))
    names(param.hat) <- c("alpha", "mu.b", "lsigmasq.b", "lsigmasq.bI", 
                          paste0("lsigmasq.e", method.names))
  }
  
  if (J > 2) {
    fit <- lme(meas ~ method - 1, data = dat, weights = varIdent(form = ~1 | 
                                                                   method), method = "ML", random = pdBlocked(list(~1, pdDiag(~method - 
                                                                                                                                1))), control = lmeControl(maxIter = 1000, msMaxIter = 5000, niterEM = 2000, 
                                                                                                                                                           msMaxEval = 1500))
    
    # extract random effect variances
    sigma.hat <- fit$sigma
    mat1 <- diag(as.matrix(fit$modelStruct$reStruct[[1]]) * sigma.hat^2)
    sigmasq.b.hat <- mat1["(Intercept)"]
    sigmasq.bI.hat <- mat1[paste0(rep("method", J), method.names)]
    names(sigmasq.bI.hat) <- method.names
    
    # extract the remaining estimates
    mu.hat <- fixef(fit)[paste0(rep("method", J), method.names)]
    names(mu.hat) <- method.names
    
    sigma.ratio.hat <- coef(fit$modelStruct$varStruct, unconstrained = FALSE, 
                            allCoef = TRUE)[method.names]
    sigmasq.e.hat <- (sigma.hat * sigma.ratio.hat)^2
    
    # cast estimates in the desired form
    mu.b.hat <- mu.hat[method.names[1]]
    alpha.hat <- mu.hat[method.names[2:J]] - mu.hat[method.names[1]] # a (J-1)-vector
    
    param.hat <- c(alpha.hat, mu.b.hat, log(sigmasq.b.hat), log(sigmasq.bI.hat), 
                   log(sigmasq.e.hat))
    names(param.hat) <- c(paste0("alpha", method.names[2:J]), "mu.b", "lsigmasq.b", 
                          paste0("lsigmasq.bI", method.names), paste0("lsigmasq.e", method.names))
  }
  
  # estimated variance parameters
  varpar.hat <- list(sigmasq.b = sigmasq.b.hat, sigmasq.bI = sigmasq.bI.hat, 
                     sigmasq.e = sigmasq.e.hat)
  
  # fitted J-variate distribution 
  Z <- cbind(rep(1, J), diag(J))
  
  if (J == 2) {
    var.hat <- Z %*% diag(c(sigmasq.b.hat, rep(sigmasq.bI.hat, J))) %*% 
      t(Z) + diag(sigmasq.e.hat)
  }
  
  if (J > 2) {
    var.hat <- Z %*% diag(c(sigmasq.b.hat, sigmasq.bI.hat)) %*% t(Z) + 
      diag(sigmasq.e.hat)
  }
  rownames(var.hat) <- method.names
  colnames(var.hat) <- method.names
  
  ydist.hat <- list(ymean = mu.hat, ycov = var.hat, ycor = cov2cor(var.hat), 
                    yvar = diag(var.hat))
  
  # fitted distribution of difference vector of interest 
  contrast.matrix <- make.contrast.matrix(method.names, comp.type)
  dmean.hat <- as.numeric(contrast.matrix %*% mu.hat)
  names(dmean.hat) <- rownames(contrast.matrix)
  dvar.hat <- contrast.matrix %*% var.hat %*% t(contrast.matrix)
  ddist.hat <- list(dmean = dmean.hat, dcov = dvar.hat, dvar = diag(dvar.hat))
  
  # reliability
  rlbty.hat <- (sigmasq.b.hat + sigmasq.bI.hat)/diag(var.hat)
  names(rlbty.hat) <- method.names
  
  return(list(fit = fit, method.names = method.names, param.hat = param.hat, 
              varpar.hat = varpar.hat, ydist.hat = ydist.hat, contrast.matrix = contrast.matrix, 
              ddist.hat = ddist.hat, reliability.hat = rlbty.hat))
}

# ######################
# Function for creating a contrast matrix for multiple comparisons of methods
# Arguments: 
#		method.names = a character vector containing levels of the method factor
#		comp.type = a character with two possible values ("all.pairwise" or 
#			"reference") indicating the type of comparisons [first level of method
#		factor serves as the reference]
# Value: A contrast matrix with nmethods = length(method.names) columns 
#		and number of rows equal to
#		either (nmethods*(nmethods-1))/2 ["all.pairwise"] or (nmethods-1) ["reference"]
#


make.contrast.matrix <- function(method.names, comp.type) {
  nmethods <- length(method.names)
  
  if (!(comp.type %in% c("all.pairwise", "reference"))) {
    stop("either missing or invalid comparison type in make.contrast.matrix")
  }
  
  if (comp.type == "all.pairwise") {
    Cmat <- matrix(0, nrow = (nmethods * (nmethods - 1)/2), ncol = nmethods)
    row.label <- character(nrow(Cmat))
    row.number <- 1
    for (i in 1:(nmethods - 1)) {
      for (j in (i + 1):nmethods) {
        Cmat[row.number, c(i, j)] <- c(-1, 1)
        row.label[row.number] <- paste0(method.names[j], "-", method.names[i])
        row.number <- row.number + 1
      }
    }
    rownames(Cmat) <- row.label
  }
  
  if (comp.type == "reference") {
    Cmat <- cbind(rep(-1, nmethods - 1), diag(nmethods - 1))
    rownames(Cmat) <- paste0(method.names[2:nmethods], "-", method.names[1])
  }
  return(Cmat)
}

# ######################
# Function for creating a matrix of plots with scatterplots in lower diagonals
# and B-A plots in upper diagonals
# Arguments: 
#		ymat = a n x J matrix containing one measurement per method on every subject, 
#			where n = number of subjects and J = number of methods
# Value: A plot and the ranges of observations and their differences used to 
#			to set up axes


plot.sba.matrix <- function(ymat) {
  J <- ncol(ymat)
  orig.yrange <- range(ymat) 
  yrange <- orig.yrange + c(-1,1)*diff(orig.yrange)*0.05
  orig.drange <- range(sapply(1:(J - 1), function(i) {
    ymat[, (i + 1):J] - ymat[, i] })) 
  drange <- orig.drange + c(-1,1)*diff(orig.drange)*0.05 
  pairs(ymat, upper.panel = function(x, y, ...) {
    par(usr = c(yrange, drange))
    points((x + y)/2,  x-y)
    abline(h = 0)
  }, lower.panel = function(x, y, ...) {
    par(usr = rep(yrange, 2))
    points(x, y)
    abline(a = 0, b = 1)
  }, gap = 0, xaxt = "n", yaxt = "n", labels = colnames(ymat), cex.labels = 1.5)
  return(list(yrange = orig.yrange, drange = orig.drange))
}
# ######################


####################################################################################
# chpater 6 #
####################################################################################

##################################

#
# Function for computing estimates, their SEs and confidence intervals
#		 for various parameters and functions for heteroscedastic paired
#		 measurements data
# Arguments: 
#		ymat = a matrix with n rows and two columns --- method1 and method2
#		param.hat = vector of ML estimates of transformed parameters
#		vtilde = a vector of length nrow(dat) containing the subject-specific 
#				value of the variance covariate 
#		varfunc = a character with one of two possible values---
#				"power" or "exponential"--- indicating the common variance function
#				model for the two methods
#		vgrid = a grid of values on which to evaluate the fitted variance 
#					functions and the paramerters that involve them
# Value: 
#		standard errors of estimates of parameters and their confidence intervals
#		estimate, se and ci of log(var ratio), and estimate and ucb of var ratio,
#		estimate, se and ucb of log(TDI), and estimate and ucb of TDI
#		estimate, se and lcb of z(CCC), and estimate and lcb of CCC
#

conf.measures.hetero.bvn <- function(ymat, param.hat, vtilde, varfunc, vgrid, 
                                     conf.level = 0.95, prob = 0.9, ulim=500) {
  require(numDeriv)
  
  if (!(varfunc %in% c("power", "exponential"))) {
    stop(paste(varfunc, "is not a valid option in conf.measures.hetero.bvn"))
  }
  
  z.twosided <- qnorm((1 + conf.level)/2)
  z.onesided <- qnorm(conf.level)
  
  # get Hessian and invert it
  hess <- -hessian(loglik.hetero.paired.data, param.hat, ymat = ymat, vtilde = vtilde, 
                   varfunc = varfunc)
  inv.hess <- solve(hess)
  rownames(inv.hess) <- names(param.hat)
  colnames(inv.hess) <- names(param.hat)
  
  
  # SEs of param.hat and confidence intevals
  se.param.hat <- sqrt(diag(inv.hess))
  ci.param <- cbind(param.hat - z.twosided * se.param.hat, param.hat + z.twosided * 
                      se.param.hat)
  results.param <- cbind(param.hat, se.param.hat, ci.param)
  rownames(results.param) <- names(param.hat)
  colnames(results.param) <- c("estimate", "se", "lcl", "ucl")
  
  # estimate and confidence interval for mean difference 
  meandiff.hat <- param.hat["mu.2"] - param.hat["mu.1"]
  se.meandiff.hat <- sqrt(inv.hess["mu.2", "mu.2"] + inv.hess["mu.1", "mu.1"] - 
                            2 * inv.hess["mu.2", "mu.1"])
  ci.meandiff <- meandiff.hat + c(-1, 1) * z.twosided * se.meandiff.hat
  meandiff <- c(meandiff.hat, se.meandiff.hat, ci.meandiff)
  names(meandiff) <- c("est", "se", "lcl", "ucl")
  
  
  # estimate and poitwise confidence interval for precision ratio
  lvratio.hat <- sapply(vgrid, lvratio.hetero.bvn.fun, param = param.hat, 
                        varfunc = varfunc)
  lvratio.grad <- sapply(vgrid, function(vv) {
    grad(lvratio.hetero.bvn.fun, param.hat, vval = vv, varfunc = varfunc)
  })
  se.lvratio.hat <- apply(lvratio.grad, MAR = 2, function(gg) {
    sqrt(drop(t(gg) %*% inv.hess %*% gg))
  })
  lcl.lvratio <- lvratio.hat - z.twosided * se.lvratio.hat
  ucl.lvratio <- lvratio.hat + z.twosided * se.lvratio.hat
  
  results.vratio <- cbind(vgrid, lvratio.hat, se.lvratio.hat, lcl.lvratio, 
                          ucl.lvratio, exp(lvratio.hat), exp(lcl.lvratio), exp(ucl.lvratio))
  rownames(results.vratio) <- NULL
  colnames(results.vratio) <- c("vgrid", "lest", "lest.se", "lcl.lpar", "ucl.lpar", 
                                "est", "lcl.par", "ucl.par")
  
  # estimate and pointwise upper confidence bound of TDI
  ltdi.hat <- sapply(vgrid, ltdi.hetero.bvn.fun, param = param.hat, varfunc = varfunc, 
                     prob = prob, ulim=ulim)
  ltdi.grad <- sapply(vgrid, function(vv) {
    grad(ltdi.hetero.bvn.fun, param.hat, vval = vv, varfunc = varfunc, 
         prob = prob, ulim=ulim)
  })
  se.ltdi.hat <- apply(ltdi.grad, MAR = 2, function(gg) {
    sqrt(drop(t(gg) %*% inv.hess %*% gg))
  })
  ucb.ltdi <- ltdi.hat + z.onesided * se.ltdi.hat
  
  results.tdi <- cbind(vgrid, ltdi.hat, se.ltdi.hat, ucb.ltdi, exp(ltdi.hat), 
                       exp(ucb.ltdi))
  rownames(results.tdi) <- NULL
  colnames(results.tdi) <- c("vgrid", "ltdi.hat", "se.ltdi.hat", "ucb.ltdi", 
                             "tdi.hat", "ucb.tdi")
  
  # estimate and pointwise lower confidence bound of CCC
  zccc.hat <- sapply(vgrid, zccc.hetero.bvn.fun, param = param.hat, varfunc = varfunc)
  zccc.grad <- sapply(vgrid, function(vv) {
    grad(zccc.hetero.bvn.fun, param.hat, vval = vv, varfunc = varfunc)
  })
  se.zccc.hat <- apply(zccc.grad, MAR = 2, function(gg) {
    sqrt(drop(t(gg) %*% inv.hess %*% gg))
  })
  lcb.zccc <- zccc.hat - z.onesided * se.zccc.hat
  
  results.ccc <- cbind(vgrid, zccc.hat, se.zccc.hat, lcb.zccc, tanh(zccc.hat), 
                       tanh(lcb.zccc))
  rownames(results.ccc) <- NULL
  colnames(results.ccc) <- c("vgrid", "zccc.hat", "se.zccc.hat", "lcb.zccc", 
                             "ccc.hat", "lcb.ccc")
  
  
  results <- list(param = results.param, meandiff=meandiff, vratio = results.vratio, tdi = results.tdi, 
                  ccc = results.ccc)
  return(results)
}
#######################
#
# Function for computing z(CCC) as a function of model parameters and variance
# 	covariate for paired measurements data
# Arguments: 
#		param = transformed model parameters
#		vval = value of variance covariate at which to evaluate log (TDI)
#		varfunc = a character with one of two possible values---
#				"power" or "exponential"--- indicating the common variance function
#				model for the two methods
#
#	Value:
#		z(CCC) at vval
#
zccc.hetero.bvn.fun <- function(param, vval, varfunc) {
  # extract the parameters
  mu.1 <- param["mu.1"]
  mu.2 <- param["mu.2"]	
  sigmasq.1 <- exp(param["lsigmasq.1"])
  sigmasq.2 <- exp(param["lsigmasq.2"])
  delta.1 <- param["delta.1"]
  delta.2 <- param["delta.2"]
  rho <- tanh(param["zrho"])
  
  # variances 	
  if (varfunc == "power") {
    var.1 <- sigmasq.1 * (abs(vval)^(2 * delta.1))
    var.2 <- sigmasq.2 * (abs(vval)^(2 * delta.2))
  }
  if (varfunc == "exponential") {
    var.1 <- sigmasq.1 * exp(2 * delta.1 * vval)
    var.2 <- sigmasq.2 * exp(2 * delta.2 * vval)
  }
  
  cov.12 <- rho * sqrt(var.1 * var.2)
  cov.mat <- matrix(c(var.1, rep(cov.12, 2), var.2), byrow = T, ncol = 2)
  
  result <- zccc(mu.2-mu.1, cov.mat[1, 1], cov.mat[2, 2], cov.mat[1, 2])
  return(unname(result))
}


##################
#
# Function for computing log(TDI) as a function of model parameters and variance
# 	covariate for paired measurements data
# Arguments: 
#		param = transformed model parameters
#		vval = value of variance covariate at which to evaluate log (TDI)
#		varfunc = a character with one of two possible values---
#				"power" or "exponential"--- indicating the common variance function
#				model for the two methods
# 		prob = probability
#		ulim = upper limit of the interval that is expectd to contain the root
#
#	Value:
#		log(TDI) at vval
#
ltdi.hetero.bvn.fun <- function(param, vval, varfunc, prob, ulim=500) {
  # extract the parameters
  mu.1 <- param["mu.1"]
  mu.2 <- param["mu.2"]	
  sigmasq.1 <- exp(param["lsigmasq.1"])
  sigmasq.2 <- exp(param["lsigmasq.2"])
  delta.1 <- param["delta.1"]
  delta.2 <- param["delta.2"]
  rho <- tanh(param["zrho"])
  
  # variances 	
  if (varfunc == "power") {
    var.1 <- sigmasq.1 * (abs(vval)^(2 * delta.1))
    var.2 <- sigmasq.2 * (abs(vval)^(2 * delta.2))
  }
  if (varfunc == "exponential") {
    var.1 <- sigmasq.1 * exp(2 * delta.1 * vval)
    var.2 <- sigmasq.2 * exp(2 * delta.2 * vval)
  }
  
  var.D <- var.1 + var.2 - 2 * rho * sqrt(var.1 * var.2)
  
  result <- ltdi(mu.2-mu.1, sqrt(var.D), prob, ulim=ulim)
  return(result)
}


##################
#
# Function for computing log(var2/var1) as a function of model parameters and variance
# 	covariate for paired measurements data
# Arguments: 
#		param = transformed model parameters
#		vval = value of variance covariate at which to evaluate log (lambda)
#		varfunc = a character with one of two possible values---
#				"power" or "exponential"--- indicating the common variance function
#				model for the two methods
#
#	Value:
#		log(var2/var1) at vval
#
lvratio.hetero.bvn.fun <- function(param, vval, varfunc) {
  # extract the parameters
  sigmasq.1 <- exp(param["lsigmasq.1"])
  sigmasq.2 <- exp(param["lsigmasq.2"])
  delta.1 <- param["delta.1"]
  delta.2 <- param["delta.2"]
  
  # variances 	
  if (varfunc == "power") {
    var.1 <- sigmasq.1 * (abs(vval)^(2 * delta.1))
    var.2 <- sigmasq.2 * (abs(vval)^(2 * delta.2))
  }
  if (varfunc == "exponential") {
    var.1 <- sigmasq.1 * exp(2 * delta.1 * vval)
    var.2 <- sigmasq.2 * exp(2 * delta.2 * vval)
  }
  result <- log(var.2) - log(var.1)
  return(unname(result))
}


##################
#
# Function for fitting a BVN model to the paired data
# Argument: 
# 		ymat = data matrix with two columns: method1 and method2
#		vtilde = variance covariate, a vector with length nrow(ymat)
#		varfunc = a character with one of two possible values---
#				"power" or "exponential"--- indicating the common variance function
#				model for the two methods
#		vgrid = a grid of values of variance covariate over which to evaluate 
#				its various functions 
# Value: A list with following components
#		param.hat = ML estimates of transformed parameters
#		bvn.hat = parameters of fitted bivariate normal distribution
#		and several others							
#
bvn.hetero.fit <- function(ymat, vtilde, varfunc, vgrid, llim.delta = -2, ulim.delta = 2) {
  require(numDeriv)
  require(mvtnorm)
  
  if (!(varfunc %in% c("power", "exponential"))) {
    stop(paste(varfunc, "is not a valid option in bvn.hetero.fit"))
  }
  
  n <- nrow(ymat)
  
  # homoscedastic fit to generate initial values
  fit.homo <- bvn.fit(ymat)
  param.init <- c(fit.homo$param.hat, rep(0, 2))
  names(param.init) <- c(names(fit.homo$param.hat), "delta.1", "delta.2")
  
  # maximize the likelihood
  fit.hetero <- nlminb(param.init, function(param) {
    -loglik.hetero.paired.data(param, ymat, vtilde, varfunc)
  }, lower = c(rep(-Inf, 5), rep(llim.delta, 2)), upper = c(rep(Inf, 5), 
                                                            rep(ulim.delta, 2)), control = list(eval.max = 1000, iter.max = 750))
  
  if (fit.hetero$convergence != 0) {
    stop(paste("nlminb in bvn.hetero.fit failed to converge while maximizing likelihood"))
  }
  
  param.hat <- fit.hetero$par
  loglik <- -fit.hetero$objective
  aic <- -2 * loglik + 2 * length(param.hat)
  bic <- -2 * loglik + log(2 * n) * length(param.hat)
  
  # compare with homoscedasitic model
  lrt.stat <- 2 * (loglik - fit.homo$loglik)
  pval <- 1 - pchisq(lrt.stat, df = 2)
  
  anova.results <- cbind(c(fit.homo$loglik, loglik), c(fit.homo$aic, aic), 
                         c(fit.homo$bic, bic), c(NA, pval))
  rownames(anova.results) <- c("homoscedastic", "heteroscedastic")
  colnames(anova.results) <- c("loglik", "aic", "bic", "pval")
  
  
  # extract the estimates
  mu.1.hat <- unname(param.hat["mu.1"])
  mu.2.hat <- unname(param.hat["mu.2"])
  sigmasq.1.hat <- unname(exp(param.hat["lsigmasq.1"]))
  sigmasq.2.hat <- unname(exp(param.hat["lsigmasq.2"]))
  rho.hat <- unname(tanh(param.hat["zrho"]))
  delta.1.hat <- unname(param.hat["delta.1"])
  delta.2.hat <- unname(param.hat["delta.2"])
  
  
  # fitted variances	
  if (varfunc == "power") {
    fitted.var.1.hat <- sigmasq.1.hat * (abs(vtilde)^(2 * delta.1.hat))
    fitted.var.2.hat <- sigmasq.2.hat * (abs(vtilde)^(2 * delta.2.hat))
  }
  
  if (varfunc == "exponential") {
    fitted.var.1.hat <- sigmasq.1.hat * exp(2 * delta.1.hat * vtilde)
    fitted.var.2.hat <- sigmasq.2.hat * exp(2 * delta.2.hat * vtilde)
  }
  
  fitted.var.D.hat <- fitted.var.1.hat + fitted.var.2.hat - 2 * rho.hat * sqrt(fitted.var.1.hat * fitted.var.2.hat)
  
  fitted.var <- cbind(method1=fitted.var.1.hat, method2=fitted.var.2.hat, D=fitted.var.D.hat)
  
  # fitted bivariate distribution on vgrid
  if (varfunc == "power") {
    var.1.hat <- sigmasq.1.hat * (abs(vgrid)^(2 * delta.1.hat))
    var.2.hat <- sigmasq.2.hat * (abs(vgrid)^(2 * delta.2.hat))
  }
  
  if (varfunc == "exponential") {
    var.1.hat <- sigmasq.1.hat * exp(2 * delta.1.hat * vgrid)
    var.2.hat <- sigmasq.2.hat * exp(2 * delta.2.hat * vgrid)
  }
  
  bvn.hat <- list(mean1 = mu.1.hat, mean2 = mu.2.hat, var1 = var.1.hat, var2 = var.2.hat, 
                  corr = rho.hat)
  
  # fitted distribution of the difference
  
  mean.D.hat <- mu.2.hat - mu.1.hat
  var.D.hat <- var.1.hat + var.2.hat - 2 * rho.hat * sqrt(var.1.hat * var.2.hat)
  D.hat <- list(mean = mean.D.hat, var = var.D.hat)
  
  # limits of agreement
  loa <- matrix(0, nrow = length(vgrid), ncol = 2)
  colnames(loa) <- c("lower", "upper")
  loa[, 1] <- mean.D.hat - 1.96 * sqrt(var.D.hat)
  loa[, 2] <- mean.D.hat + 1.96 * sqrt(var.D.hat)
  
  return(list(fit.homo = fit.homo, vtilde = vtilde, varfunc = varfunc, fit.hetero = fit.hetero, 
              param.hat = param.hat, anova.results = anova.results, 
              fitted.var=fitted.var,
              vgrid = vgrid, bvn.hat = bvn.hat, D.hat = D.hat, loa = loa))
  
}

######################
#
# Function for computing log(L) for heteroscedastic paired data 
# Arguments: 
#		param = vector of transformed model parameters
# 		ymat = data matrix with two columns: method1 and method2
#		vtilde = variance covariate, a vector with length nrow(ymat)
#		varfunc = a character with one of two possible values---
#				"power" or "exponential"--- indicating the common variance function
#				model for the two methods
# Value: The loglikelihood function
#


loglik.hetero.paired.data <- function(param, ymat, vtilde, varfunc) {
  require(mvtnorm)
  n <- nrow(ymat)
  
  # function for computing loglik for one subject
  loglik.hetero.paired.one <- function(param, yy, vv, varfunc) {
    # extract the parameters
    mu.1 <- param["mu.1"]
    mu.2 <- param["mu.2"]
    sigmasq.1 <- exp(param["lsigmasq.1"])
    sigmasq.2 <- exp(param["lsigmasq.2"])
    rho <- tanh(param["zrho"])
    delta.1 <- param["delta.1"]
    delta.2 <- param["delta.2"]
    
    # variances	
    if (varfunc == "power") {
      var.1 <- sigmasq.1 * (abs(vv)^(2 * delta.1))
      var.2 <- sigmasq.2 * (abs(vv)^(2 * delta.2))
    }
    if (varfunc == "exponential") {
      var.1 <- sigmasq.1 * exp(2 * delta.1 * vv)
      var.2 <- sigmasq.2 * exp(2 * delta.2 * vv)
    }
    
    mean.vec <- c(mu.1, mu.2)
    cov.mat <- matrix(c(var.1, rho * sqrt(var.1 * var.2), rho * sqrt(var.1 * 
                                                                       var.2), var.2), byrow = T, ncol = 2)
    
    return(dmvnorm(yy, mean = mean.vec, sigma = cov.mat, log = T))
  }
  
  if (!(varfunc %in% c("power", "exponential"))) {
    stop(paste(varfunc, "is not a valid option in loglik.hetero.paired.data"))
  }
  
  # get log-likelihood
  loglik <- sapply(1:n, function(ind) {
    loglik.hetero.paired.one(param, ymat[ind, ], vtilde[ind], varfunc = varfunc)
  })
  
  return(sum(loglik))
}

##################
#
# Function for computing log(lambda) as a function of model parameters and variance
# 	covariate
# Arguments: 
#		param = transformed model parameters
#		vval = value of variance covariate at which to evaluate log (lambda)
#		varfunc = a character with one of two possible values---
#				"power" or "exponential"--- indicating the common variance function
#				model for the two methods
#
#	Value:
#		log(lambda) at vval
#
llambda.hetero.repmeas.fun <- function(param, vval, varfunc) {
  # extract the parameters
  sigmasq.e1 <- exp(param["lsigmasq.e1"])
  sigmasq.e2 <- exp(param["lsigmasq.e2"])
  delta.1 <- param["delta.1"]
  delta.2 <- param["delta.2"]
  
  # error variance 	
  if (varfunc == "power") {
    err.var.1 <- sigmasq.e1 * (abs(vval)^(2 * delta.1))
    err.var.2 <- sigmasq.e2 * (abs(vval)^(2 * delta.2))
  }
  if (varfunc == "exponential") {
    err.var.1 <- sigmasq.e1 * exp(2 * delta.1 * vval)
    err.var.2 <- sigmasq.e2 * exp(2 * delta.2 * vval)
  }
  result <- log(err.var.1) - log(err.var.2)
  return(result)
}
##################
#
# Function for computing log(TDI) as a function of model parameters and variance
# 	covariate
# Arguments: 
#		param = transformed model parameters
#		vval = value of variance covariate at which to evaluate log (TDI)
#		varfunc = a character with one of two possible values---
#				"power" or "exponential"--- indicating the common variance function
#				model for the two methods
# 		prob = probability
#
#	Value:
#		log(TDI) at vval
#
ltdi.hetero.repmeas.fun <- function(param, vval, varfunc, prob) {
  # extract the parameters
  alpha <- param["alpha"]
  sigmasq.bI <- exp(param["lsigmasq.bI"])
  sigmasq.e1 <- exp(param["lsigmasq.e1"])
  sigmasq.e2 <- exp(param["lsigmasq.e2"])
  delta.1 <- param["delta.1"]
  delta.2 <- param["delta.2"]
  
  # error variance 		
  if (varfunc == "power") {
    err.var.1 <- sigmasq.e1 * (abs(vval)^(2 * delta.1))
    err.var.2 <- sigmasq.e2 * (abs(vval)^(2 * delta.2))
  }
  if (varfunc == "exponential") {
    err.var.1 <- sigmasq.e1 * exp(2 * delta.1 * vval)
    err.var.2 <- sigmasq.e2 * exp(2 * delta.2 * vval)
  }
  result <- ltdi(alpha, sqrt(2 * sigmasq.bI + err.var.1 + err.var.2), prob)
  return(result)
}
##################
#
# Function for computing z(CCC) as a function of model parameters and variance
# 	covariate
# Arguments: 
#		param = transformed model parameters
#		vval = value of variance covariate at which to evaluate log (TDI)
#		varfunc = a character with one of two possible values---
#				"power" or "exponential"--- indicating the common variance function
#				model for the two methods
#
#	Value:
#		z(CCC) at vval
#
zccc.hetero.repmeas.fun <- function(param, vval, varfunc) {
  # extract the parameters
  alpha <- param["alpha"]
  sigmasq.b <- exp(param["lsigmasq.b"])
  sigmasq.bI <- exp(param["lsigmasq.bI"])
  sigmasq.e1 <- exp(param["lsigmasq.e1"])
  sigmasq.e2 <- exp(param["lsigmasq.e2"])
  delta.1 <- param["delta.1"]
  delta.2 <- param["delta.2"]
  
  # error variance 		
  if (varfunc == "power") {
    err.var.1 <- sigmasq.e1 * (abs(vval)^(2 * delta.1))
    err.var.2 <- sigmasq.e2 * (abs(vval)^(2 * delta.2))
  }
  if (varfunc == "exponential") {
    err.var.1 <- sigmasq.e1 * exp(2 * delta.1 * vval)
    err.var.2 <- sigmasq.e2 * exp(2 * delta.2 * vval)
  }
  cov.mat <- matrix(c(sigmasq.b + sigmasq.bI + err.var.1, sigmasq.b, sigmasq.b, 
                      sigmasq.b + sigmasq.bI + err.var.2), byrow = T, ncol = 2)
  result <- zccc(alpha, cov.mat[1, 1], cov.mat[2, 2], cov.mat[1, 2])
  return(result)
}
##################
#
# Function for computing log(TDI.1) as a function of model parameters and variance
# 	covariate
# Arguments: 
#		param = transformed model parameters
#		vval = value of variance covariate at which to evaluate log (TDI.1)
#		varfunc = a character with one of two possible values---
#				"power" or "exponential"--- indicating the common variance function
#				model for the two methods
# 		q.chisq = percentile of a chisquare (1) distribution that corresponds to probability
#
#	Value:
#		log(TDI.1) at vval
#
ltdi.1.hetero.repmeas.fun <- function(param, vval, varfunc, q.chisq) {
  # extract the parameters
  sigmasq.e1 <- exp(param["lsigmasq.e1"])
  delta.1 <- param["delta.1"]
  
  # error variance 		
  if (varfunc == "power") {
    err.var.1 <- sigmasq.e1 * (abs(vval)^(2 * delta.1))
  }
  if (varfunc == "exponential") {
    err.var.1 <- sigmasq.e1 * exp(2 * delta.1 * vval)
  }
  result <- log(sqrt(2 * err.var.1 * q.chisq))
  return(result)
}
##################
#
# Function for computing log(TDI.2) as a function of model parameters and variance
# 	covariate
# Arguments: 
#		param = transformed model parameters
#		vval = value of variance covariate at which to evaluate log (TDI.2)
#		varfunc = a character with one of two possible values---
#				"power" or "exponential"--- indicating the common variance function
#				model for the two methods
# 		q.chisq = percentile of a chisquare (1) distribution that corresponds to probability
#
#	Value:
#		log(TDI.2) at vval
#
ltdi.2.hetero.repmeas.fun <- function(param, vval, varfunc, q.chisq) {
  # extract the parameters
  sigmasq.e2 <- exp(param["lsigmasq.e2"])
  delta.2 <- param["delta.2"]
  
  # error variance 		
  if (varfunc == "power") {
    err.var.2 <- sigmasq.e2 * (abs(vval)^(2 * delta.2))
  }
  if (varfunc == "exponential") {
    err.var.2 <- sigmasq.e2 * exp(2 * delta.2 * vval)
  }
  result <- log(sqrt(2 * err.var.2 * q.chisq))
  return(result)
}

##################
#
# Function for computing z(CCC.1) as a function of model parameters and variance
# 	covariate
# Arguments: 
#		param = transformed model parameters
#		vval = value of variance covariate at which to evaluate z (CCC.1)
#		varfunc = a character with one of two possible values---
#				"power" or "exponential"--- indicating the common variance function
#				model for the two methods
#	Value:
#		z(CCC.1) at vval
#
zccc.1.hetero.repmeas.fun <- function(param, vval, varfunc) {
  # extract the parameters
  sigmasq.b <- exp(param["lsigmasq.b"])
  sigmasq.bI <- exp(param["lsigmasq.bI"])
  sigmasq.e1 <- exp(param["lsigmasq.e1"])
  delta.1 <- param["delta.1"]
  
  # error variance 		
  if (varfunc == "power") {
    err.var.1 <- sigmasq.e1 * (abs(vval)^(2 * delta.1))
  }
  if (varfunc == "exponential") {
    err.var.1 <- sigmasq.e1 * exp(2 * delta.1 * vval)
  }
  cov.mat <- matrix(c(sigmasq.b + sigmasq.bI + err.var.1, sigmasq.b + sigmasq.bI, 
                      sigmasq.b + sigmasq.bI, sigmasq.b + sigmasq.bI + err.var.1), byrow = T, 
                    ncol = 2)
  result <- zccc(0, cov.mat[1, 1], cov.mat[2, 2], cov.mat[1, 2])
  return(result)
}

##################
#
# Function for computing z(CCC.2) as a function of model parameters and variance
# 	covariate
# Arguments: 
#		param = transformed model parameters
#		vval = value of variance covariate at which to evaluate z (CCC.2)
#		varfunc = a character with one of two possible values---
#				"power" or "exponential"--- indicating the common variance function
#				model for the two methods
#	Value:
#		z(CCC.2) at vval
#
zccc.2.hetero.repmeas.fun <- function(param, vval, varfunc) {
  # extract the parameters
  sigmasq.b <- exp(param["lsigmasq.b"])
  sigmasq.bI <- exp(param["lsigmasq.bI"])
  sigmasq.e2 <- exp(param["lsigmasq.e2"])
  delta.2 <- param["delta.2"]
  
  # error variance 		
  if (varfunc == "power") {
    err.var.2 <- sigmasq.e2 * (abs(vval)^(2 * delta.2))
  }
  if (varfunc == "exponential") {
    err.var.2 <- sigmasq.e2 * exp(2 * delta.2 * vval)
  }
  cov.mat <- matrix(c(sigmasq.b + sigmasq.bI + err.var.2, sigmasq.b + sigmasq.bI, 
                      sigmasq.b + sigmasq.bI, sigmasq.b + sigmasq.bI + err.var.2), byrow = T, 
                    ncol = 2)
  result <- zccc(0, cov.mat[1, 1], cov.mat[2, 2], cov.mat[1, 2])
  return(result)
}

##################################

#
# Function for computing estimates, their SEs and confidence intervals
#		 for various parameters and functions for heteroscedastic unlinked 
# 			repeated measurements data
# Arguments: 
#		dat = a data frame with four columns --- subject (a factor), 
#			  method (a factor with two levels method1 and method2), 
#			rep (replication number, unused) and meas
#		param.hat = vector of ML estimates of transformed parameters
#		vtilde = a vector of length nrow(dat) containing the subject-specific 
#				value of the variance covariate 
#		varfunc = a character with one of two possible values---
#				"power" or "exponential"--- indicating the common variance function
#				model for the two methods
#		vgrid = a grid of values on which to evaluate the fitted error variance 
#					functions and the paramerters that involve them
# Value: 
#		standard errors of estimates of parameters and their confidence intervals
#		estimate, se and ucb of log(TDI), and estimate and ucb of TDI
#		estimate, se and lcb of z(CCC), and estimate and lcb of CCC
#		estimates, se and ucb of log(TDI.1) and log(TDI.2), and estimates and ucb for 
#				TDI.1 and TDI.2
#		estimates, se and lcb for z(CCC.1) and z(CCC.2), and estimates and lcb for 
#				CCC.1 and CCC.2
#
##################################


conf.measures.hetero.repmeas <- function(dat, param.hat, vtilde, varfunc, vgrid, 
                                         conf.level = 0.95, prob = 0.9) {
  require(numDeriv)
  
  if (!(varfunc %in% c("power", "exponential"))) {
    stop(paste(varfunc, "is not a valid option in conf.measures.hetero.repmeas"))
  }
  
  z.twosided <- qnorm((1 + conf.level)/2)
  z.onesided <- qnorm(conf.level)
  q.chisq <- qchisq(prob, 1)
  
  # get Hessian and invert it
  hess <- -hessian(loglik.hetero.repmeas.data, param.hat, dat = dat, vtilde = vtilde, 
                   varfunc = varfunc)
  inv.hess <- solve(hess)
  rownames(inv.hess) <- names(param.hat)
  colnames(inv.hess) <- names(param.hat)
  
  
  # SEs of param.hat and confidence intevals
  se.param.hat <- sqrt(diag(inv.hess))
  ci.param <- cbind(param.hat - z.twosided * se.param.hat, param.hat + z.twosided * 
                      se.param.hat)
  results.param <- cbind(param.hat, se.param.hat, ci.param)
  rownames(results.param) <- names(param.hat)
  colnames(results.param) <- c("estimate", "se", "lcl", "ucl")
  
  # estimate and poitwise confidence interval for precision ratio
  llambda.hat <- sapply(vgrid, llambda.hetero.repmeas.fun, param = param.hat, 
                        varfunc = varfunc)
  llambda.grad <- sapply(vgrid, function(vv) {
    grad(llambda.hetero.repmeas.fun, param.hat, vval = vv, varfunc = varfunc)
  })
  se.llambda.hat <- apply(llambda.grad, MAR = 2, function(gg) {
    sqrt(drop(t(gg) %*% inv.hess %*% gg))
  })
  lcl.llambda <- llambda.hat - z.twosided * se.llambda.hat
  ucl.llambda <- llambda.hat + z.twosided * se.llambda.hat
  
  results.lambda <- cbind(vgrid, llambda.hat, se.llambda.hat, lcl.llambda, 
                          ucl.llambda, exp(llambda.hat), exp(lcl.llambda), exp(ucl.llambda))
  rownames(results.lambda) <- NULL
  colnames(results.lambda) <- c("vgrid", "lest", "lest.se", "lcl.lpar", "ucl.lpar", 
                                "est", "lcl.par", "ucl.par")
  
  # estimate and pointwise upper confidence bound of TDI
  ltdi.hat <- sapply(vgrid, ltdi.hetero.repmeas.fun, param = param.hat, varfunc = varfunc, 
                     prob = prob)
  ltdi.grad <- sapply(vgrid, function(vv) {
    grad(ltdi.hetero.repmeas.fun, param.hat, vval = vv, varfunc = varfunc, 
         prob = prob)
  })
  se.ltdi.hat <- apply(ltdi.grad, MAR = 2, function(gg) {
    sqrt(drop(t(gg) %*% inv.hess %*% gg))
  })
  ucb.ltdi <- ltdi.hat + z.onesided * se.ltdi.hat
  
  results.tdi <- cbind(vgrid, ltdi.hat, se.ltdi.hat, ucb.ltdi, exp(ltdi.hat), 
                       exp(ucb.ltdi))
  rownames(results.tdi) <- NULL
  colnames(results.tdi) <- c("vgrid", "ltdi.hat", "se.ltdi.hat", "ucb.ltdi", 
                             "tdi.hat", "ucb.tdi")
  
  # estimate and pointwise lower confidence bound of CCC
  zccc.hat <- sapply(vgrid, zccc.hetero.repmeas.fun, param = param.hat, varfunc = varfunc)
  zccc.grad <- sapply(vgrid, function(vv) {
    grad(zccc.hetero.repmeas.fun, param.hat, vval = vv, varfunc = varfunc)
  })
  se.zccc.hat <- apply(zccc.grad, MAR = 2, function(gg) {
    sqrt(drop(t(gg) %*% inv.hess %*% gg))
  })
  lcb.zccc <- zccc.hat - z.onesided * se.zccc.hat
  
  results.ccc <- cbind(vgrid, zccc.hat, se.zccc.hat, lcb.zccc, tanh(zccc.hat), 
                       tanh(lcb.zccc))
  rownames(results.ccc) <- NULL
  colnames(results.ccc) <- c("vgrid", "zccc.hat", "se.zccc.hat", "lcb.zccc", 
                             "ccc.hat", "lcb.ccc")
  
  # evaluation of repeatability
  # tdi
  # method 1
  ltdi.1.hat <- sapply(vgrid, ltdi.1.hetero.repmeas.fun, param = param.hat, 
                       varfunc = varfunc, q.chisq = q.chisq)
  ltdi.1.grad <- sapply(vgrid, function(vv) {
    grad(ltdi.1.hetero.repmeas.fun, param.hat, vval = vv, varfunc = varfunc, 
         q.chisq = q.chisq)
  })
  se.ltdi.1.hat <- apply(ltdi.1.grad, MAR = 2, function(gg) {
    sqrt(drop(t(gg) %*% inv.hess %*% gg))
  })
  ucb.ltdi.1 <- ltdi.1.hat + z.onesided * se.ltdi.1.hat
  
  # method 2
  ltdi.2.hat <- sapply(vgrid, ltdi.2.hetero.repmeas.fun, param = param.hat, 
                       varfunc = varfunc, q.chisq = q.chisq)
  ltdi.2.grad <- sapply(vgrid, function(vv) {
    grad(ltdi.2.hetero.repmeas.fun, param.hat, vval = vv, varfunc = varfunc, 
         q.chisq = q.chisq)
  })
  se.ltdi.2.hat <- apply(ltdi.2.grad, MAR = 2, function(gg) {
    sqrt(drop(t(gg) %*% inv.hess %*% gg))
  })
  ucb.ltdi.2 <- ltdi.2.hat + z.onesided * se.ltdi.2.hat
  
  results.rep.tdi <- cbind(vgrid, ltdi.1.hat, se.ltdi.1.hat, ucb.ltdi.1, 
                           exp(ltdi.1.hat), exp(ucb.ltdi.1), ltdi.2.hat, se.ltdi.2.hat, ucb.ltdi.2, 
                           exp(ltdi.2.hat), exp(ucb.ltdi.2))
  rownames(results.rep.tdi) <- NULL
  colnames(results.rep.tdi) <- c("vgrid", "ltdi.1.hat", "se.ltdi.1.hat", 
                                 "ucb.ltdi.1", "tdi.1.hat", "ucb.tdi.1", "ltdi.2.hat", "se.ltdi.2.hat", 
                                 "ucb.ltdi.2", "tdi.2.hat", "ucb.tdi.2")
  
  # ccc
  # method 1
  zccc.1.hat <- sapply(vgrid, zccc.1.hetero.repmeas.fun, param = param.hat, 
                       varfunc = varfunc)
  zccc.1.grad <- sapply(vgrid, function(vv) {
    grad(zccc.1.hetero.repmeas.fun, param.hat, vval = vv, varfunc = varfunc)
  })
  se.zccc.1.hat <- apply(zccc.1.grad, MAR = 2, function(gg) {
    sqrt(drop(t(gg) %*% inv.hess %*% gg))
  })
  lcb.zccc.1 <- zccc.1.hat - z.onesided * se.zccc.1.hat
  
  # method 2
  zccc.2.hat <- sapply(vgrid, zccc.2.hetero.repmeas.fun, param = param.hat, 
                       varfunc = varfunc)
  zccc.2.grad <- sapply(vgrid, function(vv) {
    grad(zccc.2.hetero.repmeas.fun, param.hat, vval = vv, varfunc = varfunc)
  })
  se.zccc.2.hat <- apply(zccc.2.grad, MAR = 2, function(gg) {
    sqrt(drop(t(gg) %*% inv.hess %*% gg))
  })
  lcb.zccc.2 <- zccc.2.hat - z.onesided * se.zccc.2.hat
  
  results.rep.ccc <- cbind(vgrid, zccc.1.hat, se.zccc.1.hat, lcb.zccc.1, 
                           tanh(zccc.1.hat), tanh(lcb.zccc.1), zccc.2.hat, se.zccc.2.hat, lcb.zccc.2, 
                           tanh(zccc.2.hat), tanh(lcb.zccc.2))
  rownames(results.rep.ccc) <- NULL
  colnames(results.rep.ccc) <- c("vgrid", "zccc.1.hat", "se.zccc.1.hat", 
                                 "lcb.zccc.1", "ccc.1.hat", "lcb.ccc.1", "zccc.2.hat", "se.zccc.2.hat", 
                                 "lcb.zccc.2", "ccc.2.hat", "lcb.ccc.2")
  
  results <- list(param = results.param, lambda = results.lambda, tdi = results.tdi, 
                  ccc = results.ccc, rep.tdi = results.rep.tdi, rep.ccc = results.rep.ccc)
  return(results)
}


######################
#
# Function for computing loglikelihood function as a function of model parameters
# for unlinked repeated measures data
# Arguments: 
#		param = transformed model parameters
#		dat = a data frame with four columns --- subject (a factor), 
#			  method (a factor with two levels method1 and method2), 
#			rep (replication number, unused) and meas
#		vtilde = a vector of length nrow(dat) containing the subject-specific 
#				values of the variance covariate 
#		varfunc = a character with one of two possible values---
#				"power" or "exponential"--- indicating the common variance function
#				model for the two methods
# Value: The loglikelihood function

loglik.hetero.repmeas.data <- function(param, dat, vtilde, varfunc) {
  require(mvtnorm)
  
  # function for computing loglik for one subject
  loglik.unlinked.one.hetero <- function(beta, G, yy, vv, varfunc) {
    # data
    y.1 <- yy[names(yy) == "method1"]
    y.2 <- yy[names(yy) == "method2"]
    m.1 <- length(y.1)
    m.2 <- length(y.2)
    
    # variance covariate
    vv <- unique(vv)
    if (length(vv) != 1) {
      stop("variance covariate in log.hetero.data is not unique to a subject")
    }
    vv.1 <- rep(vv, m.1)
    vv.2 <- rep(vv, m.2)
    
    # design matrices
    X <- cbind(c(rep(1, m.1), rep(0, m.2)), c(rep(0, m.1), rep(1, m.2)))
    Z <- cbind(rep(1, m.1 + m.2), X)
    
    # error variance matrix		
    if (varfunc == "power") {
      err.var.1 <- sigmasq.e1 * (abs(vv.1)^(2 * delta.1))
      err.var.2 <- sigmasq.e2 * (abs(vv.2)^(2 * delta.2))
    }
    if (varfunc == "exponential") {
      err.var.1 <- sigmasq.e1 * exp(2 * delta.1 * vv.1)
      err.var.2 <- sigmasq.e2 * exp(2 * delta.2 * vv.2)
    }
    R <- diag(c(err.var.1, err.var.2))
    
    # loglikelihood
    return(dmvnorm(c(y.1, y.2), mean = (X %*% beta), sigma = (Z %*% G %*% 
                                                                t(Z) + R), log = T))
  }
  
  if (!(varfunc %in% c("power", "exponential"))) {
    stop(paste(varfunc, "is not a valid option in loglik.hetero.repmeas.data"))
  }
  
  # parameters
  alpha <- param["alpha"]
  mu.b <- param["mu.b"]
  sigmasq.b <- exp(param["lsigmasq.b"])
  sigmasq.bI <- exp(param["lsigmasq.bI"])
  sigmasq.e1 <- exp(param["lsigmasq.e1"])
  sigmasq.e2 <- exp(param["lsigmasq.e2"])
  delta.1 <- param["delta.1"]
  delta.2 <- param["delta.2"]
  
  beta <- c(mu.b, alpha + mu.b)
  G <- diag(c(sigmasq.b, sigmasq.bI, sigmasq.bI))
  
  # get data in wide format
  ylist <- long2wide.unlinked(dat)
  vlist <- split(vtilde, dat$subject)
  
  # compute loglik for each subject
  loglik <- mapply(loglik.unlinked.one.hetero, ylist, vlist, MoreArgs = list(beta = beta, 
                                                                             G = G, varfunc = varfunc))
  
  return(sum(as.numeric(loglik)))
}


###########################
#
# Function for fitting a heteroscedastic mixed model to unlinked repeated 
# measurements data
# Arguments: 
#		dat = a data frame with four columns --- subject (a factor), 
#			  method (a factor with two levels method1 and method2), 
#			rep (replication number, unused) and meas
#		vtilde = a vector of length nrow(dat) containing the subject-specific 
#				value of the variance covariate 
#		varfunc = a character with one of two possible values---
#				"power" or "exponential"--- indicating the common variance function
#				model for the two methods
#		vgrid = a grid of values on which to evaluate the fitted error variance 
#					functions and the paramerters that involve them
# Value: The fitted object from the lme fit
#
lme.repmeas.hetero.fit <- function(dat, vtilde, varfunc, vgrid) {
  
  if (!(varfunc %in% c("power", "exponential"))) {
    stop(paste(varfunc, "is not a valid option in lme.unlinked.repmeas.hetero.fit"))
  }
  
  # fit model
  fit <- fit.hetero.lme(dat, vtilde, varfunc)
  
  # extract estimates
  mu.1.hat <- fixef(fit)["methodmethod1"]
  mu.2.hat <- fixef(fit)["methodmethod2"]
  sigma.hat <- fit$sigma
  
  mat1 <- diag(as.matrix(fit$modelStruct$reStruct[[1]]) * sigma.hat^2)
  sigmasq.b.hat <- mat1["(Intercept)"]
  sigmasq.bI.hat <- mat1["methodmethod1"]
  
  sigma.ratio.hat <- exp(as.numeric(fit$modelStruct$varStruct$A[1]))
  
  # determine which level of method is taken as reference in variance function
  wt <- attr(fit$modelStruct$varStruct$A, "weights")
  ref.level <- unique(names(which((wt == 1))))
  
  # get appropriate variance estimates
  if (ref.level == "method1") {
    sigmasq.e1.hat <- sigma.hat^2
    sigmasq.e2.hat <- (sigma.hat * sigma.ratio.hat)^2
  } else {
    sigmasq.e2.hat <- sigma.hat^2
    sigmasq.e1.hat <- (sigma.hat * sigma.ratio.hat)^2
  }
  
  # get the heteroscedasticity parameters
  delta.hat <- as.numeric(as.numeric(fit$modelStruct$varStruct$B))
  names(delta.hat) <- attr(fit$modelStruct$varStruct$A, "groupNames")
  delta.1.hat <- delta.hat["method1"]
  delta.2.hat <- delta.hat["method2"]
  
  # cast estimates in desired form
  mu.b.hat <- mu.1.hat
  alpha.hat <- mu.2.hat - mu.1.hat
  
  param.hat <- c(alpha.hat, mu.b.hat, log(c(sigmasq.b.hat, sigmasq.bI.hat, 
                                            sigmasq.e1.hat, sigmasq.e2.hat)), delta.1.hat, delta.2.hat)
  names(param.hat) <- c("alpha", "mu.b", "lsigmasq.b", "lsigmasq.bI", "lsigmasq.e1", 
                        "lsigmasq.e2", "delta.1", "delta.2")
  
  # estimated variance parameters
  var.hat <- c(sigmasq.b.hat, sigmasq.bI.hat, sigmasq.e1.hat, sigmasq.e2.hat)
  names(var.hat) <- c("sigmasq.b", "sigmasq.bI", "sigmasq.e1", "sigmasq.e2")
  
  # estimated error variances
  if (varfunc == "power") {
    err.var.1.hat <- sigmasq.e1.hat * (abs(vgrid)^(2 * delta.1.hat))
    err.var.2.hat <- sigmasq.e2.hat * (abs(vgrid)^(2 * delta.2.hat))
  }
  
  if (varfunc == "exponential") {
    err.var.1.hat <- sigmasq.e1.hat * exp(2 * delta.1.hat * vgrid)
    err.var.2.hat <- sigmasq.e2.hat * exp(2 * delta.2.hat * vgrid)
  }
  
  err.var.hat <- cbind(err.var.1.hat, err.var.2.hat)
  colnames(err.var.hat) <- c("method1", "method2")
  
  # fitted bivariate distribution 
  var.1.hat <- sigmasq.b.hat + sigmasq.bI.hat + err.var.1.hat
  var.2.hat <- sigmasq.b.hat + sigmasq.bI.hat + err.var.2.hat
  
  corr.hat <- sigmasq.b.hat/sqrt(var.1.hat * var.2.hat)
  bvn.hat <- list(mean1 = mu.1.hat, mean2 = mu.2.hat, var1 = var.1.hat, var2 = var.2.hat, 
                  corr = corr.hat)
  
  # fitted distribution to the difference
  var.D.hat <- 2 * sigmasq.bI.hat + err.var.1.hat + err.var.2.hat
  
  # reliabilities
  rlbty.1.hat <- (sigmasq.b.hat + sigmasq.bI.hat)/var.1.hat
  rlbty.2.hat <- (sigmasq.b.hat + sigmasq.bI.hat)/var.2.hat
  rlbty.hat <- cbind(rlbty.1.hat, rlbty.2.hat)
  colnames(rlbty.hat) <- c("method1", "method2")
  
  # limits of agreement
  loa.inter <- matrix(0, nrow = length(vgrid), ncol = 2)
  colnames(loa.inter) <- c("lower", "upper")
  loa.inter[, 1] <- alpha.hat - 1.96 * sqrt(var.D.hat)
  loa.inter[, 2] <- alpha.hat + 1.96 * sqrt(var.D.hat)
  
  loa.intra.1 <- matrix(0, nrow = length(vgrid), ncol = 2)
  colnames(loa.intra.1) <- c("lower", "upper")
  loa.intra.1[, 1] <- -1.96 * sqrt(2 * err.var.1.hat)
  loa.intra.1[, 2] <- 1.96 * sqrt(2 * err.var.1.hat)
  
  loa.intra.2 <- matrix(0, nrow = length(vgrid), ncol = 2)
  colnames(loa.intra.2) <- c("lower", "upper")
  loa.intra.2[, 1] <- -1.96 * sqrt(2 * err.var.2.hat)
  loa.intra.2[, 2] <- 1.96 * sqrt(2 * err.var.2.hat)
  
  # compare with the homoscedastic model
  fit0 <- lme.unlinked.repmeas.fit(dat)
  anova.results <- anova(fit0$fit, fit)
  
  return(list(varfunc = varfunc, fit = fit, param.hat = param.hat, anova.results = anova.results, 
              vgrid = vgrid, err.var.hat = err.var.hat, var.hat = var.hat, bvn.hat = bvn.hat, 
              var.D.hat = var.D.hat, reliability.hat = rlbty.hat, loa = list(loa.inter = loa.inter, loa.intra.1 = loa.intra.1, 
                                                                             loa.intra.2 = loa.intra.2)))
}


###########################
#
# Function for fitting a heteroscedastic mixed model to unlinked repeated 
# measurements data
# Arguments: 
#		dat = a data frame with four columns --- subject (a factor), 
#			  method (a factor with two levels method1 and method2), 
#			rep (replication number, unused) and meas
#		vtilde = a vector of length nrow(dat) containing the subject-specific 
#				value of the variance covariate 
#		varfunc = a character with one of two possible values---
#				"power" or "exponential"--- indicating the common variance function
#				model for the two methods
# Value: The fitted object from the lme fit
#


fit.hetero.lme <- function(dat, vtilde, varfunc) {
  require(nlme)
  
  if (!(varfunc %in% c("power", "exponential"))) {
    stop(paste(varfunc, "is not a valid option in fit.hetero.lme"))
  }
  
  # fit model
  dat <- data.frame(dat, vtilde)
  dat <- groupedData(meas ~ method | subject, data = dat, order.groups = FALSE)
  
  if (varfunc == "power") {
    fit.power <- lme(meas ~ method - 1, data = dat, random = pdBlocked(list(~1, 
                                                                            pdIdent(~method - 1))), method = "ML", weights = varComb(varIdent(form = ~1 | 
                                                                                                                                                method), varPower(form = ~vtilde | method)), control = lmeControl(maxIter = 1000, 
                                                                                                                                                                                                                  msMaxIter = 5000, niterEM = 2000, msMaxEval = 1500))
    return(fit.power)
  }
  
  if (varfunc == "exponential") {
    fit.exp <- lme(meas ~ method - 1, data = dat, random = pdBlocked(list(~1, 
                                                                          pdIdent(~method - 1))), method = "ML", weights = varComb(varIdent(form = ~1 | 
                                                                                                                                              method), varExp(form = ~vtilde | method)), control = lmeControl(maxIter = 1000, 
                                                                                                                                                                                                              msMaxIter = 5000, niterEM = 2000, msMaxEval = 1500))
    return(fit.exp)
  }
  
}
###########################
#
# Function for obtaining variance covariate for heteroscedastic fit 
# Arguments: 
#		dat = a data frame with four columns --- subject (a factor), 
#			  method (a factor with levels method1 and method2), 
#			  rep (replication number, not used), and meas
#		option= a character with one of two possible values--- 
#				"method1.mean" or "avg.mean"
#
get.var.covariate <- function(dat, option) {
  
  g <- with(dat, interaction(subject, method))
  mean.mat <- with(dat, by(meas, list(subject, method), mean))
  
  if (!(option %in% c("method1.mean", "avg.mean"))) {
    stop(paste(option, "is not a valid option in get.var.covariate"))
  }
  if (option == "method1.mean") {
    v <- mean.mat[, "method1"]
  }
  if (option == "avg.mean") {
    v <- rowMeans(mean.mat)
  }
  return(unsplit(v, g))
}

###########################

####################################################################################
# chapter 5
####################################################################################

#
# Function for computing estimates, their SEs and confidence intervals
#		 for various parameters and functions
# Arguments: 
#		dat = a data frame with four columns --- subject (a factor), 
#			  method (a factor with two levels method1 and method2), 
#			rep (replication number, a factor) and meas
#		param.hat = vector of ML estimates of transformed parameters
# Value: 
#		standard errors of estimates of parameters and their confidence intervals
#		estimate, se and ucb of log(TDI), and estimate and ucb of TDI
#		estimate, se and lcb of z(CCC), and estimate and lcb of CCC
#		estimates, se and ucb of log(TDI.1) and log(TDI.2), and estimates and ucb for 
#				TDI.1 and TDI.2
#		estimates, se and lcb for z(CCC.1) and z(CCC.2), and estimates and lcb for 
#				CCC.1 and CCC.2
#

######################
conf.measures.linked <- function(dat, param.hat, conf.level = 0.95, prob = 0.9) {
  require(numDeriv)
  
  z.twosided <- qnorm((1 + conf.level)/2)
  z.onesided <- qnorm(conf.level)
  q.chisq <- qchisq(prob, 1)
  
  # get Hessian and invert it
  hess <- -hessian(loglik.linked.data, param.hat, dat = dat)
  inv.hess <- solve(hess)
  rownames(inv.hess) <- names(param.hat)
  colnames(inv.hess) <- names(param.hat)
  
  # SEs of param.hat and confidence intevals
  se.param.hat <- sqrt(diag(inv.hess))
  ci.param <- cbind(param.hat - z.twosided * se.param.hat, param.hat + z.twosided * 
                      se.param.hat)
  results1 <- cbind(param.hat, se.param.hat, ci.param)
  rownames(results1) <- names(param.hat)
  colnames(results1) <- c("estimate", "se", "lcl", "ucl")
  
  # estimate and confidence interval for precision ratio
  llambda.hat <- param.hat["lsigmasq.e1"] - param.hat["lsigmasq.e2"]
  se.llambda.hat <- sqrt(inv.hess["lsigmasq.e1", "lsigmasq.e1"] + inv.hess["lsigmasq.e2", 
                                                                           "lsigmasq.e2"] - 2 * inv.hess["lsigmasq.e1", "lsigmasq.e2"])
  ci.llambda <- llambda.hat + c(-1, 1) * z.twosided * se.llambda.hat
  
  # estimate and confidence interval for new precision ratio
  lnewprecratio.hat <- lnewprec.ratio.linked.fun(param.hat)
  lnewprecratio.grad <- grad(lnewprec.ratio.linked.fun, param.hat)
  se.lnewprecratio.hat <- sqrt(drop(t(lnewprecratio.grad) %*% inv.hess %*% 
                                      lnewprecratio.grad))
  ci.lnewprecratio <- lnewprecratio.hat + c(-1, 1) * z.twosided * se.lnewprecratio.hat
  
  # combine the last two results
  
  results2 <- rbind(c(llambda.hat, se.llambda.hat, ci.llambda, exp(llambda.hat), 
                      exp(ci.llambda)), c(lnewprecratio.hat, se.lnewprecratio.hat, ci.lnewprecratio, 
                                          exp(lnewprecratio.hat), exp(ci.lnewprecratio)))
  rownames(results2) <- c("lambda", "newprecratio")
  colnames(results2) <- c("lest", "lest.se", "lcl.lpar", "ucl.lpar", "est", 
                          "lcl.par", "ucl.par")
  
  # estimate and upper confidence bound of TDI
  ltdi.hat <- ltdi.linked.fun(param.hat, prob)
  ltdi.grad <- grad(ltdi.linked.fun, param.hat, prob = prob)
  se.ltdi.hat <- sqrt(drop(t(ltdi.grad) %*% inv.hess %*% ltdi.grad))
  ucb.ltdi <- ltdi.hat + z.onesided * se.ltdi.hat
  results3 <- c(ltdi.hat, se.ltdi.hat, ucb.ltdi, exp(ltdi.hat), exp(ucb.ltdi))
  names(results3) <- c("ltdi.hat", "se.ltdi.hat", "ucb.ltdi", "tdi.hat", "ucb.tdi")
  
  # estimate and lower confidence bound of CCC
  zccc.hat <- zccc.linked.fun(param.hat)
  zccc.grad <- grad(zccc.linked.fun, param.hat)
  se.zccc.hat <- sqrt(drop(t(zccc.grad) %*% inv.hess %*% zccc.grad))
  lcb.zccc <- zccc.hat - z.onesided * se.zccc.hat
  results4 <- c(zccc.hat, se.zccc.hat, lcb.zccc, tanh(zccc.hat), tanh(lcb.zccc))
  names(results4) <- c("zccc.hat", "se.zccc.hat", "lcb.zccc", "ccc.hat", "lcb.ccc")
  
  
  # evaluation of repetability
  # tdi
  ltdi.1.hat <- ltdi.1.linked.fun(param.hat, prob)
  ltdi.1.grad <- grad(ltdi.1.linked.fun, param.hat, prob = prob)
  se.ltdi.1.hat <- sqrt(drop(t(ltdi.1.grad) %*% inv.hess %*% ltdi.1.grad))
  ucb.ltdi.1 <- ltdi.1.hat + z.onesided * se.ltdi.1.hat
  ltdi.2.hat <- ltdi.2.linked.fun(param.hat, prob)
  ltdi.2.grad <- grad(ltdi.2.linked.fun, param.hat, prob = prob)
  se.ltdi.2.hat <- sqrt(drop(t(ltdi.2.grad) %*% inv.hess %*% ltdi.2.grad))
  ucb.ltdi.2 <- ltdi.2.hat + z.onesided * se.ltdi.2.hat
  results5 <- rbind(c(ltdi.1.hat, se.ltdi.1.hat, ucb.ltdi.1, exp(ltdi.1.hat), 
                      exp(ucb.ltdi.1)), c(ltdi.2.hat, se.ltdi.2.hat, ucb.ltdi.2, exp(ltdi.2.hat), 
                                          exp(ucb.ltdi.2)))
  rownames(results5) <- c("method1", "method2")
  colnames(results5) <- c("ltdi.hat", "se.ltdi.hat", "ucb.ltdi", "tdi.hat", 
                          "ucb.tdi")
  
  # ccc
  zccc.1.hat <- zccc.1.linked.fun(param.hat)
  zccc.1.grad <- grad(zccc.1.linked.fun, param.hat)
  se.zccc.1.hat <- sqrt(drop(t(zccc.1.grad) %*% inv.hess %*% zccc.1.grad))
  lcb.zccc.1 <- zccc.1.hat - z.onesided * se.zccc.1.hat
  
  zccc.2.hat <- zccc.2.linked.fun(param.hat)
  zccc.2.grad <- grad(zccc.2.linked.fun, param.hat)
  se.zccc.2.hat <- sqrt(drop(t(zccc.2.grad) %*% inv.hess %*% zccc.2.grad))
  lcb.zccc.2 <- zccc.2.hat - z.onesided * se.zccc.2.hat
  
  results6 <- rbind(c(zccc.1.hat, se.zccc.1.hat, lcb.zccc.1, tanh(zccc.1.hat), 
                      tanh(lcb.zccc.1)), c(zccc.2.hat, se.zccc.2.hat, lcb.zccc.2, tanh(zccc.2.hat), 
                                           tanh(lcb.zccc.2)))
  rownames(results6) <- c("method1", "method2")
  colnames(results6) <- c("zccc.hat", "se.zccc.hat", "lcb.zccc", "ccc.hat", 
                          "lcb.ccc")
  
  results <- list(param = results1, similarity = results2, tdi = results3, 
                  ccc = results4, rep.tdi = results5, rep.ccc = results6)
  return(results)
}


######################
#
# Function for computing log-variance ratio as a function of mixed-effects model parameters
# Arguments: 
#		lme.param = transformed model parameters
#
lnewprec.ratio.linked.fun <- function(param) {
  numer <- exp(param["lsigmasq.bI"]) + exp(param["lsigmasq.e1"])
  denom <- exp(param["lsigmasq.bI"]) + exp(param["lsigmasq.e2"])
  return(log(numer/denom))
}

######################
#
# Function for computing log(TDI.1) as a function of model parameters
# Arguments: 
#		param = transformed model parameters
# 		prob = probability
#
ltdi.1.linked.fun <- function(param, prob) {
  # extract the parameters
  sigmasq.bstar <- exp(param["lsigmasq.bstar"])
  sigmasq.e1 <- exp(param["lsigmasq.e1"])
  result <- ltdi(0, sqrt(2 * (sigmasq.bstar + sigmasq.e1)), prob)
  return(result)
}

######################
#
# Function for computing log(TDI.1) as a function of model parameters
# Arguments: 
#		param = transformed model parameters
# 		prob = probability
#
ltdi.2.linked.fun <- function(param, prob) {
  # extract the parameters
  sigmasq.bstar <- exp(param["lsigmasq.bstar"])
  sigmasq.e2 <- exp(param["lsigmasq.e2"])
  result <- ltdi(0, sqrt(2 * (sigmasq.bstar + sigmasq.e2)), prob)
  return(result)
}


######################
#
# Function for computing log(TDI) as a function of model parameters
# Arguments: 
#		param = transformed model parameters
# 		prob = probability
#
ltdi.linked.fun <- function(param, prob) {
  # extract the parameters
  alpha <- param["alpha"]
  sigmasq.bI <- exp(param["lsigmasq.bI"])
  sigmasq.e1 <- exp(param["lsigmasq.e1"])
  sigmasq.e2 <- exp(param["lsigmasq.e2"])
  result <- ltdi(alpha, sqrt(2 * sigmasq.bI + sigmasq.e1 + sigmasq.e2), prob)
  return(result)
}
##################
#
# Function for computing Fisher's z-transformation of CCC.2 as a function of model parameters
# Arguments: 
#		param = transformed model parameters
#
zccc.2.linked.fun <- function(param) {
  # extract the parameters
  sigmasq.b <- exp(param["lsigmasq.b"])
  sigmasq.bI <- exp(param["lsigmasq.bI"])
  sigmasq.bstar <- exp(param["lsigmasq.bstar"])
  sigmasq.e2 <- exp(param["lsigmasq.e2"])
  cov.mat <- matrix(c(sigmasq.b + sigmasq.bI + sigmasq.bstar + sigmasq.e2, 
                      sigmasq.b + sigmasq.bI, sigmasq.b + sigmasq.bI, sigmasq.b + sigmasq.bI + 
                        sigmasq.bstar + sigmasq.e2), byrow = T, ncol = 2)
  result <- zccc(0, cov.mat[1, 1], cov.mat[2, 2], cov.mat[1, 2])
  return(result)
}

##################
#
# Function for computing Fisher's z-transformation of CCC.1 as a function of model parameters
# Arguments: 
#		param = transformed model parameters
#
zccc.1.linked.fun <- function(param) {
  # extract the parameters
  sigmasq.b <- exp(param["lsigmasq.b"])
  sigmasq.bI <- exp(param["lsigmasq.bI"])
  sigmasq.bstar <- exp(param["lsigmasq.bstar"])
  sigmasq.e1 <- exp(param["lsigmasq.e1"])
  cov.mat <- matrix(c(sigmasq.b + sigmasq.bI + sigmasq.bstar + sigmasq.e1, 
                      sigmasq.b + sigmasq.bI, sigmasq.b + sigmasq.bI, sigmasq.b + sigmasq.bI + 
                        sigmasq.bstar + sigmasq.e1), byrow = T, ncol = 2)
  result <- zccc(0, cov.mat[1, 1], cov.mat[2, 2], cov.mat[1, 2])
  return(result)
}

##################
#
# Function for computing Fisher's z-transformation of CCC as a function of model parameters
# Arguments: 
#		param = transformed model parameters
#
zccc.linked.fun <- function(param) {
  # extract the parameters
  alpha <- param["alpha"]
  sigmasq.b <- exp(param["lsigmasq.b"])
  sigmasq.bI <- exp(param["lsigmasq.bI"])
  sigmasq.bstar <- exp(param["lsigmasq.bstar"])
  sigmasq.e1 <- exp(param["lsigmasq.e1"])
  sigmasq.e2 <- exp(param["lsigmasq.e2"])
  cov.mat <- matrix(c(sigmasq.b + sigmasq.bI + sigmasq.bstar + sigmasq.e1, 
                      sigmasq.b + sigmasq.bstar, sigmasq.b + sigmasq.bstar, sigmasq.b + sigmasq.bI + 
                        sigmasq.bstar + sigmasq.e2), byrow = T, ncol = 2)
  result <- zccc(alpha, cov.mat[1, 1], cov.mat[2, 2], cov.mat[1, 2])
  return(result)
}

##################


######################
#
# Function for computing loglikelihood function as a function of model parameters
# for linked data
# Arguments: 
#		param = transformed model parameters
# 		dat = data in long format, used for model fitting
#
loglik.linked.data <- function(param, dat) {
  require(mvtnorm)
  # function for computing loglik for one subject
  loglik.linked.one <- function(yy) {
    # data
    yy1 <- subset(yy, method == "method1")
    yy2 <- subset(yy, method == "method2")
    yy1 <- yy1[order(yy1$rep), ]
    yy2 <- yy2[order(yy2$rep), ]
    ind1 <- yy1$rep
    ind2 <- yy2$rep
    if (!identical(ind1, ind2)) {
      stop("the methods don't have identical repeat numbers on a subject")
    }
    y <- c(yy1$meas, yy2$meas)
    m <- nrow(yy1)
    # design matrices
    X <- cbind(c(rep(1, m), rep(0, m)), c(rep(0, m), rep(1, m)))
    Z <- rbind(cbind(rep(1, m), rep(1, m), rep(0, m), diag(m)), cbind(rep(1, 
                                                                          m), rep(0, m), rep(1, m), diag(m)))
    # error variance matrix
    R <- diag(c(rep(sigmasq.e1, m), rep(sigmasq.e2, m)))
    # random effect covariance matrix
    G <- diag(c(sigmasq.b, sigmasq.bI, sigmasq.bI, rep(sigmasq.bstar, m)))
    # loglikelihood
    return(dmvnorm(y, mean = (X %*% beta), sigma = (Z %*% G %*% t(Z) + R), 
                   log = T))
  }
  
  # parameters
  alpha <- param["alpha"]
  mu.b <- param["mu.b"]
  sigmasq.b <- exp(param["lsigmasq.b"])
  sigmasq.bI <- exp(param["lsigmasq.bI"])
  sigmasq.bstar <- exp(param["lsigmasq.bstar"])
  sigmasq.e1 <- exp(param["lsigmasq.e1"])
  sigmasq.e2 <- exp(param["lsigmasq.e2"])
  beta <- c(mu.b, alpha + mu.b)
  
  # split the data in a list by subject id
  ylist <- split(dat, dat$subject)
  
  # compute loglik for each subject
  loglik <- lapply(ylist, loglik.linked.one)
  
  return(sum(as.numeric(loglik)))
}
##################
#
# Function for fitting a mixed model to linked repeated measurements data
# Arguments: 
#		dat = a data frame with four columns --- subject (a factor), 
#			  method (a factor with two levels method1 and method2), 
#			rep (replication number, a factor) and meas
# Value: A list with following components
#		fit = fitted lme object
#		param.hat = ML estimates of transformed parameters
#		var.hat = ML estimates of variance components and variance of difference
#		bvn.hat = parameters of fitted bivariate normal distribution
#		reliability.hat = ML estimates of reliabilities of the two methods
#

lme.linked.repmeas.fit <- function(dat) {
  require(nlme)
  # fit model
  dat <- groupedData(meas ~ method | subject, data = dat, order.groups = FALSE)
  fit <- lme(meas ~ method - 1, data = dat, weights = varIdent(form = ~1 | 
                                                                 method), method = "ML", random = pdBlocked(list(~1, pdIdent(~method - 
                                                                                                                               1), pdIdent(~rep - 1))))
  
  
  # extract estimates
  mu.1.hat <- fixef(fit)["methodmethod1"]
  mu.2.hat <- fixef(fit)["methodmethod2"]
  sigma.hat <- fit$sigma
  
  mat1 <- diag(as.matrix(fit$modelStruct$reStruct[[1]]) * sigma.hat^2)
  sigmasq.b.hat <- mat1["(Intercept)"]
  sigmasq.bI.hat <- mat1["methodmethod1"]
  sigmasq.bstar.hat <- mat1["rep1"]
  
  sigma.ratio.hat <- exp(as.numeric(fit$modelStruct$varStruct[1]))
  
  # determine which level of method is taken as reference in variance function
  wt <- attr(fit$modelStruct$varStruct, "weights")
  ref.level <- unique(names(which((wt == 1))))
  
  # get appropriate variance estimates
  if (ref.level == "method1") {
    sigmasq.e1.hat <- sigma.hat^2
    sigmasq.e2.hat <- (sigma.hat * sigma.ratio.hat)^2
  } else {
    sigmasq.e2.hat <- sigma.hat^2
    sigmasq.e1.hat <- (sigma.hat * sigma.ratio.hat)^2
  }
  # cast estimates in desired form
  mu.b.hat <- mu.1.hat
  alpha.hat <- mu.2.hat - mu.1.hat
  
  param.hat <- c(alpha.hat, mu.b.hat, log(c(sigmasq.b.hat, sigmasq.bI.hat, 
                                            sigmasq.bstar.hat, sigmasq.e1.hat, sigmasq.e2.hat)))
  names(param.hat) <- c("alpha", "mu.b", "lsigmasq.b", "lsigmasq.bI", "lsigmasq.bstar", 
                        "lsigmasq.e1", "lsigmasq.e2")
  
  # estimated variances
  var.hat <- c(sigmasq.b.hat, sigmasq.bI.hat, sigmasq.bstar.hat, sigmasq.e1.hat, 
               sigmasq.e2.hat, 2 * sigmasq.bI.hat + sigmasq.e1.hat + sigmasq.e2.hat)
  names(var.hat) <- c("sigmasq.b", "sigmasq.bI", "sigmasq.bstar", "sigmasq.e1", 
                      "sigmasq.e2", "sigmasq")
  
  # fitted bivariate distribution 
  var.1.hat <- sigmasq.e1.hat + sigmasq.bstar.hat + sigmasq.bI.hat + sigmasq.b.hat
  var.2.hat <- sigmasq.e2.hat + sigmasq.bstar.hat + sigmasq.bI.hat + sigmasq.b.hat
  corr.hat <- (sigmasq.b.hat + sigmasq.bstar.hat)/sqrt(var.1.hat * var.2.hat)
  bvn.hat <- c(mu.1.hat, mu.2.hat, var.1.hat, var.2.hat, corr.hat)
  names(bvn.hat) <- c("mean1", "mean2", "var1", "var2", "corr")
  
  # reliabilities
  
  rlbty.1.hat <- (sigmasq.b.hat + sigmasq.bstar.hat + sigmasq.bI.hat)/var.1.hat
  rlbty.2.hat <- (sigmasq.b.hat + sigmasq.bstar.hat + sigmasq.bI.hat)/var.2.hat
  rlbty.hat <- c(rlbty.1.hat, rlbty.2.hat)
  names(rlbty.hat) <- c("method1", "method2")
  
  return(list(fit = fit, param.hat = param.hat, var.hat = var.hat, bvn.hat = bvn.hat, 
              reliability.hat = rlbty.hat))
}

##################
#
# Function for computing log-variance ratio as a function of mixed-effects model parameters
# Arguments: 
#		lme.param = transformed model parameters
#
lnewprec.ratio.unlinked.fun <- function(param) {
  numer <- exp(param["lsigmasq.bI"]) + exp(param["lsigmasq.e1"])
  denom <- exp(param["lsigmasq.bI"]) + exp(param["lsigmasq.e2"])
  return(log(numer/denom))
}


##################
#
# Function for computing estimates, their SEs and confidence intervals
#		 for various parameters and functions
# Arguments: 
#		dat = a data frame with four columns --- subject (a factor), 
#			  method (a factor with two levels method1 and method2), 
#			rep (replication number, unused) and meas
#		param.hat = vector of ML estimates of transformed parameters
# Value: 
#		standard errors of estimates of parameters and their confidence intervals
#		estimate, se and ucb of log(TDI), and estimate and ucb of TDI
#		estimate, se and lcb of z(CCC), and estimate and lcb of CCC
#		estimates, se and ucb of log(TDI.1) and log(TDI.2), and estimates and ucb for 
#				TDI.1 and TDI.2
#		estimates, se and lcb for z(CCC.1) and z(CCC.2), and estimates and lcb for 
#				CCC.1 and CCC.2
#
##################################

conf.measures.unlinked <- function(dat, param.hat, conf.level = 0.95, prob = 0.9) {
  require(numDeriv)
  
  z.twosided <- qnorm((1 + conf.level)/2)
  z.onesided <- qnorm(conf.level)
  q.chisq <- qchisq(prob, 1)
  
  # get Hessian and invert it
  hess <- -hessian(loglik.unlinked.data, param.hat, dat = dat)
  inv.hess <- solve(hess)
  rownames(inv.hess) <- names(param.hat)
  colnames(inv.hess) <- names(param.hat)
  
  # SEs of param.hat and confidence intevals
  se.param.hat <- sqrt(diag(inv.hess))
  ci.param <- cbind(param.hat - z.twosided * se.param.hat, param.hat + z.twosided * 
                      se.param.hat)
  results1 <- cbind(param.hat, se.param.hat, ci.param)
  rownames(results1) <- names(param.hat)
  colnames(results1) <- c("estimate", "se", "lcl", "ucl")
  
  # estimate and confidence interval for precision ratio
  llambda.hat <- param.hat["lsigmasq.e1"] - param.hat["lsigmasq.e2"]
  se.llambda.hat <- sqrt(inv.hess["lsigmasq.e1", "lsigmasq.e1"] + inv.hess["lsigmasq.e2", 
                                                                           "lsigmasq.e2"] - 2 * inv.hess["lsigmasq.e1", "lsigmasq.e2"])
  ci.llambda <- llambda.hat + c(-1, 1) * z.twosided * se.llambda.hat
  
  # estimate and confidence interval for new precision ratio
  lnewprecratio.hat <- lnewprec.ratio.unlinked.fun(param.hat)
  lnewprecratio.grad <- grad(lnewprec.ratio.unlinked.fun, param.hat)
  se.lnewprecratio.hat <- sqrt(drop(t(lnewprecratio.grad) %*% inv.hess %*% 
                                      lnewprecratio.grad))
  ci.lnewprecratio <- lnewprecratio.hat + c(-1, 1) * z.twosided * se.lnewprecratio.hat
  
  # combine the last two results
  
  results2 <- rbind(c(llambda.hat, se.llambda.hat, ci.llambda, exp(llambda.hat), 
                      exp(ci.llambda)), c(lnewprecratio.hat, se.lnewprecratio.hat, ci.lnewprecratio, 
                                          exp(lnewprecratio.hat), exp(ci.lnewprecratio)))
  rownames(results2) <- c("lambda", "newprecratio")
  colnames(results2) <- c("lest", "lest.se", "lcl.lpar", "ucl.lpar", "est", 
                          "lcl.par", "ucl.par")
  
  # estimate and upper confidence bound of TDI
  ltdi.hat <- ltdi.unlinked.fun(param.hat, prob)
  ltdi.grad <- grad(ltdi.unlinked.fun, param.hat, prob = prob)
  se.ltdi.hat <- sqrt(drop(t(ltdi.grad) %*% inv.hess %*% ltdi.grad))
  ucb.ltdi <- ltdi.hat + z.onesided * se.ltdi.hat
  results3 <- c(ltdi.hat, se.ltdi.hat, ucb.ltdi, exp(ltdi.hat), exp(ucb.ltdi))
  names(results3) <- c("ltdi.hat", "se.ltdi.hat", "ucb.ltdi", "tdi.hat", "ucb.tdi")
  
  # estimate and lower confidence bound of CCC
  zccc.hat <- zccc.unlinked.fun(param.hat)
  zccc.grad <- grad(zccc.unlinked.fun, param.hat)
  se.zccc.hat <- sqrt(drop(t(zccc.grad) %*% inv.hess %*% zccc.grad))
  lcb.zccc <- zccc.hat - z.onesided * se.zccc.hat
  results4 <- c(zccc.hat, se.zccc.hat, lcb.zccc, tanh(zccc.hat), tanh(lcb.zccc))
  names(results4) <- c("zccc.hat", "se.zccc.hat", "lcb.zccc", "ccc.hat", "lcb.ccc")
  
  
  # evaluation of repetability
  # tdi
  ltdi.1.hat <- (1/2) * (log(2 * q.chisq) + param.hat["lsigmasq.e1"])
  se.ltdi.1.hat <- (1/2) * se.param.hat["lsigmasq.e1"]
  ucb.ltdi.1 <- ltdi.1.hat + z.onesided * se.ltdi.1.hat
  ltdi.2.hat <- (1/2) * (log(2 * q.chisq) + param.hat["lsigmasq.e2"])
  se.ltdi.2.hat <- (1/2) * se.param.hat["lsigmasq.e2"]
  ucb.ltdi.2 <- ltdi.2.hat + z.onesided * se.ltdi.2.hat
  results5 <- rbind(c(ltdi.1.hat, se.ltdi.1.hat, ucb.ltdi.1, exp(ltdi.1.hat), 
                      exp(ucb.ltdi.1)), c(ltdi.2.hat, se.ltdi.2.hat, ucb.ltdi.2, exp(ltdi.2.hat), 
                                          exp(ucb.ltdi.2)))
  rownames(results5) <- c("method1", "method2")
  colnames(results5) <- c("ltdi.hat", "se.ltdi.hat", "ucb.ltdi", "tdi.hat", 
                          "ucb.tdi")
  
  # ccc
  zccc.1.hat <- zccc.1.unlinked.fun(param.hat)
  zccc.1.grad <- grad(zccc.1.unlinked.fun, param.hat)
  se.zccc.1.hat <- sqrt(drop(t(zccc.1.grad) %*% inv.hess %*% zccc.1.grad))
  lcb.zccc.1 <- zccc.1.hat - z.onesided * se.zccc.1.hat
  
  zccc.2.hat <- zccc.2.unlinked.fun(param.hat)
  zccc.2.grad <- grad(zccc.2.unlinked.fun, param.hat)
  se.zccc.2.hat <- sqrt(drop(t(zccc.2.grad) %*% inv.hess %*% zccc.2.grad))
  lcb.zccc.2 <- zccc.2.hat - z.onesided * se.zccc.2.hat
  
  results6 <- rbind(c(zccc.1.hat, se.zccc.1.hat, lcb.zccc.1, tanh(zccc.1.hat), 
                      tanh(lcb.zccc.1)), c(zccc.2.hat, se.zccc.2.hat, lcb.zccc.2, tanh(zccc.2.hat), 
                                           tanh(lcb.zccc.2)))
  rownames(results6) <- c("method1", "method2")
  colnames(results6) <- c("zccc.hat", "se.zccc.hat", "lcb.zccc", "ccc.hat", 
                          "lcb.ccc")
  
  results <- list(param = results1, similarity = results2, tdi = results3, 
                  ccc = results4, rep.tdi = results5, rep.ccc = results6)
  return(results)
}

######################
#
# Function for computing loglikelihood function as a function of model parameters
# for unlinked data
# Arguments: 
#		param = transformed model parameters
# 		dat = data in long format, used for model fitting
#

loglik.unlinked.data <- function(param, dat) {
  require(mvtnorm)
  # function for computing loglik for one subject
  loglik.unlinked.one <- function(beta, G, yy) {
    # data
    y.1 <- yy[names(yy) == "method1"]
    y.2 <- yy[names(yy) == "method2"]
    m.1 <- length(y.1)
    m.2 <- length(y.2)
    # design matrices
    X <- cbind(c(rep(1, m.1), rep(0, m.2)), c(rep(0, m.1), rep(1, m.2)))
    Z <- cbind(rep(1, m.1 + m.2), X)
    # error variance matrix
    R <- diag(c(rep(sigmasq.e1, m.1), rep(sigmasq.e2, m.2)))
    # loglikelihood
    return(dmvnorm(c(y.1, y.2), mean = (X %*% beta), sigma = (Z %*% G %*% 
                                                                t(Z) + R), log = T))
  }
  
  # parameters
  alpha <- param["alpha"]
  mu.b <- param["mu.b"]
  sigmasq.b <- exp(param["lsigmasq.b"])
  sigmasq.bI <- exp(param["lsigmasq.bI"])
  sigmasq.e1 <- exp(param["lsigmasq.e1"])
  sigmasq.e2 <- exp(param["lsigmasq.e2"])
  beta <- c(mu.b, alpha + mu.b)
  G <- diag(c(sigmasq.b, sigmasq.bI, sigmasq.bI))
  
  # get data in wide format
  ylist <- long2wide.unlinked(dat)
  
  # compute loglik for each subject
  loglik <- lapply(ylist, loglik.unlinked.one, beta = beta, G = G)
  
  return(sum(as.numeric(loglik)))
}



######################
#
# Function for coverting data in long format, used for model fitting, to wide
# format useful for computing log-likelihood
# Arguments: 
#		dat = data in long format
# Value:
#	data in wide format, one list for each subject
#
long2wide.unlinked <- function(dat) {
  assign.method <- function(obs, method.ind) {
    names(obs) <- method.ind
    return(obs)
  }
  # split the meas and method columns	
  obs <- split(dat$meas, dat$subject)
  method.ind <- split(dat$method, dat$subject)
  return(mapply(assign.method, obs, method.ind, SIMPLIFY = FALSE))
}

##################
#
# Function for computing Fisher's z-transformation of CCC.2 as a function of model parameters
# Arguments: 
#		param = transformed model parameters
#
zccc.2.unlinked.fun <- function(param) {
  # extract the parameters
  sigmasq.b <- exp(param["lsigmasq.b"])
  sigmasq.bI <- exp(param["lsigmasq.bI"])
  sigmasq.e2 <- exp(param["lsigmasq.e2"])
  cov.mat <- matrix(c(sigmasq.b + sigmasq.bI + sigmasq.e2, sigmasq.b + sigmasq.bI, 
                      sigmasq.b + sigmasq.bI, sigmasq.b + sigmasq.bI + sigmasq.e2), byrow = T, 
                    ncol = 2)
  result <- zccc(0, cov.mat[1, 1], cov.mat[2, 2], cov.mat[1, 2])
  return(result)
}

##################
#
# Function for computing Fisher's z-transformation of CCC.1 as a function of model parameters
# Arguments: 
#		param = transformed model parameters
#
zccc.1.unlinked.fun <- function(param) {
  # extract the parameters
  sigmasq.b <- exp(param["lsigmasq.b"])
  sigmasq.bI <- exp(param["lsigmasq.bI"])
  sigmasq.e1 <- exp(param["lsigmasq.e1"])
  cov.mat <- matrix(c(sigmasq.b + sigmasq.bI + sigmasq.e1, sigmasq.b + sigmasq.bI, 
                      sigmasq.b + sigmasq.bI, sigmasq.b + sigmasq.bI + sigmasq.e1), byrow = T, 
                    ncol = 2)
  result <- zccc(0, cov.mat[1, 1], cov.mat[2, 2], cov.mat[1, 2])
  return(result)
}

######################
#
# Function for computing log(TDI) as a function of model parameters
# Arguments: 
#		param = transformed model parameters
# 		prob = probability
#
ltdi.unlinked.fun <- function(param, prob) {
  # extract the parameters
  alpha <- param["alpha"]
  sigmasq.bI <- exp(param["lsigmasq.bI"])
  sigmasq.e1 <- exp(param["lsigmasq.e1"])
  sigmasq.e2 <- exp(param["lsigmasq.e2"])
  result <- ltdi(alpha, sqrt(2 * sigmasq.bI + sigmasq.e1 + sigmasq.e2), prob)
  return(result)
}

##################
#
# Function for computing Fisher's z-transformation of CCC as a function of model parameters
# Arguments: 
#		param = transformed model parameters
#
zccc.unlinked.fun <- function(param) {
  # extract the parameters
  alpha <- param["alpha"]
  sigmasq.b <- exp(param["lsigmasq.b"])
  sigmasq.bI <- exp(param["lsigmasq.bI"])
  sigmasq.e1 <- exp(param["lsigmasq.e1"])
  sigmasq.e2 <- exp(param["lsigmasq.e2"])
  cov.mat <- matrix(c(sigmasq.b + sigmasq.bI + sigmasq.e1, sigmasq.b, sigmasq.b, 
                      sigmasq.b + sigmasq.bI + sigmasq.e2), byrow = T, ncol = 2)
  result <- zccc(alpha, cov.mat[1, 1], cov.mat[2, 2], cov.mat[1, 2])
  return(result)
}

##################
#
# Function for fitting a mixed model to unlinked repeated measurements data
# Arguments: 
#		dat = a data frame with four columns --- subject (a factor), 
#			  method (a factor with two levels method1 and method2), 
#			rep (replication number, unused) and meas
# Value: A list with following components
#		fit = fitted lme object
#		param.hat = ML estimates of transformed parameters
#		var.hat = ML estimates of variance components and variance of difference
#		bvn.hat = parameters of fitted bivariate normal distribution
#		reliability.hat = ML estimates of reliabilities of the two methods
#


lme.unlinked.repmeas.fit <- function(dat) {
  require(nlme)
  # fit model
  dat <- groupedData(meas ~ method | subject, data = dat, order.groups = FALSE)
  fit <- lme(meas ~ method - 1, data = dat, weights = varIdent(form = ~1 | 
                                                                 method), method = "ML", random = pdBlocked(list(~1, pdIdent(~method - 
                                                                                                                               1))))
  
  # extract estimates
  mu.1.hat <- fixef(fit)["methodmethod1"]
  mu.2.hat <- fixef(fit)["methodmethod2"]
  sigma.hat <- fit$sigma
  
  mat1 <- diag(as.matrix(fit$modelStruct$reStruct[[1]]) * sigma.hat^2)
  sigmasq.b.hat <- mat1["(Intercept)"]
  sigmasq.bI.hat <- mat1["methodmethod1"]
  
  sigma.ratio.hat <- exp(as.numeric(fit$modelStruct$varStruct[1]))
  
  # determine which level of method is taken as reference in variance function
  wt <- attr(fit$modelStruct$varStruct, "weights")
  ref.level <- unique(names(which((wt == 1))))
  
  # get appropriate variance estimates
  if (ref.level == "method1") {
    sigmasq.e1.hat <- sigma.hat^2
    sigmasq.e2.hat <- (sigma.hat * sigma.ratio.hat)^2
  } else {
    sigmasq.e2.hat <- sigma.hat^2
    sigmasq.e1.hat <- (sigma.hat * sigma.ratio.hat)^2
  }
  # cast estimates in desired form
  mu.b.hat <- mu.1.hat
  alpha.hat <- mu.2.hat - mu.1.hat
  
  param.hat <- c(alpha.hat, mu.b.hat, log(c(sigmasq.b.hat, sigmasq.bI.hat, 
                                            sigmasq.e1.hat, sigmasq.e2.hat)))
  names(param.hat) <- c("alpha", "mu.b", "lsigmasq.b", "lsigmasq.bI", "lsigmasq.e1", 
                        "lsigmasq.e2")
  
  # estimated variances
  var.hat <- c(sigmasq.b.hat, sigmasq.bI.hat, sigmasq.e1.hat, sigmasq.e2.hat, 
               2 * sigmasq.bI.hat + sigmasq.e1.hat + sigmasq.e2.hat)
  names(var.hat) <- c("sigmasq.b", "sigmasq.bI", "sigmasq.e1", "sigmasq.e2", 
                      "sigmasq")
  
  # fitted bivariate distribution 
  var.1.hat <- sigmasq.e1.hat + sigmasq.bI.hat + sigmasq.b.hat
  var.2.hat <- sigmasq.e2.hat + sigmasq.bI.hat + sigmasq.b.hat
  corr.hat <- sigmasq.b.hat/sqrt(var.1.hat * var.2.hat)
  bvn.hat <- c(mu.1.hat, mu.2.hat, var.1.hat, var.2.hat, corr.hat)
  names(bvn.hat) <- c("mean1", "mean2", "var1", "var2", "corr")
  
  # reliabilities
  
  rlbty.1.hat <- (sigmasq.b.hat + sigmasq.bI.hat)/var.1.hat
  rlbty.2.hat <- (sigmasq.b.hat + sigmasq.bI.hat)/var.2.hat
  rlbty.hat <- c(rlbty.1.hat, rlbty.2.hat)
  names(rlbty.hat) <- c("method1", "method2")
  
  return(list(fit = fit, param.hat = param.hat, var.hat = var.hat, bvn.hat = bvn.hat, 
              reliability.hat = rlbty.hat))
}

###########################
# Function for making a B-A plot of averaged unlinked repeated measurementds data
# Arguments: 
#		dat = a data frame with four columns --- subject (a factor), 
#			  method (a factor with levels method1 and method2), 
#			  rep (replication number, not used), and meas
# 		ylim = range for y axis
# 		axes.label = a character vector of length two containing 
#			  labels for x and y axes
#

plot.ba.avg <- function(dat, ylim) {
  yy <- unlinked2avg(dat)
  # plot them	
  plot(1, type = "n", xlim = range(dat$meas), ylim = ylim, 
       xlab = paste0("(", colnames(yy)[1], "+", colnames(yy)[2], ")/2"),
       ylab = paste(colnames(yy)[2], "-", colnames(yy)[1]),
       panel.first = grid(lty = "solid"))
  points((yy[, 1] + yy[, 2])/2, yy[, 2] - yy[, 1])
  abline(h = 0)
}
###########################
# Function for making scatterplot of averaged unlinked repeated measurementds data
# Arguments: 
#		dat = a data frame with four columns --- subject (a factor), 
#			  method (a factor with levels method1 and method2), 
#			  rep (replication number, not used), and meas
# 		axes.label = a character vector of length two containing 
#			  labels for x and y axes
#


plot.sp.avg <- function(dat) {
  # get the averages
  yy <- unlinked2avg(dat)
  # plot them
  plot(1, type = "n", xlim = range(dat$meas), ylim = range(dat$meas), 
       xlab = colnames(yy)[1], ylab = colnames(yy)[2], 
       panel.first = grid(lty = "solid"))
  points(yy[, 1], yy[, 2])
  abline(a = 0, b = 1)
}



###########################
#
# Function for plotting B-A plot of linked repeated measurementds
# data by using subject id as plotting symbol
# Arguments: 
#		dat = a data frame with four columns --- subject (a factor), 
#			  method (a factor with levels method1 and method2), 
#			  rep (replication number, not used), and meas
# 		ylim = range for y axis
# 		axes.label = a character vector of length two containing 
#			  labels for x and y axes
#

plot.ba.subid <- function(dat, ylim, axes.labels = NULL) {
  # make an empty plot
  if (is.null(axes.labels)) {
    axes.labels <- c("(method1 + method2)/2", "method2 - method1")
  }
  plot(1, type = "n", xlim = range(dat$meas), ylim = ylim, xlab = axes.labels[1], 
       ylab = axes.labels[2], panel.first = grid(lty = "solid"))
  # split data by subject id
  dat.split <- split(dat, dat$subject)
  for (i in 1:length(dat.split)) {
    # further split data by rep number
    yy <- split(dat.split[[i]], dat.split[[i]]$rep, drop = TRUE)
    yy1 <- numeric(length(yy))
    yy2 <- numeric(length(yy))
    for (j in 1:length(yy)) {
      # extract the measurements by the two methods
      yy1[j] <- yy[[j]][yy[[j]]$method == "method1", "meas"]
      yy2[j] <- yy[[j]][yy[[j]]$method == "method2", "meas"]
    }
    # finally plot
    text((yy1 + yy2)/2, yy2 - yy1, as.character(i), cex = 0.65)
  }
  abline(h = 0)
}
###########################
#
# Function for plotting B-A plot of linked repeated measurementds
# data by joining the points from the same subject
# Arguments: 
#		dat = a data frame with four columns --- subject (a factor), 
#			  method (a factor with levels method1 and method2), 
#			  rep (replication number, not used), and meas
# 		ylim = range for y axis
# 		axes.label = a character vector of length two containing 
#			  labels for x and y axes
#

plot.ba.join <- function(dat, ylim, axes.labels = NULL) {
  # make an empty plot
  if (is.null(axes.labels)) {
    axes.labels <- c("(method1 + method2)/2", "method2 - method1")
  }
  plot(1, type = "n", xlim = range(dat$meas), ylim = ylim, xlab = axes.labels[1], 
       ylab = axes.labels[2], panel.first = grid(lty = "solid"))
  # split data by subject id
  dat.split <- split(dat, dat$subject)
  for (i in 1:length(dat.split)) {
    # further split data by rep number
    yy <- split(dat.split[[i]], dat.split[[i]]$rep, drop = TRUE)
    yy1 <- numeric(length(yy))
    yy2 <- numeric(length(yy))
    for (j in 1:length(yy)) {
      # extract the measurements by the two methods
      yy1[j] <- yy[[j]][yy[[j]]$method == "method1", "meas"]
      yy2[j] <- yy[[j]][yy[[j]]$method == "method2", "meas"]
    }
    # finally plot
    points((yy1 + yy2)/2, yy2 - yy1, type = "o", lty = 2)
  }
  abline(h = 0)
}
###########################
#
# Function for plotting scatterplot of linked repeated measurementds
# data by using subject id as plotting symbol
# Arguments: 
#		dat = a data frame with four columns --- subject (a factor), 
#			  method (a factor with levels method1 and method2), 
#			  rep (replication number, not used), and meas
# 		axes.label = a character vector of length two containing 
#			  labels for x and y axes
#
plot.sp.subid <- function(dat, axes.labels = NULL) {
  # make an empty plot
  if (is.null(axes.labels)) {
    axes.labels <- c("method1", "method2")
  }
  plot(1, type = "n", xlim = range(dat$meas), ylim = range(dat$meas), xlab = axes.labels[1], 
       ylab = axes.labels[2], panel.first = grid(lty = "solid"))
  # split data by subject id
  dat.split <- split(dat, dat$subject)
  for (i in 1:length(dat.split)) {
    # further split data by rep number
    yy <- split(dat.split[[i]], dat.split[[i]]$rep, drop = TRUE)
    yy1 <- numeric(length(yy))
    yy2 <- numeric(length(yy))
    for (j in 1:length(yy)) {
      # extract the measurements by the two methods
      yy1[j] <- yy[[j]][yy[[j]]$method == "method1", "meas"]
      yy2[j] <- yy[[j]][yy[[j]]$method == "method2", "meas"]
    }
    # finally plot
    text(yy1, yy2, as.character(i), cex = 0.65)
  }
  abline(a = 0, b = 1)
}

###########################
#
# Function for plotting scatterplot of linked repeated measurementds
# data by joining the points from the same subject
# Arguments: 
#		dat = a data frame with four columns --- subject (a factor), 
#			  method (a factor with levels method1 and method2), 
#			  rep (replication number, not used), and meas
# 		axes.label = a character vector of length two containing 
#			  labels for x and y axes
#
plot.sp.join <- function(dat, axes.labels = NULL) {
  # make an empty plot
  if (is.null(axes.labels)) {
    axes.labels <- c("method1", "method2")
  }
  plot(1, type = "n", xlim = range(dat$meas), ylim = range(dat$meas), xlab = axes.labels[1], 
       ylab = axes.labels[2], panel.first = grid(lty = "solid"))
  # split data by subject id
  dat.split <- split(dat, dat$subject)
  for (i in 1:length(dat.split)) {
    # further split data by rep number
    yy <- split(dat.split[[i]], dat.split[[i]]$rep, drop = TRUE)
    yy1 <- numeric(length(yy))
    yy2 <- numeric(length(yy))
    for (j in 1:length(yy)) {
      # extract the measurements by the two methods
      yy1[j] <- yy[[j]][yy[[j]]$method == "method1", "meas"]
      yy2[j] <- yy[[j]][yy[[j]]$method == "method2", "meas"]
    }
    # finally plot
    points(yy1, yy2, type = "o", lty = 2)
  }
  abline(a = 0, b = 1)
}

###########################
#
# Function for computing averages of measurements from each method on 
# every subject from unlinked repeated measurements data
# Arguments: 
#		dat = a data frame with four columns --- subject (a factor), 
#			  method (a factor), 
#			  rep (replication number, not used), and meas
# Value:
#	A n x J dataframe of averages where n is the number of subjects and 
# 	J is the number of methods
unlinked2avg <- function(dat) {
  newdat <- with(dat, by(meas, list(subject, method), mean))
  n <- length(attr(newdat, "dimnames")[[1]])
  J <- length(attr(newdat, "dimnames")[[2]])
  return(newdat[1:n, 1:J])
}

###########################
#
# Function for obtaining one randomly select measurement from each method on 
# every subject from unlinked repeated measurements data
# Arguments: 
#		dat = a data frame with four columns --- subject (a factor), 
#			  method (a factor), 
#			  rep (replication number, not used), and meas
# Value:
#	A n x J dataframe of measurements where n is the number of subjects and 
# 	J is the number of methods
unlinked2paired <- function(dat) {
  resample <- function(x, ...) x[sample.int(length(x), ...)] # needed to avoid problems with sample of size 
  newdat <- with(dat, by(meas, list(subject, method), resample, size = 1))
  n <- length(attr(newdat, "dimnames")[[1]])
  J <- length(attr(newdat, "dimnames")[[2]])
  return(newdat[1:n, 1:J])
}
###########################


####################################################################################
# chapter 4 #
####################################################################################

##################
#
# Function for computing estimates and their SEs for various parameters and functions
# Arguments: 
# 		ymat = data matrix with two columns: method1 and method2
#		param.hat = vector of ML estimates of transformed parameters
#		fit.method = lme (mixed-effects model) or bvn (BVN model)
# 		fit.flag = yes (bvn fit is same as lme fit) or no (the two fits are not the same)
# Value: 
#		standard errors of estimates of parameters and their confidence intervals
#		estimate, se and ucb of log(TDI), and estimate and ucb of TDI
#		estimate, se and lcb of z(CCC), and estimate and lcb of CCC
#
conf.measures <- function(ymat, param.hat, fit.method, bvn.fit.flag, conf.level = 0.95, 
                          prob = 0.9) {
  require(numDeriv)
  
  if (!((fit.method == "lme") || (fit.method == "bvn"))) {
    stop("invalid 'fit.method' argument")
  }
  
  if (!((bvn.fit.flag == "yes") || (bvn.fit.flag == "no"))) {
    stop("invalid 'bvn.fit.flag' argument")
  }
  
  z.twosided <- qnorm((1 + conf.level)/2)
  z.onesided <- qnorm(conf.level)
  
  # get Hessian and invert it
  hess <- -hessian(loglik.paired.data, param.hat, ymat = ymat, fit.method = fit.method)
  inv.hess <- solve(hess)
  rownames(inv.hess) <- names(param.hat)
  colnames(inv.hess) <- names(param.hat)
  
  # SEs of param.hat and confidence intevals
  se.param.hat <- sqrt(diag(inv.hess))
  ci.param <- cbind(param.hat - z.twosided * se.param.hat, param.hat + z.twosided * 
                      se.param.hat)
  results1 <- cbind(param.hat, se.param.hat, ci.param)
  rownames(results1) <- names(param.hat)
  colnames(results1) <- c("estimate", "se", "lcl", "ucl")
  
  if (fit.method == "lme") {
    # estimate and confidence interval for precision ratio
    llambda.hat <- param.hat["lsigmasq.e1"] - param.hat["lsigmasq.e2"]
    se.llambda.hat <- sqrt(inv.hess["lsigmasq.e1", "lsigmasq.e1"] + inv.hess["lsigmasq.e2", 
                                                                             "lsigmasq.e2"] - 2 * inv.hess["lsigmasq.e1", "lsigmasq.e2"])
    ci.llambda <- llambda.hat + c(-1, 1) * z.twosided * se.llambda.hat
    
    # estimate and confidence interval for variance ratio
    lvarratio.hat <- lvar.ratio.fun(param.hat)
    lvarratio.grad <- grad(lvar.ratio.fun, param.hat)
    se.lvarratio.hat <- sqrt(drop(t(lvarratio.grad) %*% inv.hess %*% lvarratio.grad))
    ci.lvarratio <- lvarratio.hat + c(-1, 1) * z.twosided * se.lvarratio.hat
    
    results2 <- rbind(c(llambda.hat, se.llambda.hat, ci.llambda, exp(llambda.hat), 
                        exp(ci.llambda)), c(lvarratio.hat, se.lvarratio.hat, ci.lvarratio, 
                                            exp(lvarratio.hat), exp(ci.lvarratio)))
    rownames(results2) <- c("lambda", "varratio")
    colnames(results2) <- c("lest", "lest.se", "lcl.lpar", "ucl.lpar", "est", 
                            "lcl.par", "ucl.par")
  }
  
  if (fit.method == "bvn") {
    # estimate and confidence interval for mean difference 
    meandiff.hat <- param.hat["mu.2"] - param.hat["mu.1"]
    se.meandiff.hat <- sqrt(inv.hess["mu.2", "mu.2"] + inv.hess["mu.1", "mu.1"] - 
                              2 * inv.hess["mu.2", "mu.1"])
    ci.meandiff <- meandiff.hat + c(-1, 1) * z.twosided * se.meandiff.hat
    meandiff <- c(meandiff.hat, se.meandiff.hat, ci.meandiff)
    names(meandiff) <- c("est", "se", "lcl", "ucl")
    
    # estimate and confidence interval for variance ratio
    
    lvarratio.hat <- param.hat["lsigmasq.2"] - param.hat["lsigmasq.1"]
    se.lvarratio.hat <- sqrt(inv.hess["lsigmasq.2", "lsigmasq.2"] + inv.hess["lsigmasq.1", 
                                                                             "lsigmasq.1"] - 2 * inv.hess["lsigmasq.2", "lsigmasq.1"])
    ci.lvarratio <- lvarratio.hat + c(-1, 1) * z.twosided * se.lvarratio.hat
    
    # estimate and confidence interval for precision ratio
    
    if (bvn.fit.flag == "yes") {
      llambda.hat <- lprec.ratio.fun(param.hat)
      llambda.grad <- grad(lprec.ratio.fun, param.hat)
      se.llambda.hat <- sqrt(drop(t(llambda.grad) %*% inv.hess %*% llambda.grad))
      ci.llambda <- llambda.hat + c(-1, 1) * z.twosided * se.llambda.hat
      
      results2 <- rbind(c(llambda.hat, se.llambda.hat, ci.llambda, exp(llambda.hat), 
                          exp(ci.llambda)), c(lvarratio.hat, se.lvarratio.hat, ci.lvarratio, 
                                              exp(lvarratio.hat), exp(ci.lvarratio)))
      rownames(results2) <- c("lambda", "varratio")
      colnames(results2) <- c("lest", "lest.se", "lcl.lpar", "ucl.lpar", 
                              "est", "lcl.par", "ucl.par")
    } else {
      results2 <- c(lvarratio.hat, se.lvarratio.hat, ci.lvarratio, exp(lvarratio.hat), 
                    exp(ci.lvarratio))
      names(results2) <- c("lvarratio.est", "lvarratio.est.se", "lcl.lvarratio", 
                           "ucl.lvarratio", "varratio.est", "lcl.varratio", "ucl.varratio")
    }
    
    results2 <- list(meandiff = meandiff, varpar = results2)
  }
  
  # estimate and upper confidence bound of TDI
  ltdi.hat <- ltdi.fun(param.hat, prob, fit.method)
  ltdi.grad <- grad(ltdi.fun, param.hat, prob = prob, fit.method = fit.method)
  se.ltdi.hat <- sqrt(drop(t(ltdi.grad) %*% inv.hess %*% ltdi.grad))
  ucb.ltdi <- ltdi.hat + z.onesided * se.ltdi.hat
  results3 <- c(ltdi.hat, se.ltdi.hat, ucb.ltdi, exp(ltdi.hat), exp(ucb.ltdi))
  names(results3) <- c("ltdi.hat", "se.ltdi.hat", "ucb.ltdi", "tdi.hat", "ucb.tdi")
  
  # estimate and lower confidence bound of CCC
  zccc.hat <- zccc.fun(param.hat, fit.method)
  zccc.grad <- grad(zccc.fun, param.hat, fit.method = fit.method)
  se.zccc.hat <- sqrt(drop(t(zccc.grad) %*% inv.hess %*% zccc.grad))
  lcb.zccc <- zccc.hat - z.onesided * se.zccc.hat
  results4 <- c(zccc.hat, se.zccc.hat, lcb.zccc, tanh(zccc.hat), tanh(lcb.zccc))
  names(results4) <- c("zccc.hat", "se.zccc.hat", "lcb.zccc", "ccc.hat", "lcb.ccc")
  
  results <- list(param = results1, similarity = results2, tdi = results3, 
                  ccc = results4)
  return(results)
}
##################
#
# Function for computing log-precision ratio as a function of BVN model parameters
# Arguments: 
#		bvn.param = transformed model parameters
#
lprec.ratio.fun <- function(bvn.param) {
  sigmasq.1 <- exp(bvn.param["lsigmasq.1"])
  sigmasq.2 <- exp(bvn.param["lsigmasq.2"])
  rho <- tanh(bvn.param["zrho"])
  sigmasq.b <- rho * sqrt(sigmasq.1 * sigmasq.2)
  sigmasq.e1 <- sigmasq.1 - sigmasq.b
  sigmasq.e2 <- sigmasq.2 - sigmasq.b
  return(log(sigmasq.e1) - log(sigmasq.e2))
}

##################
#
# Function for computing log-variance ratio as a function of mixed-effects model parameters
# Arguments: 
#		lme.param = transformed model parameters
#
lvar.ratio.fun <- function(lme.param) {
  sigmasq.1 <- exp(lme.param["lsigmasq.e1"]) + exp(lme.param["lsigmasq.b"])
  sigmasq.2 <- exp(lme.param["lsigmasq.e2"]) + exp(lme.param["lsigmasq.b"])
  return(log(sigmasq.2) - log(sigmasq.1))
}


##################
#
# Function for computing Fisher's z-transformation of CCC as a function of model parameters
# Arguments: 
#		param = transformed model parameters
#		fit.method = lme (mixed-effects model) or bvn (BVN model)
#
zccc.fun <- function(param, fit.method) {
  if (fit.method == "lme") {
    # extract the parameters
    alpha <- param["alpha"]
    sigmasq.b <- exp(param["lsigmasq.b"])
    sigmasq.e1 <- exp(param["lsigmasq.e1"])
    sigmasq.e2 <- exp(param["lsigmasq.e2"])
    cov.mat <- matrix(c(sigmasq.b + sigmasq.e1, sigmasq.b, sigmasq.b, sigmasq.b + 
                          sigmasq.e2), byrow = T, ncol = 2)
    result <- zccc(alpha, cov.mat[1, 1], cov.mat[2, 2], cov.mat[1, 2])
  }
  
  if (fit.method == "bvn") {
    # extract the parameters
    mu.1 <- param["mu.1"]
    mu.2 <- param["mu.2"]
    sigmasq.1 <- exp(param["lsigmasq.1"])
    sigmasq.2 <- exp(param["lsigmasq.2"])
    rho <- tanh(param["zrho"])
    cov.mat <- matrix(c(sigmasq.1, rho * sqrt(sigmasq.1 * sigmasq.2), rho * 
                          sqrt(sigmasq.1 * sigmasq.2), sigmasq.2), byrow = T, ncol = 2)
    result <- zccc(mu.2 - mu.1, cov.mat[1, 1], cov.mat[2, 2], cov.mat[1, 
                                                                      2])
  }
  return(result)
}

####################
#
# Function for computing Fisher's z-transformation of CCC
# Arguments: 
#		dmean = mean of difference
#		var1, var2, cov12 = variances and covariance of measurements
#
zccc <- function(dmean, var1, var2, cov12) {
  ccc <- (2 * cov12)/(var1 + var2 + dmean^2)
  return(atanh(ccc))
}

######################
#
# Function for computing log(TDI) as a function of model parameters
# Arguments: 
#		param = transformed model parameters
# 		prob = probability
#		fit.method = lme (mixed-effects model) or bvn (BVN model)
#
ltdi.fun <- function(param, prob, fit.method) {
  if (fit.method == "lme") {
    # extract the parameters
    alpha <- param["alpha"]
    sigmasq.e1 <- exp(param["lsigmasq.e1"])
    sigmasq.e2 <- exp(param["lsigmasq.e2"])
    result <- ltdi(alpha, sqrt(sigmasq.e1 + sigmasq.e2), prob)
  }
  if (fit.method == "bvn") {
    # extract the parameters
    mu.1 <- param["mu.1"]
    mu.2 <- param["mu.2"]
    sigmasq.1 <- exp(param["lsigmasq.1"])
    sigmasq.2 <- exp(param["lsigmasq.2"])
    rho <- tanh(param["zrho"])
    result <- ltdi(mu.2 - mu.1, sqrt(sigmasq.1 + sigmasq.2 - 2 * rho * sqrt(sigmasq.1 * 
                                                                              sigmasq.2)), prob)
  }
  return(result)
}

########################
#
# Function for computing log(TDI)
# Arguments: 
#		dmean, dsd = mean and sd of difference
# 		prob = probability
# 		llim, ulim = range in which the true value lies 
#
ltdi <- function(dmean, dsd, prob, llim = 0.001, ulim = 200) {
  # set up the equation whose solution is tdi
  cdf.abs.diff <- function(q, dmean, dsd, prob) {
    result <- pnorm((q - dmean)/dsd) - pnorm((-q - dmean)/dsd)
    return(result - prob)
  }
  result <- uniroot(cdf.abs.diff, c(llim, ulim), extendInt="upX", dmean = dmean, dsd = dsd, 
                    prob = prob)$root
  return(log(result))
}
##################
#
# Function for fitting a BVN model to the paired data
# Argument: 
#		ymat = data matrix with two columns: method1 and method2
# Value: A list with following components
#		param.hat = ML estimates of transformed parameters
#		bvn.hat = parameters of fitted bivariate normal distribution
#		fit.flag = "yes" if bvn.fit gives mixed-effects model fit as well,
#					"no", otherwise
#		param.hat.lme = ML estimates of parameters of mixed-effects model, 
#					fit.flag is "yes"									
#
bvn.fit <- function(ymat) {
  require(mvtnorm)
  n <- nrow(ymat)
  # ML estimates of bvn parameters
  mu.hat <- colMeans(ymat)
  Sigma.hat <- var(ymat) * ((n - 1)/n)
  mu.1.hat <- mu.hat[1]
  mu.2.hat <- mu.hat[2]
  sigmasq.1.hat <- Sigma.hat[1, 1]
  sigmasq.2.hat <- Sigma.hat[2, 2]
  rho.hat <- cov2cor(Sigma.hat)[1, 2]
  
  # cast estimates in desired form
  param.hat <- c(mu.hat, log(c(sigmasq.1.hat, sigmasq.2.hat)), atanh(rho.hat))
  names(param.hat) <- c("mu.1", "mu.2", "lsigmasq.1", "lsigmasq.2", "zrho")
  
  # fitted bivariate distribution 
  bvn.hat <- c(mu.hat, sigmasq.1.hat, sigmasq.2.hat, rho.hat)
  names(bvn.hat) <- c("mean1", "mean2", "var1", "var2", "corr")
  
  # fitted distribution of difference
  mean.D.hat <- mu.2.hat - mu.1.hat
  var.D.hat <- sigmasq.1.hat + sigmasq.2.hat - 2*rho.hat*sqrt(sigmasq.1.hat*sigmasq.2.hat)
  D.hat <- c(mean.D.hat, var.D.hat)
  names(D.hat) <- c("mean", "var")
  
  # limits of agreement
  loa <- mean.D.hat + c(-1,1)*1.96*sqrt(var.D.hat)
  
  
  # see if bvn fit gives mixed-effects model fit as well
  sigmasq.b.hat <- Sigma.hat[1, 2]
  sigmasq.e1.hat <- sigmasq.1.hat - sigmasq.b.hat
  sigmasq.e2.hat <- sigmasq.2.hat - sigmasq.b.hat
  
  if (((sigmasq.b.hat > 0) && (sigmasq.e1.hat > 0) && (sigmasq.e2.hat > 0))) {
    fit.flag <- "yes"
    param.hat.lme <- c(mu.2.hat - mu.1.hat, mu.1.hat, log(c(sigmasq.b.hat, 
                                                            sigmasq.e1.hat, sigmasq.e2.hat)))
    names(param.hat.lme) <- c("alpha", "mu.b", "lsigmasq.b", "lsigmasq.e1", 
                              "lsigmasq.e2")
  } else {
    fit.flag <- "no"
    param.hat.lme <- NULL
  }
  
  # log-likelihood of fitted model
  loglik <- sum(dmvnorm(ymat, mean=mu.hat, sigma=Sigma.hat, log=TRUE))
  aic <- -2*loglik + 2*length(param.hat)
  bic <- -2*loglik + log(2*n)*length(param.hat)		
  
  
  return(list(param.hat = param.hat, bvn.hat = bvn.hat, loglik=loglik,
              aic=aic, bic=bic, D.hat=D.hat, loa=loa,
              fit.flag = fit.flag, param.hat.lme = param.hat.lme))
}

##################
#
# Function for fitting a mixed model to the paired data
# Arguments: 
#		dat = a data frame with three columns --- subject (a factor), 
#			  method (a factor with levels method1 and method2), and meas
# Value: A list with following components
#		fit = fitted lme object
#		param.hat = ML estimates of transformed parameters
#		var.hat = ML estimates of variance components and variance of difference
#		bvn.hat = parameters of fitted bivariate normal distribution
#		reliability.hat = ML estimates of reliabilities of the two methods
#
lme.fit <- function(dat) {
  # fit model
  require(nlme)
  dat <- groupedData(meas ~ method | subject, data = dat, order.groups = FALSE)
  fit <- lme(meas ~ method - 1, data = dat, weights = varIdent(form = ~1 | 
                                                                 method), method = "ML", random = ~1)
  
  # extract estimates
  mu.1.hat <- fixef(fit)[1]
  mu.2.hat <- fixef(fit)[2]
  sigma.hat <- fit$sigma
  sigmasq.b.hat <- exp(2 * as.numeric(fit$modelStruct$reStruct[[1]])) * sigma.hat^2
  sigma.ratio.hat <- exp(as.numeric(fit$modelStruct$varStruct[1]))
  
  # determine which level of method is taken as reference in variance function
  wt <- attr(fit$modelStruct$varStruct, "weights")
  ref.level <- unique(names(which((wt == 1))))
  
  # get appropriate variance estimates
  if (ref.level == "method1") {
    sigmasq.e1.hat <- sigma.hat^2
    sigmasq.e2.hat <- (sigma.hat * sigma.ratio.hat)^2
  } else {
    sigmasq.e2.hat <- sigma.hat^2
    sigmasq.e1.hat <- (sigma.hat * sigma.ratio.hat)^2
  }
  # cast estimates in desired form
  mu.b.hat <- mu.1.hat
  alpha.hat <- mu.2.hat - mu.1.hat
  
  param.hat <- c(alpha.hat, mu.b.hat, log(c(sigmasq.b.hat, sigmasq.e1.hat, 
                                            sigmasq.e2.hat)))
  names(param.hat) <- c("alpha", "mu.b", "lsigmasq.b", "lsigmasq.e1", "lsigmasq.e2")
  
  # estimated variances
  var.hat <- c(sigmasq.b.hat, sigmasq.e1.hat, sigmasq.e2.hat, sigmasq.e1.hat + 
                 sigmasq.e2.hat)
  names(var.hat) <- c("sigmasq.b", "sigmasq.e1", "sigmasq.e2", "sigmasq")
  
  # fitted bivariate distribution 
  var.1.hat <- sigmasq.e1.hat + sigmasq.b.hat
  var.2.hat <- sigmasq.e2.hat + sigmasq.b.hat
  corr.hat <- sigmasq.b.hat/sqrt(var.1.hat * var.2.hat)
  bvn.hat <- c(mu.1.hat, mu.2.hat, var.1.hat, var.2.hat, corr.hat)
  names(bvn.hat) <- c("mean1", "mean2", "var1", "var2", "corr")
  
  # reliabilities
  
  rlbty.1.hat <- sigmasq.b.hat/var.1.hat
  rlbty.2.hat <- sigmasq.b.hat/var.2.hat
  rlbty.hat <- c(rlbty.1.hat, rlbty.2.hat)
  names(rlbty.hat) <- c("method1", "method2")
  
  return(list(fit = fit, param.hat = param.hat, var.hat = var.hat, bvn.hat = bvn.hat, 
              reliability.hat = rlbty.hat))
}

######################
#
# Function for computing log(L) using paired data
# Arguments: 
#		param = vector of transformed model parameters
# 		ymat = data matrix with two columns: method1 and method2
#		fit.method = lme (mixed-effects model) or bvn (BVN model)
#


loglik.paired.data <- function(param, ymat, fit.method) {
  require(mvtnorm)
  
  if (!((fit.method == "lme") || (fit.method == "bvn"))) {
    stop("invalid 'fit.method' argument")
  }
  
  if (fit.method == "lme") {
    # extract the parameters
    alpha <- param["alpha"]
    mu.b <- param["mu.b"]
    sigmasq.b <- exp(param["lsigmasq.b"])
    sigmasq.e1 <- exp(param["lsigmasq.e1"])
    sigmasq.e2 <- exp(param["lsigmasq.e2"])
    
    # get mean vector and covariance matrix
    mean.vec <- c(mu.b, alpha + mu.b)
    cov.mat <- matrix(c(sigmasq.b + sigmasq.e1, sigmasq.b, sigmasq.b, sigmasq.b + 
                          sigmasq.e2), byrow = T, ncol = 2)
  }
  
  if (fit.method == "bvn") {
    # extract the parameters
    mu.1 <- param["mu.1"]
    mu.2 <- param["mu.2"]
    sigmasq.1 <- exp(param["lsigmasq.1"])
    sigmasq.2 <- exp(param["lsigmasq.2"])
    rho <- tanh(param["zrho"])
    
    # get mean vector and covariance matrix
    mean.vec <- c(mu.1, mu.2)
    cov.mat <- matrix(c(sigmasq.1, rho * sqrt(sigmasq.1 * sigmasq.2), rho * 
                          sqrt(sigmasq.1 * sigmasq.2), sigmasq.2), byrow = T, ncol = 2)
  }
  
  # get log-likelihood
  loglik <- sum(apply(ymat, MAR = 1, FUN = dmvnorm, mean = mean.vec, sigma = cov.mat, 
                      log = T))
  return(loglik)
}
###########################
#
# Function for making a scatterplot and a B-A plot
# Arguments: 
# 		ymat = data matrix with two columns: method1 and method2
#
plot.sba <- function(ymat) {
  y1 <- ymat[, 1]
  y2 <- ymat[, 2]
  par(mfrow = c(1, 2))
  
  # scatter plot
  plot(y1, y2, xlab = "method1", ylab = "method2", xlim = c(min(y1, y2), max(y1, 
                                                                             y2)), ylim = c(min(y1, y2), max(y1, y2)))
  abline(a = 0, b = 1)
  title(main = "(a)")
  
  # B-A plot
  plot((y2 + y1)/2, y2 - y1, ylab = "method2 - method1", xlab = "(method1 + method2)/2")
  abline(h = 0)
  title(main = "(b)")
  par(mfrow = c(1, 1))
}
###########################
#
# Function for making a Trellis plot
# Arguments: 
#		dat = a data frame with three columns --- subject (a factor), 
#			  method (a factor with levels method1 and method2), and meas
#
plot.trellis <- function(dat) {
  require(lattice)
  require(nlme)
  dat <- groupedData(meas ~ method | subject, data = dat, order.groups = TRUE)
  plot(dat, xlab = "measurement", ylab = "subject")
}
###########################
#
# Function for converting a data frame with 3 columns into a matrix with 2 columns
# Arguments: 
#		dat = a data frame with three columns --- subject (a factor), 
#			  method (a factor with levels method1 and method2), and meas
#
make.2col <- function(dat) {
  y1 <- subset(dat, subset = (method == "method1"), select = meas, drop = TRUE)
  y2 <- subset(dat, subset = (method == "method2"), select = meas, drop = TRUE)
  result <- data.frame(y1, y2)
  names(result) <- c("method1", "method2")
  return(result)
}
###########################
#
# Function for converting a matrix with 2 columns into a data frame with 3 columns
# Arguments: 
#		ymat = a data frame with three columns --- subject (a factor), 
#			  method (a factor with levels method1 and method2), and meas
#
make.3col <- function(ymat, response="meas") {
  dat <- data.frame(as.factor(rep(1:nrow(ymat), 2)), as.factor(rep(colnames(ymat), 
                                                                   each = nrow(ymat))), c(ymat[, 1], ymat[, 2]))
  names(dat) <- c("subject", "method", response)
  return(dat)
}

################################################################################
# end of MAMMA.R
################################################################################