#######################################################################################################################
## var.wt() - weighted variance
##
## As all of the x_i are drawn from the same distribution and the integer weights w_i indicate the number of occurrences
## ("repeat") of the observations (ratios of sample/control), the unbiased estimator of the weighted population variance
## is given by:
##
##  s^2 = \frac {\sum_{i=1}^N w_i \left(x_i - \mu^*\right)^2 } {\sum_{i=1}^n w_i - 1}
##
var.wt <- function(x, w, mu=NULL, na.rm = TRUE) {
  if (na.rm) {
    w <- w[i <- !is.na(x)]
    x <- x[i]
  }
  if (is.null(mu)) {
    mu <- weighted.mean(x, w, na.rm=na.rm)
  }  
  return( (sum(w*(x-mu)^2)) / (sum(w)-1) )
}

#######################################################################################################################
## weighted.t.test(x, w, mu, conf.level, alternative, na.rm) - weighted t-test 
##
## This will use weighted mean (x.w) and weighted variance (var.w) to 
## calculate a regular t-test.
##
weighted.t.test <- function(x, w, mu, conf.level = 0.95, alternative="two.sided", na.rm=TRUE) {
  
  if(!missing(conf.level) &
      (length(conf.level) != 1 || !is.finite(conf.level) ||
        conf.level < 0 || conf.level > 1))
    stop("'conf.level' must be a single number between 0 and 1")
  
  if (na.rm) {
    w <- w[i <- !is.na(x)]
    x <- x[i]
  }
  
# to achieve consistent behavior in loops, return NA-structure in case of complete missings
  if (sum(is.na(x)) == length(x)) return(list(estimate=NA, se=NA, sd=NA, conf.int=NA, statistic=NA, df=NA, p.value=NA))
  
# if only one value is present: this is the best estimate, no significance test provided
  if (sum(!is.na(x)) == 1) {
    warning("Warning weighted.t.test: only one value provided; this value is returned without test of significance!", call.=FALSE)
    return(list(estimate=x[which(!is.na(x))], se=NA, sd=NA, conf.int=NA, statistic=NA, df=NA, p.value=NA))
  }
  
  x.w   <- weighted.mean(x,w, na.rm=na.rm)
  var.w <- var.wt(x,w, x.w, na.rm=na.rm) # new version of var.w requires x.w
  n     <- length(x)
  ess   <- (sum(w)^2) / sum(w^2)           # this is the effective sample size
  t     <- sqrt(ess)*((x.w-mu)/sqrt(var.w)) 
  se    <- sqrt(var.w)/sqrt(ess)            
  
  if (alternative == "less") {
    pval <- pt(t, ess)
    cint <- c(-Inf, x.w + se*qt(conf.level, ess) )
  }
  else if (alternative == "greater") {
    pval <- pt(t, ess, lower.tail = FALSE)
    cint <- c(x.w - se * qt(conf.level, ess), Inf)
  }
  else {
    pval <- 2 * pt(-abs(t), ess)
    alpha <- 1 - conf.level
    cint <- x.w + se * qt(1 - alpha/2, ess)*c(-1,1)
  }
  
  names(t) <- "t"
  return(list(estimate=x.w, se=se, conf.int=cint, statistic=t, df=ess, p.value=pval, sd=sqrt(var.w)))
}