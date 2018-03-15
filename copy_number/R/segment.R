# segment will try to split a region into parts with homogenous log2 values in each segment
segment <- function(
    dataset,                              ## table with data, will use log2
    refcopy  = 2,                         ## the copy number of the reference
    cutoff   = config.segment.diff.pval,  ## cutoff factor to call segments
    maxscore = config.score.max,
    ## these will be used during recursion to adjust boundaries
    result   = NA, 
    lower    = 0, 
    upper    = Inf
) {
	
	##  message(sprintf("called with lower=%.0f, upper=%.0f", lower, upper))
	
	# subset to usable data first, only check only the defined region
	data <- subset(
        dataset, !is.na(dataset$log2) 
            & is.finite(dataset$log2) 
            & dataset$position >= lower 
            & dataset$position <= upper
  )
  rows <- nrow(data)
  
  # preset some variables
  best.breakpoint <- NA  
  best.var        <- NA
	best.distance   <- NA
  best.pval       <- 1
  
	# basic stats
	data.sd     <- sd(data$log2)
  data.var    <- var(data$log2)
	data.mean   <- mean(data$log2)
	data.median <- median(data$log2)
	
  # get confidence intervals for the distribution
  wtest <- weighted.t.test(
      data$log2,
      w          = data$weight,
      mu         = 0, #data.mean,
      conf.level = config.conf_level
  )
  
#  pvaltest <- weighted.t.test(
#      data$log2,
#      w          = data$weight,
#      mu         = 0,
#      conf.level = config.conf_level
#  )
#  
  
  # data frame row for current range
	current <- data.frame(
        left            = min(data$position), 
        right           = max(data$position), 
        total_amplicons = rows,
        cn              = refcopy * 2 ^ wtest$estimate,
        min             = refcopy * 2 ^ wtest$conf.int[1],
        max             = refcopy * 2 ^ wtest$conf.int[2],
        se              = wtest$se,
        p.val           = wtest$p.val,
        score           = round(min(maxscore ,-5 * log10( wtest$p.val )))
#        is.changed      = ifelse(-5 * log10(wtest$p.val) > minscore, TRUE, FALSE)
  )
	current <- current
#  print(summary(current))
  
	## message("Searching for breakpoint")
	indexes <- c(config.segment.min.amplicons : (rows - config.segment.min.amplicons))

  #for (breakpoint in data$position) {
	for (index in indexes) {
		breakpoint <- data$position[index]

    # // split the data into two parts
    #
    # left part
		left <- subset(data, data$position < breakpoint)
    left <- left
    if (config.segment.rm.out) {
      ## left$outlier  <- flag_outliers(left$log2, loop = TRUE) ## old
      left$outlier  <- flag_outliers(left$log2)
      left          <- subset(left, left$outlier < 1)
    }
    # ensure there are at least ten values at either side
    if (nrow(left)  < config.segment.min.amplicons) next

    # right part
		right <- subset(data, data$position  >= breakpoint)
    right <- right
    if (config.segment.rm.out) {
      ## right$outlier <- flag_outliers(right$log2, loop = TRUE) ## old
      right$outlier <- flag_outliers(right$log2)
      right         <- subset(right, right$outlier < 1)
    }
    # ensure there are at least ten values at either side
    if (nrow(right) < config.segment.min.amplicons) next

    # means
    left.mean <- mean(left$log2)
    left.var  <- var(left$log2)
    
    right.mean <- mean(right$log2)
    right.var <- var(right$log2)
    
    # left test
    left.test <- weighted.t.test(
        left$log2,
        w          = left$weight,
        mu         = right.mean,
        conf.level = config.conf_level
    )
    
    # right test
    right.test <- weighted.t.test(
        right$log2,
        w          = right$weight,
        mu         = left.mean,
        conf.level = config.conf_level
    )
    
    separate <- ifelse(left.test$p.val <= cutoff & right.test$p.val <= cutoff, TRUE, FALSE)
    
#        (
#        (left.test$conf.int[1] > right.test$conf.int[2] )
#      | (left.test$conf.int[2] < right.test$conf.int[1] )
#      ) & (left.test$p.val <= cutoff & right.test$p.val <= cutoff),
#        TRUE, 
#        FALSE
#    )
#
    sum.var  <- left.var + right.var
    distance <- abs(left.mean - right.mean)
    
    if(  separate & 
        (sum.var < data.var) & #distance > 0.5 & #data.sd &
        (sum.var < best.var | is.na(best.var)[1]) ) {
        best.var        <- sum.var
        best.distance   <- distance
        best.breakpoint <- breakpoint
        best.pval       <- left.test$p.val + right.test$p.val
#        print(
#            sprintf("break at: %.0f  L (V) : %.3f (%.3f) / R (V) : %.3f (%.3f). VAR = %.3f; DIST: %.3f PVALS: L = %.5e R = %.5e", 
#                breakpoint, 
#                left.mean, left.var,
#                right.mean, right.var,
#                sum.var, distance,
#                left.test$p.val, right.test$p.val
#            )
#        )
        
		}
	}
	
	if (is.na(best.breakpoint)[1]) {
		## message("no breakpoint found, consider stable")
		if (!is.na(result)[1]) {
      # current$seq.pval <- best.pval
			result<-rbind(result, current)
		}
		else {
			result<-current
		}
		if (upper < max(dataset$position)) {
			return(segment(dataset, result=result, refcopy=refcopy, cutoff=cutoff, lower=upper, upper=Inf))
		}
		return(result)
	}
	## message(best.breakpoint)
	return(segment(dataset, result=result, refcopy=refcopy, cutoff=cutoff, lower=min(data$position), upper=best.breakpoint))
}
