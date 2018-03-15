###############################################################################
# build_model
#
# This function performs a quantile regression using parameters (tau) 
# that will include one stard-deviation of the data (sigma).
#
# input:      ratio     : log2 ratios (x)
#             weight    : weights to use (y)
# 
# output:     model coefficients for q1, q2 and q3 (slope, intercept)
# 
build_model <- function(ratio, weight) {
  # do a quantile regression to find a scaling function for medians
  qreg  <- rq(ratio ~ weight, tau=c(0.15865, 0.5, 0.84135))
  return(c(coef(qreg)))  
}

###############################################################################
# rescale 
#
# input:      x         : row with data
#             col       : column in row containing the value
#             scaling   : factor (slope) of the function
#             intercept : intercept of the scaling function
# 
# output:     rescaled values (array)
# 
# this function is used to scale a value with a dynamic factor, not a constant
rescale <- function(x, col, scaling, intercept) {
  c( as.numeric(x[col]) * 2 ^ ( scaling * log2( as.numeric( x[col] ) ) + intercept) )
}

###############################################################################
# stderr
#
# input:      x         : vector of raw data
stderr <- function(x) sqrt(var(x)/length(x))


###############################################################################
# max_likely_cn
#
# input:      table     : table with amplicon data
max_likely_cn <- function(table, model) {
  
  # start with undefined results
  ml <- data.frame()
  rc <- table$refcopy[1]
  # loop through guesses
  for (guess in (1 : config.plots.max.cn)) {
    probs <- calc_model_pvalue(model, table$log2, table$weight, mu=log2(guess / rc))
    total <- sum(probs, na.rm = TRUE)
    ml <- rbind(ml, data.frame(cn = guess, prob = total))
  }
  ml <- ml[order(-ml[,2]),]
  return(list(
        most_likely = ml[1,1], 
        probability = ml[1,2] / (ml[1,2] + ml[2,2]) )
  )
  # return a named list of results
#  return(list(most_likely=ml.call, probability=ml.sum))
}


###############################################################################
# remove_outliers
#
# input:      x         : raw data
#
# output:     array of values with NA instead of outliers
#
remove_outliers <- function(x, na.rm = TRUE, range = 1.5) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm)
  H <- range * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

###############################################################################
# pull_in_range
#
# input:      x         : raw data
#
# output:     array of values with 1 (out of range) or
#
pull_in_range <- function(
    table, 
    min = config.log2.min, 
    max = config.log2.max,
    mincov = config.min_cov_for_call
  ) {
    
  x <- table$log2
  y <- x
  y[is.na(x)] <- NA
  y[is.infinite(x)] <- NA
  y[x > max] <- NA
  y[x < min] <- NA
    
  x[is.na(x)] <- min
  x[is.infinite(x)] <- max
  x[x > max] <- max
  x[x < min] <- min
  
  table$log2    <- x
  table$outlier <- y

  table
}

###############################################################################
# narrow_till_normal
#
# input:      x         : raw data
#             weights   : list of weights (one for each x)
#
# output:     list of flags with 1 (out of range) or 0 (in range)
#

narrow_till_normal <- function (x, 
    weights   = NA, 
    min.pval  = 0.05, 
    min.left  = 0,
    low.noise = 1/3,
    test      = 'shapiro'
  ) {
  
  # preset to fail the criterion (want at least 3 values)
  exit <- FALSE
  if (length(x[!is.na(x)]) < 3) {
    exit <- TRUE
  }
#  message(sprintf("Current exit status is %s", exit))
  
  if (max(x, na.rm=TRUE) == min(x, na.rm=TRUE)) {
    exit <- TRUE
  }
#  message(sprintf("Current exit status is %s", exit))
  
  # use standard weights if missing
  if (is.na(weights[1])) {
    weights <- rep(1/length(x), length(x))
  }
  
  # limits
  x.lower <- min(x, na.rm = TRUE)
  x.upper <- max(x, na.rm = TRUE)
  
  # loop until we have a normal subset - or need to exit due to other reasons
  while (!exit) {

    # determine the center of the distribution
    x.center <- weighted.mean(x, weights, na.rm=TRUE)
    x.var <- var.wt(x, weights, x.center, na.rm=TRUE)
    
    if (test == 'shapiro') {
      sw <- shapiro.test(x)
      if (sw$p.value > min.pval) {
        exit <- TRUE
      }
    } else if (test == 'ks') {
      message("ks")
      ks.result <- ks.test(x, "pnorm", mean=x.center, sd=x.sd)
#       print(x)
#       print(ks.result)
#       
      if (ks.result$p.value > min.pval) {
        exit <- TRUE
      }
    } else {
      warning("Test is not supported (none of 'shapiro' or 'ks'), will return input data as is")
      exit <- TRUE
    }
    # check for stop criterion
    if (!exit) {
      low  <- min(x[x >= x.lower], na.rm = TRUE)
      high <- max(x[x <= x.upper], na.rm = TRUE)
      
      ## debug output
      # message(sprintf("center = %.3f (low = %.3f , high = %.3f)", x.center, low, high))
      if ( abs(x.center - low) > abs(x.center - high) ) {
        x[x <= low] <- NA
        x.lower <- low
      } else {
        x[x >= high] <- NA
        x.upper <- high
      } 
      if (length(x[!is.na(x)]) <= length(x) * min.left) {
        exit <- TRUE
      }
      if (length(x[!is.na(x)]) < 3) {
        exit <- TRUE
      }
      if (!exit) {
        
        xsd<-sd(x, na.rm=TRUE)
        if ( is.na(xsd) | xsd < low.noise) {
          exit <- TRUE
        }
      }
    }
  } # end while !exit
  
  outliers <- rep(0, length(x))
  outliers[is.na(x)] <- 1
  return(outliers)
}

###############################################################################
# flag_extremes
#
# input:      x         : raw data
#
# output:     array of values with NA instead of outliers
#

flag_extremes <- function (x, against = 'median', target = 1, min = 10) {
  sd     <- sd(x)
  center <- ifelse(against == 'median', median(x), mean(x))
  lower <- min(x)
  upper <- max(x)
  
  while (!is.na(sd) & sd > target ) {
    # check current low and high values (non-outliers)
    low  <- min(x[x > lower])
    high <- max(x[x < upper])
    if (abs(center - low) > abs(center - high)) {
      # lowest is farther away from center than highest
      lower <- low
#      print(sprintf("lower %.3f", lower))
      
    } else {
      upper <- high
#      print(sprintf("upper %.3f", upper))
    }
    y <- x[x > lower & x < upper]
    center <- ifelse(against == 'median', median(y), mean(y))
    sd     <- sd(y)
  }
  outliers  <- rep(0, length(x))
  outliers[x < lower | x > upper] <- 1
  outliers
}

###############################################################################
# flag_outliers
#
# input:      x         : raw data
#
# output:     array of values with NA instead of outliers
#

flag_outliers <- function(x, na.rm = TRUE, maxnoise = 1, range = 1.5, loop = FALSE, outliers = NA) {
  qnt <- quantile(x, probs=c(.25, 0.5, .75), na.rm = na.rm)
  H <- min( range * IQR(x, na.rm = na.rm)) ## hack to force removal of extreme in noisy regions
  #message(H)
  if (is.na(outliers)[1]) {
    outliers  <- rep(0, length(x))
  }
  seen <- sum(outliers)
  
  outliers[x < (qnt[1] - H)] <- 1         # lower than range x lower quartile
  outliers[x > (qnt[3] + H)] <- 1         # higher than range x upper quartile
  outliers[abs(x - qnt[2])>maxnoise] <- 1 # bigger distance from median than maxnoise
  
  if (loop) {
    found <- sum(outliers)
    
    if (found == 0 | found == seen) return (outliers)
    # found (additional) outliers
    x[outliers == 1] <- NA ## remove values for all flagged outliers
    return(flag_outliers(x, na.rm = TRUE, range = range, loop = TRUE, outliers = outliers))
  } else {
    outliers
  }
  
}

###############################################################################
# subset_gene
#
# input:      table      : full table with amplicon data
#             gene       : name of the gene to subset (column gene)
#
# output:     subset of table
#
subset_gene <- function(table, gene) {
  genedata <- subset(table, gene == genename ) # & !is.na(log2) & is.finite(log2)
  # genedata$outlier <- flag_outliers(genedata$log2, loop = config.recurse.outliers) ## old
  
  if (config.narrow.normal) {
    genedata$outlier <- narrow_till_normal(
        genedata$log2, 
        weights   = genedata$weight, 
        test      = config.narrow.test,
        min.pval  = config.narrow.pval,
        min.left  = config.narrow.min,
        low.noise = config.narrow.lownoise
    )
  }
  else if (config.flag.outliers) {
    genedata$outlier <- flag_outliers(genedata$log2, loop = config.recurse.outliers, maxnoise = config.max.noise)
  }
  else if (config.flag.extremes) {
    genedata$outlier <- flag_extremes(genedata$log2, target = config.max.noise)
  }
  
  # additionally, flag those that fail minmal coverage in the reference
  genedata$outlier[genedata$reference < config.min_cov_for_call] <- 1
  
  genedata
}

###############################################################################
# test_for_change
#
# input:      table       : full table with amplicon data
#             alternative : test against this
#             call        : data frame to insert results
#
# output:     call with changed values
#
test_for_change <- function(
      table, call, 
      mu          =   0, 
      conf.level  =   0.95,
      checkscore  =   5,
      checkcolor  = 'gray50',
      minscore    =  10,
      maxscore    = 100,
      correction  =  NA
  ) {
  
  ## requirement: sufficient amount of data
  if (!call$sufficient) {
    ## defaults are already in place for this, just set report (?) and return
    if (config.report.insufficient) call$report <- TRUE
    return(call)
  }

  ## requirements are met: run the test
  ## //  
  ##     weighted.t.test runs a weighted t.test:
  ##
  ##     test data: column named "log2" from the table
  ##     weights:   column named "weight" from the table
  ## 
  ## No alternatives set, will test "two.sided"
  ##
  wtest <- weighted.t.test( 
    table$log2,
    w          = table$weight,
    mu         = mu,
    conf.level = conf.level
  )
  
  ## insert the results into the call data frame
  call$p.val   <- wtest$p.value
  call$cn.pval <- scientific_10( wtest$p.value )
  call$cn      <- call$refcopy * 2 ^ wtest$estimate
  call$min     <- call$refcopy * 2 ^ wtest$conf.int[1]
  call$max     <- call$refcopy * 2 ^ wtest$conf.int[2]
  call$se      <- wtest$se
  call$sd      <- wtest$sd
  call$qp      <- -10 * log( wtest$se ) # * (call$usable / call$amplicons) 
  call$log2    <- wtest$estimate
  
  ## check better alternatives for this multiplier, e.g. min confident distance
  multiplier <- 1
   if (!is.na(correction)) {
     if (correction == 'log2') {
       multiplier <- abs(wtest$estimate)
     } else if (correction == 'divsq') {
       multiplier <- (wtest$estimate / 0.8)^2
     } else if (correction == 'minconf') {
       multiplier <- min(abs(wtest$conf.int[1]), abs(wtest$conf.int[2]))
     } else if (correction == 'sqrt(sd)') {
       multiplier <- sqrt(call$sd)
     } else if (correction == 'sqrt(sd)*log2') {
       multiplier <- sqrt(call$sd) * (1 + abs( wtest$estimate) )
     } else if (correction == 'sqrt(sd)+log2') {
       multiplier <- sqrt(call$sd) + (1 + abs( wtest$estimate) )
     } else if (correction == 'sigma') {
       multiplier <- sqrt(call$se)
     } else if (correction == 'phi') {
       multiplier <- sqrt(call$max - call$min)
     } else if (correction == 'bonferroni') {
       multiplier <- 1 / call$usable
     }
  }
  # apply scoring (and correction if multiply is not == 1)
  call$score   <- min(c(maxscore, -10 * log10(wtest$p.value) * multiplier)) #  * abs(wtest$estimate)
  
  if (is.nan(call$score)) {
    call$score  <- 0
    call$usable <- 0
  }
  
  if (call$usable < config.minimal_observations) {
    call$filter <-  paste('o', config.minimal_observations, sep="")
  }
  ## score *must* be relevant
  else if ( call$score >= config.score.report ) {
      call$cn.label <- sprintf("%.1f ~ %.1f", call$min, call$max) # sprintf("%.1f (+-%.1f)", call$cn, call$se)
#      call$cn.label <- paste(ifelse(call$cn > call$refcopy, "> ", "< "), call$refcopy)
      call$cn.color <- calc_color(call$cn, call$refcopy)
      call$report   <- TRUE
      call$filter   <- 'PASS'
      #call$alt      <- ifelse( call$cn < call$refcopy , '<DEL>', '<AMP>')
    } else if ( call$score >= checkscore   & 
      ( call$min   > call$refcopy | call$max < call$refcopy ) ) {
      call$cn.label <- sprintf("%.1f ~ %.1f", call$min, call$max)
      call$cn.color <- checkcolor
      call$filter   <- paste('q', config.score.report, sep="")
      #call$alt      <- 'CHECK'
      call$report   <- TRUE
    }
  return(call)
}

###############################################################################
# set_filter
#
# Test if provided score is higher than pass or check and return a tag
# to put into the filter column of vcf
# 
# input:      score     : the score to check
#             pass      : threshold to set PASS
#             check     : threshold to set q<pass> (check < score < pass)
#             fail      : what to return for score < check (NA = q<check>)
#
# output:     filter    : example: PASS; q40, q20, FAIL
#
set_filter <- function (
    score, 
    obs,
    pass  = config.score.report, 
    check = config.score.check, 
    mino  = config.minimal_observations,
    fail  = NA
    ) {
  
    if (obs < mino ) {
      return(paste('o', mino, sep=""))
    }
    if (score >= pass) {
      return('PASS')
    } else if (score >= check) {
      return(paste('q', pass, sep=""))
    } else {
      if (is.na(fail)) {
        return(paste('q', check, sep=""))
      }
      else {
        return(fail)
      }
    }
}

###############################################################################
# test_for_segmentation
#
# input:      table       : full table with amplicon data
#             alternative : test against this
#             call        : data frame to insert results
#
# output:     call with changed values
#
test_for_segmentation <- function(
      table, call, 
      mu = 0, 
      cl = config.conf_level, 
      minscore = config.score.report,
      maxscore = config.score.max,
      segdiffp = config.segment.diff.pval
  ) {

  ## first requirement: sufficient amount of data
  if (!call$sufficient) {
    ## defaults are already in place for this, just set report (?) and return
    if (config.report.insufficient) call$report <- TRUE
    return(list(call=call, segments=NA))
  }
  if (!config.report.segmented) {
    call$sg.color <- config.color.untested
    call$sg.label <- "segmentation not tested" ## first line of right panel
    call$sg.pval  <- NA                        ## second line of right panel
    return(list(call=call, segments=NA))
  }

  ## run the segmentation routine on the data
  segments<-segment(table, refcopy = call$refcopy, maxscore = maxscore, cutoff = segdiffp)
  
  ## check number of segments
  if (nrow(segments) > 1) {

    call$segmented <- TRUE
#    call$report    <- TRUE ## this should be checkin in main script
    
    cn.list       <- paste(sprintf("%.0f", segments$cn), collapse=" | ")
    call$sg.label <- paste("segmented (cn = ", cn.list, ")", collapse="")
    call$sg.color <- config.color.segmented
#    for (segnum in 1:nrow(segments)) {
#      # assemble results for final output
#      current <- data.frame(
#          gene            = call$genename, 
#          chr             = call$chr, 
#          from            = segments$left[segnum], 
#          to              = segments$right[segnum], 
#          cn              = segments$cn[segnum], 
#          min             = segments$min[segnum], 
#          max             = segments$max[segnum], 
#          se              = segments$se[segnum], 
#          total_amplicons = segments$total_amplicons[segnum]
#      )
#    }
    return(list(call=call, segments=segments))
    #			paste(sprintf("%.1f", c(1.202,1.302,1.402)), collapse=":")
  }
  else {
    call$sg.color <- config.color.even
    call$sg.label <- "biased, but not segmented" ## first line of right panel
    #call$sg.pval  <- NA                   ## second line of right panel
    return(list(call=call, segments=NA))
  }
}

###############################################################################
# calc_pvalue
#
# Calculate the probability, that the ratio is derived from a 
# normal distrubution around *mu* with sd based on the provided model.
# The model should contain coeficients of a quantile regression.
#
# input:      model     : model params (array with 6 values)
#             ratio     : current log2 ratio
#             weight    : weight factor, e.g. coverage, same as for model
#   [opt]   ( mu        : center of the distribution to probe             )
#
# output:     probability
#
calc_pvalue <- function(ratio, weight, sd = 1, mu = 0) {
  # if not center is given, use the original model (typically around 0)
	
  return(dnorm(ratio, mean=mu, sd=model_sd)*weight)
}

###############################################################################
# calc_model_pvalue
#
# Calculate the probability, that the ratio is derived from a 
# normal distrubution around *mu* with sd based on the provided model.
# The model should contain coeficients of a quantile regression.
#
# input:      model     : model params (array with 6 values)
#             ratio     : current log2 ratio
#             weight    : weight factor, e.g. coverage, same as for model
#   [opt]   ( mu        : center of the distribution to probe             )
#
# output:     probability
#
calc_model_pvalue <- function(model, ratio, weight, mu = NA){
  # if not center is given, use the original model (typically around 0)
  if (is.na(mu)) {
    mu <- (weight * model[4]) + model[3]
  }
  model_sd <- max(0.05, ( (weight * model[6]) + model[5] ) - ( (weight * model[2]) + model[1]))
  return(dnorm(ratio, mean=mu, sd=model_sd))
}

###############################################################################
# determine.sex
#
# Return the word 'female' or 'male' based on input that can also 
# consist of the chromosome representations 'XX' or 'XY'.
#
# input:      sex       : input string
# 
# output:     either 'female' [default] or 'male'
#
determine.sex <- function(sex) {
  # default to return 'female'
  if (is.na(sex) | sex == "") {
    return('female')
  }  
  if (length(grep("y", sex, ignore.case = TRUE)) | 
      length(grep("^male", sex, ignore.case = TRUE, perl = TRUE))
  ) {
    return('male')
  }
  return('female')
}

###############################################################################
# calc_phred
#
# Simply convert an error probability to a phred score.
# 
# input:      pvalue    : error probability
#   [opt]   ( digits    : round to this accurary [0] )
#
# output:     phred
#
calc_phred <- function(pvalue, digits = 0){
	phred <- round(-10 * log10(pvalue), digits)
	phred <- max(1,min(255,phred))
	return(phred)
}

###############################################################################
# calc_color
#
# This function will provide a color for a combination of seen and expected
# copy numbers. Unchanges (similar) values will be green, and lower seen
# values will red, higher values blue to purple.
#
# input:      copy      : number of copies found
#             refcopy   : number of copies in the reference
#
# output:     color
#
mygrad <- gradient_n_pal(
    colours=c("brown", "red", "green","blue", "purple", "black"), 
    values=c(-5,-1,0,log2(3/2),1,5)
)

calc_color <- function(copy, refcopy = 2) {
  v <- mygrad(log2(copy/refcopy))
  return(v)
}

#calc_color <- function(copy, refcopy){
#  ratio<-round(log2(copy/refcopy)*20)+15
#  ratio<-min(40, ratio)
#  ratio<-max(1, ratio)
#  return(rainbow(40)[ratio])
#  #	return('#CC999999')
#}

###############################################################################
# scientific_10
#
# input:      number              # example: 0.00021
# output:     converted format    # example: 2.1 x 10^-4
#
scientific_10 <- function(x) {
	return(parse(text=gsub("e", " %*% 10^", paste(sprintf("%3.2g", x)))))
}

###############################################################################
# scientific_pval
#
# input:      number              # example: 0.00021
# output:     converted format    # example: p=2.1 x 10^-4
#
scientific_pval <- function(x) {
	return(parse(text=gsub("e", " %*% 10^", paste(sprintf("p=%3.2g", x)))))
}
###############################################################################
# calc_impure_log2
#
# Function to calculate the expected log2 for a given combination of inputs:
#
# input:      sample    : number of copies in sample
#             reference : number of copies in reference
#             purity    : purity of sample (between 0 to 1)
#
calc_impure_log2 <- function(sample, reference = 2, purity = 1) {
  # some sanity checks
  sample    <- min(9, max(0, sample))    # force sample into allowed range
  reference <- min(9, max(1, reference)) # force reference into allowed range
  purity    <- min(1, max(0, purity))    # force purity into allowed range
  
  if (sample == 0) {
    return(-Inf)
  }
  admix <- (sample * purity) + (reference * (1 - purity) )
  return( log2( admix / reference ) )
}

###############################################################################
# vcf_header
#
# input:      assembly          : assmebly name (id), examples: hg19, GRCh38
#             filename          : name (path) of the file to read data from
#             date              : date to put in VCF in YYYYMMDD (20140331)
#             version           : version of VCF to put in header (4.1)
#
# Function to calculate the expected log2 for a given combination of inputs:
#
vcf_header <- function(
    assembly = config.assembly.version,
    filename = config.assembly.filename,
    date     = format(Sys.Date(), format="%Y%m%d"),
    version  = '4.1'
) {
  
  current <- read_chromosomes(assembly, filename)
  
  lines   <- c(
      paste("##fileformat=VCFv", version, sep=""), 
      paste("##fileDate=", date, sep=""),
      paste("##source=", config.software.name, 'V', config.software.version, sep=""),
      paste("##reference=", assembly, ".fas", sep="")  # this is a fake name!  
  )
  
  shift.by <- length(lines)
  for (row in 1:nrow(current)) {
    lines[row + shift.by] <- sprintf("##contig=<ID=%s,assembly=%s,length=%s>", 
        current$chromosome[row], 
        current$assembly[row], 
        current$total_length[row]
    )
  }  
  return(paste(lines, sep="", collapse="\n"))
}

###############################################################################
# read_chromsomes
#
# Function to read chromosome data from a file
#
# input:      assembly          : assmebly name (id), examples: hg19, GRCh38
#             filename          : name (path) of the file to read data from
#
#             * the input file can be either a custom tabular format 
#               (there is an example shipped with this code) or a .fai
#               index file generated by samtools. The type is choosen based 
#               on the presence of '.tsv' or '.fai' in the filename!
#
# output:     Data table with headers defined by the first line in file.
#
read_chromosomes <- function(
    assembly = config.assembly.version, 
    filename = config.assembly.filename
) {
  
  if(length(grep(".tsv", filename))) {
  assemblies <- read.table(filename, header=TRUE)
  
  if (length(grep("^hg", assembly))) {
    names(assemblies)[4] <- 'assembly'
  } else {
    names(assemblies)[5] <- 'assembly'
  }
  
  assemblies <- assemblies[assemblies$assembly == assembly, ]
  assemblies <- assemblies
  return(assemblies)
  }  else if(length(grep(".fai", filename))) {
    assemblies <- read.table(
          filename, 
          header=FALSE, 
          sep="\t", 
          col.names = c("chromosome", "total_length", "offset", "line", "bytes")
    )
    assemblies$assembly <- rep(assembly, nrow(assemblies))
    assemblies <- assemblies
    return(assemblies)
  }
  # croak if we get here!
  warning("No assembly table generated - input data missing or not recognized!")
}


###############################################################################
# plot_overview
#
# Function to plot calls from all regions sorted by location
#
# input:      data             : the table that is also exported to csv
#             filename         : name (path) of the file to wrtie to
#
# output:     This will write the PDF plot to the filename
#


plot_overview <- function(csv, filename=NA, title=NA) {
  
  data <- read.csv(csv)
  data <- subset(data, !is.na(p.val)) # subset to keep only valid regions
  # replace 0 scores with check score for plotting (gray, not green)
  data$score[is.na(data$min)] <- config.score.check 
  
  yrange <- range(data$min, data$max, 0, 4, na.rm=TRUE)
  
  if (nrow(data)> 0) {
    
    # collect brakes, labels and lines for chromosome axis
    breaks<-c()
    labels<-c()
    lines<-c(0) # start with left line
    
    last.pos <- 0
    temp <- data.frame()
    last.row <- 0
    
    # get all chromosomes in natural order
    chromosomes <- naturalsort(unique(data$chromosome), decreasing=FALSE)
    
    # loop through all chromosomes
    for(chr in chromosomes) {    
      sub <- subset(data, chromosome==chr) # subset to keep only current chr
      
      last.row <- nrow(sub)
      if(last.row==0) next
      
      for (target in sub$locus[order(sub$start)]) {
        region <- subset(sub, locus==target)
        region <- region
        last.pos <- last.pos +1 # move on axis, first wil be 1
        region$position <- last.pos
        temp<-rbind(temp, region)
      }
      # add one more for the ending line on the right
      breaks <- c(breaks, last.pos - ((last.row-1)/2))
      last.pos <- last.pos +1
      labels <- c(labels, chr)
      lines <- c(lines, last.pos)
    }
    
    p <- ggplot(temp, aes(x=position, y=copies, colour=score, ymin=min, ymax=max)) #+ theme_overview()
    p <- p + geom_vline(xintercept = lines, color="white", size=0.5)
    p <- p + geom_pointrange(size=0.5)
    p <- p + scale_y_continuous('copies', limits=yrange)
    p <- p + scale_x_continuous('chromosome', breaks=breaks, labels=labels, expand=c(0.01,0.01), limits=c(0,last.pos))
    p <- p + scale_colour_gradientn(
      "score", 
      values = c(
        0,                                           # min
        0.5 * config.score.check / config.score.max, # 0.5 x check score
        config.score.check       / config.score.max, # 1.0 x check score
        config.score.report      / config.score.max, # 1.0 x call score
        config.score.report * 2  / config.score.max, # 2.0 x call score
        1                                            # max
      ), 
      colours = c(
        "green",  # min
        "green",  # 0.5 x check
        "gray",   # 1.0 x check
        "orange",    # 1.0 x call
        "red",   # max
        "brown"
      ),
      limits = c(0,config.score.max),
      guide  = guide_colorbar(label.theme = element_text(size = 7, angle = 90))
    )
    
    if(!is.na(title)) p <- p + labs(title = title)
    
  #  p <- p + theme(
  #      axis.text.x = element_text(angle=90, hjust = 0.5, vjust = 0.5),
  #      panel.grid.major.x = element_blank(), 
  #      panel.grid.minor.x = element_line(colour="white", size=0.5), 
  #      panel.grid.major.y = element_line(colour="grey50", size=0.5, linetype="dotted"),
  #      panel.grid.minor.y = element_blank(), 
  #      legend.position = "bottom"
  #  )
    if (!is.na(filename)) {
      pdf(filename, paper="US", colormodel = config.report.pdf.colormodel)
      print(p)
      dev.off()
    }
    else {
      print(p)
    }
  }
  else {
    if(!is.na(filename)) {
    pdf(filename, paper="US", colormodel = config.report.pdf.colormodel )
    plot(1, type="n", axes=F, xlab="", ylab="")
    mtext("no data for overview plot")
    dev.off()    
    }
    else {
      print("no data for overview plot")
    }
  }
  #
  ## finished
}
