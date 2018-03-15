# enable debugging output
options(error = quote({
  sink(file="R_traceback.txt");
  dump.frames();
  print(attr(last.dump,"error.message"));
  traceback();
  sink();
  q()}))

# save current directory, change to code directory (need to remove hard-code here)
dirsave <- getwd()
setwd("/srv/qgen/code/qiaseq-dna/copy_number/R/")

# import libraries and own sources
#suppressPackageStartupMessages(
  source("libraries.R")
#)

# set directory back to start directory
setwd(dirsave)

##' overriding general settings to collect performance statistics
##' source("config_woowt.R") # TODO: COMMENT OR REMOVE THIS LINE !!!
##'

# process command line - provide coverage of reference + sample
args<-commandArgs()

# read the reference and sample coverage files,
# and provide the expected counts for X and Y
depth <- import_counts(
    sample      = args[4],
    sample.x    = as.numeric(args[5]),
    sample.y    = as.numeric(args[6]),
    reference   = args[7], 
    reference.x = as.numeric(args[8]),
    reference.y = as.numeric(args[9])
)

## get the list of chromosomes with data from the table
chromosomes <- unique(naturalsort(depth$chr))


## // START OUTPUT - This is for visual tracking of what is being done
##
## path and filenames
output.dir       <- args[10]
output.basename  <- args[11]
if (!is.na(args[12]) & args[12] != "") {
  #message("setting genome file to ", args[12])  
  config.assembly.filename   <- args[12] ## defaults to: './chromosome_lengths.tsv'
  
}
if (!is.na(args[13]) & args[13] != "") {
  #message("setting assembly version to ", args[13])
  config.assembly.version   <- args[13] ## defaults to: 'hg19'
}

options(error = quote(
    {
      sink(file=paste(output.dir , '/', output.basename, '.err', sep=""));
      dump.frames();
      print(attr(last.dump,"error.message"));
      traceback();
      sink();
      q()
    }
  )
)

## chromosome handling
## 
## read the total length (used for fully scaled plots only)
## calling read_chromosomes will use defaults for unset params:
## -> filename = config.assemblies.filename
## -> assembly = config.assemblies.assembly
chrlen <- read_chromosomes() 

## pdf
pdf( paste(output.dir , '/', output.basename, '.pdf', sep=""),
     paper       = config.report.pdf.paper,
     colormodel  = config.report.pdf.colormodel,
     title = paste(
         "QIAGEN ",
         config.software.name, ' v', config.software.version,
         sep=""
     )
)

###############################################################################
##
## plot the data as is - if requested in config
##
if (config.plot.rawdata) {
  p <- ggplot(depth, aes(reference,sample))
  p <- p + geom_point(aes(color=chr), alpha=1/3, size=1, na.rm=TRUE) + geom_abline()
  p <- p + scale_x_log10("raw reference") + scale_y_log10("raw sample") 
  p <- p + labs(title="Raw counts scatterplot")
  p
}

###############################################################################
##
##  normalize counts
##
## This will create the log2 from the ratios and shift by the median 
## to normalize the whole set (median =defined=> 0).
# calculate the raw ratios
depth$log2  <- log2(depth$sample / depth$reference)
depth$mincov  <- (depth$sample + depth$reference)

# filter for minimal depth in either ref or sample
usable <- subset(depth, 
    !is.na(log2) & is.finite(log2) & mincov >= config.min_cov_for_scale
)
depth$log2  <- depth$log2 - median(usable$log2, na.rm = TRUE)

# the next step will replace log2 values < min and > max by min and max,
depth <- pull_in_range( 
      depth, 
      min    = config.log2.min,         # < min = min
      max    = config.log2.max
)

## store the weights that will be used later to perform:
## 
## initialize weights to be all 1
depth$weight <- rep(1, nrow(depth))

## use real weights (if requested)
if (!is.na(config.use.weights)) {
    
  if  (config.use.weights == 'squared') {
    # this will divide the greater of the two normalized values by the mininmal
    # required count, then square the result (bias towards use of mainly big values)
    depth$weight  <- apply(depth, 1, function(x) as.numeric(max(x[13], x[14])))
    depth$weight <- (depth$weight / config.min_cov_for_call)^2
  } else if (config.use.weights == 'raw_max') {
    # the greater of the two raw counts
    depth$weight  <- apply(depth, 1, function(x) as.numeric(max(x[6], x[7]))) # raw
  } else if (config.use.weights == 'norm_max') {
    # the greater of the two normalized counts
    depth$weight  <- apply(depth, 1, function(x) as.numeric(max(x[13], x[14]))) # normalized
  } else if (config.use.weights == 'log2_norm_max') {
    # the greater of the two normalized counts, but log2 converted
    depth$weight  <- apply(depth, 1, function(x) as.numeric(max(x[13], x[14]))) # normalized
    depth$weight <- log2(depth$weight)
  } else if (config.use.weights == 'sum(raw)') {
    # raw counts summed
    depth$weight <- (depth$sample + depth$reference)     # 
  } else if (config.use.weights == 'log2(sum(raw))') {
    # log2 of raw counts summed
    depth$weight <- log2(depth$sample + depth$reference) # 
  } else if (config.use.weights == 'norm_sum') {
    # normalized counts summed
    depth$weight <- (depth$sample_norm + depth$reference_norm) / 2 ## check best solution for this
  }

}

#if (config.plot.normalized) {
#  xlimits <- quantile(c(depth$rawlog2, depth$log), c(0.01, 0.98), na.rm = TRUE)
#  mv      <- max(abs(xlimits))
#  xlimits <- c(- mv, mv) # symmetric limits
#  
#  # plot the rawlog2/log2 comparison
#  cd <- ggplot(depth)
#  cd <- cd + geom_density( aes( x=log2(sample/reference) ), adjust=1/2, fill='red',   alpha=1/3,  na.rm = TRUE)
#  cd <- cd + geom_density( aes( x=log2 ),    adjust=1/2, fill='green', alpha=1/3,  na.rm = TRUE)
#  cd <- cd + labs(title="Raw and normalized log2-ratios")
#  cd <- cd + scale_x_continuous("log2(sample/reference)", limits=xlimits)
#  cd
#}
#
## filter for minimal depth in either ref or sample after scaling
#depth   <- subset(depth, sample >= config.min_cov_for_call | reference >= config.min_cov_for_call )
#
## calculate the copy number per row
depth$cn <- depth$refcopy * (2 ^ depth$log2)
maxcov   <- max(depth$sample, depth$reference, na.rm = TRUE)
#
#
#if (config.plot.normalized | 1 == 1) {
#  # plot the usable, normalized data after scaling
#  normplot <- ggplot(depth, aes(reference_norm,sample_norm))
#  normplot <- normplot + geom_point(aes(colour=chr), alpha=1/2, size=1, na.rm=TRUE) + geom_abline()
#  normplot <- normplot + scale_x_log10("normalized reference") + scale_y_log10("normalized sample")
#  normplot <- normplot + labs(title="Scaled datasets scatterplot")
#  normplot
#}

# calculate the absolute max

# remove undefined and infinite values for creating the model
#defined_depth <- subset(depth, !is.na(log2) & is.finite(log2))
#defined_depth <- subset(defined_depth, mincov >= config.min_cov_for_call )

## calculate data to plot the model
#xlimits     <- c(config.min_cov_for_call, min(config.min_cov_for_call * 100, maxcov))
#coverage    <- seq(config.min_cov_for_call, maxcov, length=nrow(defined_depth))
#coverages   <- c(1,25,50,100) * config.min_cov_for_call

# plot chromsome boxplots
#if (config.plot.chrom_boxes) {
#  # try to plot all chromosome boxplots at once  
#  allchr <- ggplot(depth)
#  allchr <- allchr + geom_jitter(aes(y=cn, x=factor(chr), size=mincov, colour=log2), alpha=0.25, na.rm=TRUE)
#  allchr <- allchr + geom_boxplot(aes(x=factor(chr), y=cn), outlier.size=0, fill = "#ffffff00", na.rm=TRUE)
#  allchr <- allchr + scale_colour_gradientn(colours=c("black", "red", "green", "blue", "purple"), breaks=c(-2,-1,0,1,2), limits=c(-2,2))
#  allchr <- allchr + scale_size("coverage", trans="log", range=c(1, 3), limits=c(config.min_cov_for_call,maxcov), breaks=c(config.min_cov_for_call,config.min_cov_for_call*10,config.min_cov_for_call*100))
#  allchr <- allchr + scale_y_continuous(limits=config.plots.ylimits, breaks=config.plots.ybreaks, labels=config.plots.ylabels)  
#  allchr <- allchr + scale_x_discrete(limits=chromosomes)
#  allchr <- allchr + theme(legend.position="bottom")
#  allchr <- allchr + labs(list(x="copies"))
#  allchr
#}


# plot every chromosome individually
#if (config.plot.chromosomes) {
#  
#  grid.newpage()
#  
#  vprow     <- 1
#  vprow_max <- 4
#  maxchrlen <- max(chrlen$total_length)
#  vph       <- (1/vprow_max)
#  vpx       <- 0
#  vpy       <- 1+(vph*0.5)
#  
#  for(chromosome in chromosomes) { #c(1,2,6,22,'X')) {
#    
#    # print(sprintf("Working on chromosome %s ...", chromosome))
#    
#    chrdata <- data.frame(subset(depth, chr == chromosome))
#    chrdata <- chrdata
#    refcopy <- unique(chrdata$refcopy)
#    thislen <- chrlen[ which(chrlen$chromosome == chromosome), ]
#    xlimits <- c(0, thislen$total_length / 1000000)
#    
#    # precalc the color settings for this (needs refcopy)
#    color.scale <- list()
#    for (p in c(0:config.plots.max.cn)) {
#      color.scale[p+1] <- calc_color(p, refcopy)
#    }
#    # start plotting ...
#    chromosome.plot <- ggplot(chrdata, aes(x=position / 1000000, y=cn, colour=log2)) + theme_chromosome()
#    
#    # add layers
#    chromosome.plot <- chromosome.plot + geom_point(alpha=0.05, na.rm=TRUE) 
#    chromosome.plot <- chromosome.plot + scale_size("coverage", trans="log", range=c(0.5, 2), limits=c(config.min_cov_for_call,maxcov), breaks=coverages)
#    chromosome.plot <- chromosome.plot + scale_colour_gradientn(colours=c("black", "red", "green", "blue", "purple"), breaks=c(-2,-1,0,1,2), limits=c(-2,2))
#    chromosome.plot <- chromosome.plot + xlim(xlimits)	
#    chromosome.plot <- chromosome.plot + scale_y_continuous(limits=c(0,config.plots.max.cn))
#    chromosome.plot <- chromosome.plot + labs(list(x = sprintf("chr%s [Mbp]",chromosome ), y = "copies"))
#    
#    vpw   <- 0.125 + ( (thislen$total_length / maxchrlen) * 0.875)
#    vpy   <- vpy - vph
#    
#    # start a new page if the last is full
#    if(vprow > vprow_max) {
#      vprow <- 1
#      vpy   <- 1 + (vph * 0.5) -vph
#      grid.newpage()
#    }
#    # print to the viewport
#    print(chromosome.plot, vp=viewport(x=vpx, y=vpy, width=vpw, height=vph, just=0))
#    vprow <- vprow + 1
#  }
#}
#
plots   <- list()
counter <- 0
## // create lists to store results
##    for exporting after the loop

## regions
genes.all <- data.frame()

# vcf
vcf <- data.frame()

## segments
segments.all <- data.frame()

#print(unique(depth$gene))

# sort the vcf table by chromosome / position
depth <- depth[naturalorder(paste(depth$chr, depth$position)),] 
allgenes <- depth$gene
lastgene <- '_NOT_A_GENE_'

## loop all regions
for (genename in unique(depth$gene)) {
#  if (genename == lastgene) next
#  lastgene <- genename
  message(genename)
#  next;
  # subset gene data
  genedata <- subset_gene(depth, gene = genename)
  n.total  <- nrow(genedata)
  
  # smoothed data (no outliers)
  majority  <- subset(genedata, outlier < 1)
  majority  <- majority
  n.usable  <- nrow(majority)
  #majority.dip <- dip.test(sort(majority$log2))  
  
  # use this data frame to centrally store values to report
  call <- data.frame(
      
      ## gene specific data, tunneled through for reporting
      genename    = genename,
      chr         = genedata$chr[1], 
      from        = min(genedata$position),
      to          = max(genedata$position),
      refcopy     = genedata$refcopy[1],
      
      ## summary of data, to make decisions based on these
      amplicons   = n.total,
      outliers    = n.total - n.usable,
      usable      = n.usable,
      sufficient  = ifelse( n.usable > config.minimal_observations, TRUE, FALSE),
      #bimodal     = ifelse( majority.dip$p.value < config.segment.p.value, TRUE, FALSE),

      ## in case of copy number changes, collect statistics for that
      alternative = NA,
      p.val       = NA,
      conf.level  = config.conf_level,
      score       = 0,
      qp          = 0,
      filter      = set_filter(0, n.usable, config.score.report, config.score.check, config.minimal_observations),
  ## paste('q', config.score.check, sep=""), # 'UNCH',

      ## provide an alternative copy number call from weighted.t.test
      cn          = genedata$refcopy[1],
      min         = NA,
      max         = NA,
      se          = NA,
      log2        = NA,
      sd          = NA,
      
      ## test for segmentation of the gene region
      segmented   = FALSE,
      seq.p.value = NA,#1/majority.sd,
      segments    = NA,
      
      ## collect strings and colors for cn output, defaults to "insufficient"
      cn.color    = config.color.unchanged,
      cn.label    = "unchanged",                ## first  line of left panel
      cn.pval     = 1,                         ## second line of left panel
      
      ## collect strings and colors for segmentation, defaults to "even"
      sg.color    = config.color.even,
      sg.label    = "unimodal",           ## first  line of right panel
      sg.pval     = 1, #majority.dip$p.value,  ## second line of right panel
      
      ## flag to set for PDF reporting
      report      = config.report.all
  )
  
  # check for general copy number change (different from reference)
  call <- test_for_change(
        majority, call, 
        conf.level = config.conf_level, 
        minscore   = config.score.report,
        maxscore   = config.score.max,
        correction = config.score.correction,
        checkscore = config.score.check,
        checkcolor = config.color.check
    )  
    
  segments <- NA
  
  if (!call$sufficient) {
    call$cn.color <- config.color.insufficient
    call$cn.label <- "insufficient"
    call$cn.pval  <- ''
    call$filter   <- paste('o', config.minimal_observations, sep="")
    call$score    <- 0
  }
  
  if (config.report.segmented) {
    # check for segmentation inside this region (if requested)
    segresults <- test_for_segmentation(
          genedata, call, 
          config.score.report, 
          config.score.max,
          segdiffp = config.segment.diff.pval
    )
    call       <- segresults$call        ## resolve the returned objects
    segments   <- segresults$segments     ## from the named list
    if ( call$segmented ) {
      call$report <- TRUE
    }
  }
  
  xlabel<-sprintf("%s ~ %.1f Mbp", call$chr, mean(genedata$position)/1000000)
  
  if (call$report & call$usable > 0) {
    
    ## // the left panel showing the violin plot for all amplicons of the gene
    ## * new from 0.8.7: the violin now also uses weights
    ##
    panel.left <- ggplot(majority, aes(factor(gene), cn, weight=weight/sum(weight))) # + theme_gene_left() # , label=call$cn.label
    
    ## add layers
    panel.left <- panel.left + geom_abline(intercept=call$refcopy, slope=0, colour="green", size=4, alpha=0.25)
    panel.left <- panel.left + geom_violin(fill=call$cn.color, adjust=config.plots.smooth.violin, alpha=config.plots.alpha.violin, adjust=1.5)
    panel.left <- panel.left + scale_y_continuous(limits=config.plots.ylimits, breaks=config.plots.ybreaks, labels=config.plots.ylabels)  
    panel.left <- panel.left + labs(list(y = "copies", x = sprintf("n = %.0f (%.0f)", call$usable, call$amplicons)), title=xlabel)
    panel.left <- panel.left + geom_text(x=1, y=config.plots.line1, label=paste(call$cn.label), vjust=0, size=3.5, color=call$cn.color, parse=FALSE)
    panel.left <- panel.left + geom_text(x=1, y=config.plots.line2, label=paste(round(call$score,0), " / ", round(call$qp, 0)), vjust=0, size=3.5, color=call$cn.color, parse=FALSE)
    
    ## // the right panel showing the amplicons in chromosomal context (zoom)
    ##
    panel.right <- ggplot(genedata, aes(position / 1000000, cn)) # + theme_gene_right()
    
    ## add layers
    panel.right <- panel.right + geom_abline(intercept=call$refcopy, slope=0, colour="green", size=4, alpha=0.25)
    panel.right <- panel.right + scale_size("coverage", trans="log", range=c(1, 3), limits=c(config.min_cov_for_call,maxcov), breaks=c(config.min_cov_for_call,config.min_cov_for_call*10,config.min_cov_for_call*100))
    
    # no outliers for segmented data, regardless of outlier settings
    if (call$segmented) {

      # no outliers here
      panel.right <- panel.right + geom_point(aes(size=mincov), alpha=0.25, na.rm=TRUE) 
      
      # write the segemantation text into the panel
      panel.right <- panel.right + geom_text(x = min(genedata$position)/1000000, y=config.plots.line1, label=paste(call$sg.label), vjust=0, hjust=0, size=3.5, color=call$sg.color, parse=FALSE)
      
      # draw vertical lines for segmented regions
      for (region in 2:nrow(segments)) {
        ##        if (!segments$is.changed[region]) next  #skip unchanged
        breakpoint  <- segments$left[region]/1000000
        panel.right <- panel.right + geom_vline(
            xintercept = breakpoint, 
            color      = config.plots.seg.lcolor, 
            linetype   = config.plots.seg.ltype
        )
      }
    }
    else {
      # not segmented, so outliers may be removed if requested
      if(config.flag.outliers | config.narrow.normal | config.flag.extremes) {
        # plot outliers and indicate by color change
        panel.right <- panel.right + geom_point(aes(size=mincov, colour=factor(outlier)), alpha=0.25, na.rm=TRUE) 
        panel.right <- panel.right + scale_colour_manual("outlier", values = c("black", "red"), labels=c("no", "yes"))
      } else {
        # plot all values
        panel.right <- panel.right + geom_point(aes(size=mincov), alpha=0.25, na.rm=TRUE) 
      }
    }
    # label x axis and scale y
    panel.right <- panel.right + labs(list(x = sprintf("chr%s [Mbp]", call$chr)))
    panel.right <- panel.right + scale_y_continuous(limits=config.plots.ylimits, breaks=config.plots.ybreaks, labels=config.plots.ylabels) 
    
    
    if (config.report.genes) {
      # save the individual gene to a single file
      multiplot(
          plotlist = list(panel.left, panel.right), 
          file     = paste(output.dir , '/', output.basename, '_', genename, '.pdf', sep=""),
  #        width    = config.png.width,
  #        height   = config.png.height,
          paper   = config.gene.pdf.paper,
          layout  = config.gene.layout
      )
    }
  } # end if call$report and usable > 0
  else {
    ## this is purely for regions without any data
    panel.left <- ggplot(genedata, aes(factor(gene), cn, weight=weight/sum(weight))) #+ theme_gene_left() # , label=call$cn.label
    
    ## add layers
    panel.left <- panel.left + geom_abline(intercept=call$refcopy, slope=0, colour="green", size=4, alpha=0.25)
    panel.left <- panel.left + geom_blank()
    panel.left <- panel.left + scale_y_continuous(limits=config.plots.ylimits, breaks=config.plots.ybreaks, labels=config.plots.ylabels)  
    panel.left <- panel.left + labs(list(y = "copies", x = sprintf("n = %.0f (%.0f)", call$usable, call$amplicons)), title=xlabel)
    panel.left <- panel.left + geom_text(x=1, y=config.plots.line1, label=paste("no data"), vjust=0, size=3.5, color="gray50", parse=FALSE)
    panel.left <- panel.left + geom_text(x=1, y=config.plots.line2, label="- / -", vjust=0, size=3.5, color="gray50", parse=FALSE)
    
    panel.right <- ggplot(genedata, aes(position / 1000000, cn)) #+ theme_gene_right() + geom_blank()
    
    ## add layers
    panel.right <- panel.right + geom_abline(intercept=call$refcopy, slope=0, colour="green", size=4, alpha=0.25)
    panel.right <- panel.right + scale_size("coverage", trans="log", range=c(1, 3), limits=c(config.min_cov_for_call,maxcov), breaks=c(config.min_cov_for_call,config.min_cov_for_call*10,config.min_cov_for_call*100))
    panel.right <- panel.right + labs(list(x = sprintf("chr%s [Mbp]", call$chr)))
    panel.right <- panel.right + scale_y_continuous(limits=config.plots.ylimits, breaks=config.plots.ybreaks, labels=config.plots.ylabels) 
    
  }

  # collect this plot	
  counter <- counter + 1
  plots[[counter]] <- panel.left
  counter <- counter + 1
  plots[[counter]] <- panel.right
  
  # do the multiplotting
  if (counter == 4) {
    multiplot(plotlist = plots, layout=config.report.layout)
    plots<-list()
    counter = 0
  }
  
  # assemble results for final output
  # add to lists
#  if (call$segmented) print(call)
  call$cn.pval <- NULL
  
  if (!call$segmented) {
    genes.all    <- rbind(genes.all, call)
    #print(call)
  }
#  else {
#    print(head(call))
#    q()
#  }
  # dispatch results to lists for export
  if (call$report) {
    if (!call$segmented) {
      bp      <- call$to -  call$from + 1
      refline <- subset(depth, chr == call$chr & position == call$from  )

    if (call$filter == 'PASS') {
      vcf <- rbind(
          vcf, data.frame(
              CHROM   = call$chr, #paste('chr', call$chr, sep=""),
              POS    = call$from, 
              ID      = '.',
              REF     = refline$base[1],
              ALT      = ifelse( call$cn < call$refcopy , '<DEL>', '<DUP>'),
              QUAL     = round(call$score),
              FILTER   = call$filter,
              INFO     = paste(
                  'IMPRECISE', ';',
                  'SVTYPE=', ifelse( call$cn < call$refcopy , 'DEL', 'DUP'), ';',
                  'END=', call$to, ';',
                  'SVLEN=', ifelse( call$cn < call$refcopy , -bp, bp), ';',
                  'LOCUS=', call$genename,
                  sep=""
              ),
              FORMAT    = 'CN:CNMIN:CNMAX:CNQ:PVAL',
              SAMPLE    = paste(
                  round(call$cn,  2), 
                  round(call$min, 2), 
                  round(call$max, 2),
                  round(call$qp,  0),
                  format(call$p.val, digits=2),
                  sep=":"
              )
          )
      )
    } else if (call$score >= config.score.check) {
      vcf <- rbind(
          vcf, data.frame(
              CHROM   = call$chr, #paste('chr', call$chr, sep=""),
              POS    = call$from, 
              ID      = '.',
              REF     = refline$base[1],
              ALT      = ifelse( call$cn < call$refcopy , '<DEL>', '<DUP>'),
              QUAL     = round(call$score),
              FILTER   = call$filter,
              INFO     = paste(
                  'IMPRECISE', ';',
                  'SVTYPE=', ifelse( call$cn < call$refcopy , 'DEL', 'DUP'), ';',
                  'END=', call$to, ';',
                  'SVLEN=', ifelse( call$cn < call$refcopy , -bp, bp), ';',
                  'LOCUS=', call$genename,
                  sep=""
              ),
              FORMAT    = 'CN:CNMIN:CNMAX:CNQ:PVAL',
              SAMPLE    = paste(
                  round(call$cn,  2), 
                  round(call$min, 2), 
                  round(call$max, 2),
                  round(call$qp,  0),
                  format(call$p.val, digits=2),
                  sep=":"
              )
          )
      )
    } else {
      
      info <- paste('UNCHANGED;END=', call$to,';LOCUS=', call$genename, sep="")
      if (!call$sufficient) {
        info <- paste('NOCALL;END=', call$to,';LOCUS=', call$genename, sep="")
      }
      vcf <- rbind(
          vcf, data.frame(
              CHROM   = call$chr, #paste('chr', call$chr, sep=""),
              POS     = call$from, 
              ID      = '.',
              REF     = refline$base[1],
              ALT     = '.',
              QUAL    = round(call$score),
              FILTER  = call$filter,
              INFO    = info,
              FORMAT  = 'CN:CNMIN:CNMAX:CNQ:PVAL',
              SAMPLE  = paste(
                  round(call$cn,  2), 
                  round(call$min, 2), 
                  round(call$max, 2),
                  round(call$qp,  0),
                  format(call$p.val, digits=2),
                  sep=":"
              )
          )
      )
      }
    } else { ## segmented
      
        for (region in 1:nrow(segments)) {

          bp      <- segments$right[region] - segments$left[region]
          refline <- subset(depth, chr == call$chr & position == segments$left[region] )
          info    <-  paste('UNCHANGED;SEGMENTED;END=', segments$right[region]-1,';LOCUS=', call$genename, '/', region, sep="")
          alt     <- '.'

          if (segments$score[region] >=  config.score.check) {
            alt  <- ifelse( segments$cn[region] < call$refcopy , '<DEL>', '<DUP>')
            info <- paste(
                'IMPRECISE;SEGMENTED', ';',
                'SVTYPE=CNV', ';',
                'END=', segments$right[region], ';',
                'SVLEN=', ifelse( segments$cn[region] < call$refcopy , -bp, bp), ';',
                'LOCUS=', call$genename, '/', region,
                sep=""
            )
          }
          vcf <- rbind(
              vcf, data.frame(
                  CHROM   = call$chr, #paste('chr', call$chr, sep=""),
                  POS    = segments$left[region], 
                  ID      = '.',
                  REF     = refline$base[1],
                  ALT      = alt,
                  QUAL     = round(min(c(config.score.max, segments$score[region]))),
                  FILTER   = set_filter(segments$score[region], segments$usable[region], config.score.report, config.score.check, config.minimal_observations),
                  INFO     = info,
                  FORMAT   = 'CN:CNMIN:CNMAX:CNQ:PVAL',
                  SAMPLE   = paste(
                        round(segments$cn[region],  2), 
                        round(segments$min[region], 2), 
                        round(segments$max[region], 2), 
                        round(segments$qp[region],  0),                         
                        format(segments$p.val[region], digits=2),
                        sep=":"
                  ) 
              )
          )
          genes.all <- rbind(
              genes.all, data.frame(
                  genename    = sprintf("%s/%d", call$genename, region),
                  chr         = call$chr,
                  from        = segments$left[region],
                  to          = segments$right[region],
                  refcopy     = call$refcopy,
                  amplicons   = segments$total_amplicons[region],
                  outliers    = 0, # untracked so far
                  usable      = segments$total_amplicons[region],
                  sufficient  = TRUE,
                  alternative = NA,
                  p.val       = segments$p.val[region],
                  conf.level  = config.conf_level,
                  score       = round(min(c(config.score.max, segments$score[region]))),
                  filter      = set_filter(segments$score[region], segments$usable[region], config.score.report, config.score.check, config.minimal_observations),
                  cn          = segments$cn[region],
                  min         = segments$min[region],
                  max         = segments$max[region],
                  se          = segments$se[region],
                  qp          = segments$qp[region],
                  segmented   = TRUE,
                  seq.p.value = NA,
                  segments    = NA,
                  cn.color    = NA,
                  cn.label    = NA,
                  sg.color    = NA,
                  sg.label    = NA,
                  sg.pval     = NA,
                  report      = TRUE
              )
          )
            
        }
    }
  }
}

if (counter == 2) {
  # add the last plot to the pdf (if any)
  multiplot(plotlist = plots, layout=config.report.layout)
}

if (nrow(vcf) > 0) {
  
# export the vcf
writeLines(
  paste(vcf_header(), '
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=0,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=UNCHANGED,Number=0,Type=Flag,Description="No copy number variation in this region">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=SEGMENTED,Number=0,Type=Flag,Description="Structural variation segments named locus group (e.g. gene)">
##INFO=<ID=NOCALL,Number=0,Type=Flag,Description="Copy number analysis was not possible (see FILTER for reason)">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=LOCUS,Number=1,Type=String,Description="Locus description (e.g. gene name in this region)">
##ALT=<ID=CNV,Description="Copy number variable region">
##ALT=<ID=DEL,Description="Deletion relative to reference (decreased copy number, loss)">
##ALT=<ID=DUP,Description="Duplication relative to reference (increased copy number, gain)">
##FILTER=<ID=q', config.score.check ,',Description="Region not changed in copy number: Score below ', config.score.check,'">
##FILTER=<ID=q', config.score.report ,',Description="Region potentially changed in copy number: Score below ', config.score.report,'">
##FILTER=<ID=o', config.minimal_observations , ',Description="Less than ' , config.minimal_observations ,' amplicons have suffcient read depth (<' , config.min_cov_for_call , ')">
##FORMAT=<ID=CN,Number=1,Type=Float,Description="Estimated copy number for imprecise events">
##FORMAT=<ID=CNMIN,Number=1,Type=Float,Description="Minimal copy number supported by statistical test (conf=', config.conf_level ,')">
##FORMAT=<ID=CNMAX,Number=1,Type=Float,Description="Maximal copy number supported by statistical test (conf=', config.conf_level ,')">
##FORMAT=<ID=CNQ,Number=1,Type=Integer,Description="Copy number quality score">
##FORMAT=<ID=PVAL,Number=1,Type=Float,Description="Probability for normal copy number">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE',
  sep=""),
  paste(output.dir , '/', output.basename, '.vcf', sep="")
)

# sort the vcf table by chromosome / position
vcf <- vcf[naturalorder(paste(vcf$CHROM, vcf$POS)),]
vcf$CHROM <- paste('chr', vcf$CHROM, sep="")

# export the vcf to a file
write.table(
    vcf,
    paste(output.dir , '/', output.basename, '.vcf', sep=""),
    sep       = "\t", 
    quote     = FALSE,
    append    = TRUE,
    row.names = FALSE,
    col.names = FALSE
)
}

## export excel file (tab separated for now)
##
## first select columns to include in export:
## available:
## genename	chr	from	to	refcopy	amplicons	outliers	usable	noise	sufficient	noisy	alternativep.val	conf.level	score	filter	cn	min	max	se	ml.cn	ml.prob	segmented	seq.p.value	segments	cn.color	cn.label	sg.color	sg.label	sg.pval	report
if (nrow(genes.all)> 0) {

  genes.all$readSet <- output.basename
  export <- genes.all[c("readSet", "chr", "from", "to", "genename",  "amplicons", "outliers", "usable", "refcopy", "log2", "cn", "min", "max", "qp", "sd", "score", "p.val", "filter")]
  names(export) <- c("readSet","chromosome", "start", "end", "locus", "amplicons", "outliers", "usable", "expected", "log2", "copies", "min", "max", "qp", "sd", "score", "p.val", "filter")
  
  # round to integer
  export$score <- round(export$score, 0)
  export$qp    <- round(export$qp, 0)
  
  # round to two decimal
  export$log2   <- round(export$log2,   2)
  export$copies <- round(export$copies, 2)
  export$min    <- round(export$min, 2)
  export$max    <- round(export$max, 2)
  export$sd     <- round(export$sd, 2)
  
  # format p.values
  export$p.val <- format(export$p.val, digits=2)
  
  # remove NA
  export$p.val[is.na(export$copies)]  <- ''
  export$score[is.na(export$copies)]  <- ''
  export$min[is.na(export$copies)]    <- ''
  export$max[is.na(export$copies)]    <- ''
  export$log2[is.na(export$copies)]   <- ''
  export$sd[is.na(export$copies)]     <- ''
  export$copies[is.na(export$copies)] <- ''
  
  # sort by chromosome and position
  export            <- export[naturalorder(paste(export$chromosome, export$start)),] 
  export$chromosome <- paste('chr', export$chromosome, sep="")
  
  # save the subset
  write.table(
      export,
      paste(output.dir , '/', output.basename, '.csv', sep=""), 
      sep=",", 
      quote = FALSE,
      row.names = FALSE
  )

}

# create the overview with final reported ranges (separate page)
plot_overview(
    csv      = paste(output.dir , '/', output.basename, '.csv', sep=""), 
    filename = paste(output.dir , '/', output.basename, '_overview.pdf', sep="")
)

# save the full dataset (for debugging)
if (config.write.dataset) {
  write.table(
      depth,
      paste(output.dir , '/', output.basename, '_full_data.tsv', sep=""), 
      sep="\t", 
      quote = FALSE,
      row.names = FALSE
  )
}

#message("Done")

dev.off()
summary(warnings())
q()
