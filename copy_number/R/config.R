## // GLOBAL CONFIGURATION SETTINGS //
##
config.default.ploidity    <- 2 # default to diploid
config.default.sample.x    <- 2 # default sample to be female
config.default.sample.y    <- 0 # default sample to be female
config.default.reference.x <- 2 # default reference to be female
config.default.reference.y <- 0 # default reference to be female

## // DATA READ FROM EXTERNAL FILES AND VERSIONING STUFF
##
## internal
config.software.name       <- 'quandico'
config.software.version    <- '1.0.2'
## reference assembly
#config.assembly.filename   <- './hg19.fa.fai'
#config.assembly.filename   <- '/srv/qgen/data/genome/ucsc.hg19.fasta.fai'
config.assembly.version    <- 'hg19'

## // SETTINGS AFFECTING NORMALIZATION
##
## absolute coverage required in either ref or sample, that means,
## count data where both reference and sample are below, will not be used
config.min_cov_for_scale     <- 20 ## absolute coverage for normalization
config.min_cov_for_call      <- 20 ## absolute coverage for calls
config.percentile            <- 50 ## fraction of "neutral" amplicons (TODO: unused)

## // SETTINGS AFFECTING TREATMENT OF OUTLIERS AND NOISY DATA
##
## require this many amplicons per region (gene) to justify statistic testing,
## that means, at least this number of count-values must be present per 
## region (gene); this will be checked *after* removal of outliers

config.conf_level            <-  0.95
config.log2.max              <-  6.0 # max log2 (set Inf to this)
config.log2.min              <- -6.0 # min log2 (set NaN to this)

config.minimal_observations  <- 10                  # 10
config.use.weights           <- 'sum(raw)' # FALSE = off TODO: SET!
config.max.noise             <-  1

# outlier removal strategies
config.narrow.normal         <- TRUE # TRUE  TODO: SET!
config.narrow.test           <- 'shapiro'
config.narrow.pval           <- 0.05 # remove outlier untill p >= pval
config.narrow.lownoise       <- 1/3  # do not remove value if sd <= lownoise
config.narrow.min            <- 2/3  # keep this fraction of all values in any case

config.flag.outliers         <- FALSE # TRUE  TODO: SET!
config.recurse.outliers      <- FALSE
config.flag.extremes         <- FALSE

## // SETTINGS AFFECTING SENSITIVITY / SPECIFICITY => REPORTING
##
## report all copy numbers with score above the report threshold
config.score.report   <-  50 ## report above this score             (50)
config.score.check    <-  25 ## report uncertain above this score   (25)
config.score.max      <- 250 ## maximal score (hard clip)          (250)

## // SETTINGS AFFECTING THE SCORE CALCULATIONS
##
## post-process intial score from p.values
## supported methods: NA (off), divseq, log2, sigma <- last is recommended
config.score.correction       <- 'sqrt(sd)*log2' ## method to modify scores # TODO: SET!!!

## // SETTINGS AFFECTING AMOUNT OF REPORTING
##
## report output
config.report.layout         <- matrix(c(1,2,2,2,3,4,4,4), nrow = 2, byrow = TRUE)
config.report.pdf.paper      <- 'US'
config.report.pdf.colormodel <- 'cmyk'
##
## what to report
config.report.segmented       <- FALSE ## report segmented            (FALSE)
config.report.insufficient    <- TRUE  ## report insufficient          (TRUE)
config.report.noisy           <- TRUE  ## report noisy                 (TRUE)
config.report.all             <- TRUE  ## originally for dvlpmt phase  (TRUE)
config.report.genes           <- FALSE ## save a pdf for single genes (FALSE)
##
## dump data for debugging
config.write.dataset          <- FALSE  ## saves the complete data     (FALSE)

## // SETTINGS AFFECTING SEGMENTATION
##
config.segment.min.amplicons  <- 10.0
config.segment.p.value        <-  0.7
config.segment.diff.pval      <-  1.0e-8
config.segment_max_segment_sd <-  1.0
config.segment.rm.out         <- TRUE  ## remove outliers in segments

## // SETTINGS AFFECTING OPTIONAL REPORTS
##
## plots (on / off)
config.plot.rawdata      <- FALSE ## plot raw log2 counts              (FALSE)
config.plot.normalized   <- FALSE ## plot notrmalized log2 counts      (FALSE)
config.plot.overview     <- TRUE  ## plot summary per region/locus     (TRUE)
config.plot.chrom_boxes  <- FALSE ## generate boxplot per chromosome   (FALSE)
config.plot.chromosomes  <- FALSE ## generate scaled chromosome-grpahs (FALSE)

## // SETTINGS AFFECTING THE LAYOUTS OF PLOTS (also see file: themes.R)
##
## scales
config.plots.max.cn          <- 12
config.plots.ylimits         <- c(0, config.plots.max.cn + 0.9)
config.plots.ybreaks         <- c(0, 1, 2, 3, 4, 5, 6, 8, 10, 12)
config.plots.ylabels         <- c("0", "1", "2", "3", "4", "5", "6", "8", "10", "12")
config.plots.line1           <- 12.25
config.plots.line2           <- 11.50
config.plots.seg.lcolor      <- "orange"
config.plots.seg.ltype       <- "solid"
config.plots.alpha.violin    <-  0.50
config.plots.smooth.violin   <-  1.75
##
## single gene output
config.gene.layout           <- matrix(c(1,2,2,2), nrow = 1, byrow = TRUE)
config.gene.pdf.paper        <- 'USr'
##
## png graphics
config.png.width            <- 1200
config.png.height           <- 600
##
## pdf 
##
## <-- nothing here so far -->
##
## colors
config.color.noisy           <- "brown"
config.color.even            <- "green"
config.color.unchanged       <- "green"
config.color.segmented       <- "orange"
config.color.insufficient    <- "gray50"
config.color.check           <- "gray50"
config.color.untested        <- "gray50"
