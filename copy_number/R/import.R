###############################################################################
# import_counts
#
# This function reads two tables with locations and counts, merged these into 
# a data frame and also sets the reference copy number based on the provided 
# copy numbers for X and Y in reference and sample
#
# input:      sample      = filename with sample data
#             sample.x    = expected counts for X chromosome in sample
#             sample.y    = expected counts for Y chromosome in sample
#             reference   = filename with reference data
#             reference.x = real counts (may be float) of X in reference
#             reference.y = real counts (may be float) of Y in reference
#
# output:     data.frame
#

import_counts <- function(
      sample      = NA, 
      reference   = NA, 
      sample.x    = config.default.sample.x, # default to female sample
      sample.y    = config.default.sample.y, # default to female sample
      reference.x = config.default.reference.x, # default to female reference
      reference.y = config.default.reference.y  # default to female reference
  ) {
  # check if filenames are passed
  if(is.na(sample) | is.na(reference)) {
    warn("Sample or reference input file not provided, stopping here.")
    q()
  }
  # read both files
  smp <- read.table(sample, header = FALSE)    # data for sample
  ref <- read.table(reference, header = FALSE) # data for reference
  
  # name columns (just to easy step below, and for sanity check debugging)
  colnames(smp) <- c("chr_full", "position", "direction", "base", "gene", "sample", "orgname")
  colnames(ref) <- c("chr_full", "position", "direction", "base", "gene", "reference", "orgname")
  
  # merge
  depth <- merge(smp, ref, by = intersect(names(smp), names(ref)))
  
  # replace chr to make a short version of the chromosome
  depth$chr <- gsub("^chr", "", depth$chr_full)
  
  # preset refcopy number to ploidity default
  depth$refcopy <- config.default.ploidity
  
  # treat the sex chromosomes differently, depending on input. Here we 
  # will use the counts of the *sample* (although is says reference copies)
  # this is confusing, but this will be used as "expected" copy number for
  # reporting, and we expect X=1 and Y=1 for a male, even if the reference
  # is from a female!
  depth$refcopy[depth$chr == 'X'] <- sample.x # expected X copies in sample
  depth$refcopy[depth$chr == 'Y'] <- sample.y # expected Y copies in sample
  
  # calculate correction factors for x and y (for non-matched or ref-pools)
  correct.x <- ifelse(reference.x > 0, sample.x / reference.x, 0)
  correct.y <- ifelse(reference.y > 0, sample.y / reference.y, 0)
  
  # correct counts by reducing the counts that are higher than needed
  # (compared to a a truely matched reference). this is to avoid scaling up, 
  # which may add artificial precision (higher count) than it actually has
  if (correct.x < 1) {
    # X in reference is higher than in sample, reduce reference to match
    depth$reference[depth$chr == 'X'] <- depth$reference[depth$chr == 'X'] * correct.x
  }  else if (correct.x > 1) {
    # X in sample is higher than in reference, reduce sample to match
    depth$sample[depth$chr == 'X'] <- depth$sample[depth$chr == 'X'] / correct.x
  }
  if (correct.y < 1) {
    # Y in reference is higher than in sample, reduce reference to match
    depth$reference[depth$chr == 'Y'] <- depth$reference[depth$chr == 'Y'] * correct.y
  } else if (correct.y > 1) {
    # Y in sample is higher than in reference, reduce sample to match
    depth$sample[depth$chr == 'Y'] <- depth$sample[depth$chr == 'Y'] / correct.x
  }
  
#  
## change the expected copy number
#  if (sample.sex == 'male') {
#    depth$refcopy[depth$chr == 'X'] <- 1
#    depth$refcopy[depth$chr == 'Y'] <- 1
#    if (reference.sex == 'female') {
#      # do we need to correct for too many X reads here?
#      depth$reference[depth$chr == 'X'] <- depth$reference[depth$chr == 'X'] * 0.5
#    }
#  } else {
#    depth$refcopy[depth$chr == 'Y'] <- 0
#    if (reference.sex == 'male') {
#      # do we need to correct for missing X reads here?
#      depth$reference[depth$chr == 'X'] <- depth$reference[depth$chr == 'X'] * 2
#    }
#  }
  
  # return the data.frame
  depth <- depth
  return(depth)
}

## EXAMPLE INPUT FILE FORMAT:
##
##  chr1	5923229	0	A	NPHP4	1069
##  chr1	5923341	1	T	NPHP4	1061
##  chr1	5923335	0	A	NPHP4	486
##  chr1	5923440	1	A	NPHP4	491
##  chr1	5923429	0	A	NPHP4	1064
##  chr1	5923549	1	A	NPHP4	1122
##  chr1	5923905	0	T	NPHP4	1366
