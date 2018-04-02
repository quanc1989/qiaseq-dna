import os
import string
import subprocess
import imp
from collections import defaultdict

# 3rd party modules
import pysam
import editdist

# constants for read pair accounting
NUM_PRIMER_SIDE_NOT_MAPPED = 0
NUM_RANDOM_SIDE_NOT_MAPPED = 1
NUM_R1_R2_NOT_AT_SAME_LOCUS = 2
NUM_R1_R2_SAME_ORIENTATION = 3
NUM_SPLIT_ALIGNMENT = 4
NUM_LOW_MAPQ = 5
NUM_LT_25BP_ALIGNED = 6
NUM_TOO_MUCH_SOFTCLIP = 7
NUM_NO_PRIMER_BUT_NEAR_SPE_PRIMING_SITE = 8
NUM_NO_PRIMER_NOT_NEAR_SPE_PRIMING_SITE = 9
NUM_ENDOGENOUS_BP_ALIGNED_TOO_LOW = 10
#NUM_PRIMER_FOUND_VIA_EDIT_DIST = 11 
#NUM_PRIMER_FOUND_VIA_SW_SHORTCUT = 12
NUM_PRIMER_AT_DESIGN_SITE = 13 - 2
NUM_PRIMER_NOT_AT_DESIGN_SITE = 14 - 2
NUM_METRICS_TOTAL = 15 - 2

#---------------------------------------------------------------------
# reverse complement a seq (warning: will not work with Python 3)
#---------------------------------------------------------------------
dnaComplementTranslation = string.maketrans("ATGC", "TACG")
def reverseComplement(seq):
   seq = seq[::-1]
   return seq.translate(dnaComplementTranslation)

#---------------------------------------------------------------------
# main function
#---------------------------------------------------------------------
def run(cfg,bamFileIn):
   print("umi_filter: starting...")

   # get parameters
   readSet          = cfg.readSet
   primerFile       = cfg.primerFile
   endogenousLenMin = int(cfg.endogenousLenMin)
   tagNameUmiSeq    = cfg.tagNameUmiSeq   
   deleteLocalFiles = cfg.deleteLocalFiles
   numCores         = cfg.numCores
   primer3Bases     = int(cfg.primer3Bases)
   maxSoftClipBp1   = min(primer3Bases + 4, 20) if primer3Bases != -1 else 20      # 16 + 4 = 20
   tagNamePrimer    = cfg.tagNamePrimer
   tagNamePrimerErr = cfg.tagNamePrimerErr

   # set output file prefix
   filePrefixOut = readSet + ".umi_filter"

   # put primer seqs in two dictionaries, one for each read strand
   primerSeq = {}
   primerDicts = (defaultdict(list), defaultdict(list))
   for line in open(primerFile, "r"):
      (chrom, loc3, direction, primer) = line.strip().split("\t")
      primerStrand = 0 if direction == "L" or direction == "0" else 1
      loc3 = int(loc3)
      key = (chrom,primerStrand,loc3)
      
      # debug check (need to modify code to handle same primer for multiple design loci)
      if key in primerSeq:
         raise Exception("ERROR: duplicate primer specification! " + primer)

      # hash primer sequence by site as key
      primerSeq[key] = primer      
         
   print("# of primers:", len(primerSeq))
   
   # open output files
   fileout         = open(filePrefixOut + ".alignments.txt", "w")
   fileoutNoPrimer = open(filePrefixOut + ".no-primer.txt" , "w")

   # open BAM read alignment file
   bam = pysam.Samfile(bamFileIn, "rb")

   # loop over read alignments
   primingSitesApprox = {}
   primingReadDepths = {}
   readPairCounts = [0] * NUM_METRICS_TOTAL
   for read in bam:   
      # this is dangerous, but drop these for now
      if read.is_secondary or read.is_supplementary:
         continue

      # crash if read is not paired
      if not read.is_paired:
         print(read.qname)
         raise Exception("read not paired!")
      
      # this should be R1      
      read1 = read
      
      # get mate, assuming mate is the next record in the BAM file
      while True:
         read = bam.next()
         if not read.is_secondary and not read.is_supplementary:
            break
            
      # this should be R2
      read2 = read
      
      # debug check
      if read1.qname != read2.qname:
         print(read1.qname, read2.qname)
         raise Exception("read mate is not next in BAM record order!")
         
      # debug check 
      if not read1.is_read1 or not read2.is_read2:
         raise Exception("R1/R2 mixed up!")
     
      # skip but count unmapped R1 reads, even if R2 mapped.  Need to look at these later...
      if read1.is_unmapped:
         readPairCounts[NUM_PRIMER_SIDE_NOT_MAPPED] += 1
         continue
      
      # skip but count unmapped R2 reads
      if read2.is_unmapped:
         readPairCounts[NUM_RANDOM_SIDE_NOT_MAPPED] += 1
         continue
         
      # skip reads not mapped to same chrom
      chrom1 = bam.getrname(read1.tid)       # getrname is depricated, read1.reference_name can be used instead
      chrom2 = bam.getrname(read2.tid)
      if chrom1 != chrom2:
         readPairCounts[NUM_R1_R2_NOT_AT_SAME_LOCUS] += 1
         continue
         
      # skip reads not mapped to same locus
      locRead1 = int(read1.aend) - 1 if read1.is_reverse else read1.pos
      locRead2 = int(read2.aend) - 1 if read2.is_reverse else read2.pos
      if abs(locRead1 - locRead2) > 2000:
         readPairCounts[NUM_R1_R2_NOT_AT_SAME_LOCUS] += 1
         continue

      # skip pairs with odd alignment orientation
      if read1.is_reverse == read2.is_reverse:
         readPairCounts[NUM_R1_R2_SAME_ORIENTATION] += 1
         continue
      
      # drop read pair if either end has a supplementary split alignment
      if read1.has_tag("SA") or read2.has_tag("SA"):
         readPairCounts[NUM_SPLIT_ALIGNMENT] += 1
         continue
            
      # drop read pair if R1 or R2 read has low mapq 
      if read2.mapq < 17 or read1.mapq < 17:
         readPairCounts[NUM_LOW_MAPQ] += 1
         continue

      # get 5' soft clips on both reads;  also adjust alignLoc on random R2 end
      if read2.is_reverse:
         # R1 soft clip
         (cigarType, cigarLen) = read1.cigar[0]
         softClip1 = cigarLen if cigarType == 4 else 0
         
         # R2 soft clip
         (cigarType, cigarLen) = read2.cigar[-1]
         softClip2 = cigarLen if cigarType == 4 else 0
            
         # R2 pseudo-start
         alignLocRand = read2.aend - 1 + softClip2
      else:
         # R1 soft clip
         (cigarType, cigarLen) = read1.cigar[-1]
         softClip1 = cigarLen if cigarType == 4 else 0
      
         # R2 soft clip
         (cigarType, cigarLen) = read2.cigar[0]
         softClip2 = cigarLen if cigarType == 4 else 0
            
         # R2 pseudo-start
         alignLocRand = read2.pos - softClip2
         
      # require some significant alignment to genome
      if read2.aend - read2.pos < 25 or read1.aend - read1.pos < 25:
         readPairCounts[NUM_LT_25BP_ALIGNED] += 1
         continue

      # get primer info from read tag
      primerErrBp = int(read1.get_tag(tagNamePrimerErr))
      
      # drop read if a lot of R2 soft clipping - makes barcode clustering problematic.  also minimial R1 primer side soft clipping. NOTE: changed threshold based on primer3Bases parameter
      if softClip2 > 3 or softClip1 > (maxSoftClipBp1 - primerErrBp):
         readPairCounts[NUM_TOO_MUCH_SOFTCLIP] += 1
         continue
         
      # get UMI sequence and primer info from "mi" and "pi" tags
      readId = read1.qname
      umiSeq = read1.get_tag(tagNameUmiSeq)
      primerInfo = read1.get_tag(tagNamePrimer)
      if primerInfo == "-1":
         primer = None
         primerLen = 0
      else:
         (chrom,primerStrand,loc3) = primerInfo.split("-")
         loc3 = int(loc3)
         primerStrand = int(primerStrand)
         primer = primerSeq[(chrom,primerStrand,loc3)]
         primerLen = len(primer)
         loc5 = loc3 - primerLen + 1 if primerStrand == 0 else loc3 + primerLen - 1         
      
      # bases to adjust for removed bases from 5' of primer
      primerOffset = primerLen if (primer3Bases == -1 or primer3Bases > primerLen) else primer3Bases
      if primerStrand == 1:
         primerOffset = -primerOffset
      
      # get some alignment info from R1 read (SPE primer side)
      readSeq     = read1.seq
      alignChrom  = chrom1
      alignCigar  = read1.cigar
      alignStrand = 1 if read1.is_reverse else 0
      if alignStrand == 0:
         alignLoc = read1.pos
      else:
         alignLoc = read1.aend - 1 #  0-based position of the 5' end of the read
         alignCigar.reverse()
      
      # check if near a primer site     
      primer_ = None
      for offset in range(-5,6):    #(0,1,-1):
         designLoc = alignLoc + offset + primerOffset
         key = (alignChrom, alignStrand, designLoc)
         if key in primerSeq:
            primer_ = primerSeq[key]               
            break
      
      # if primer not found
      if primer == None:               
         if primer_ == None:
            isNearSpePriming = 0
            readPairCounts[NUM_NO_PRIMER_NOT_NEAR_SPE_PRIMING_SITE] += 1         
         else:
            isNearSpePriming = 1
            readPairCounts[NUM_NO_PRIMER_BUT_NEAR_SPE_PRIMING_SITE] += 1   

         # save alignment to disk
         outvec = (alignChrom, alignStrand, alignLoc, isNearSpePriming, umiSeq, readId, read1.pos, read1.aend, read1.cigarstring, read2.pos, read2.aend, read2.cigarstring)
         outvec = (str(x) for x in outvec)
         fileoutNoPrimer.write("|".join(outvec))
         fileoutNoPrimer.write("\n")
            
         # skip to next read pair         
         continue          
      
      # drop read pair if R2 not aligned to 15 bp beyond the primer (assuming primer does not have huge indel bubbles)
      if alignStrand == 1:
         endogenousLen = alignLoc - read2.pos + 1 - len(primer)
      else:
         endogenousLen = read2.aend - alignLoc - len(primer)
      if endogenousLen < endogenousLenMin:
         readPairCounts[NUM_ENDOGENOUS_BP_ALIGNED_TOO_LOW] += 1
         continue     
      
      # check if priming was at the intended design site or not (note: possible impedance mismatch between the <5 criteria and the primer finding stringency!)
      #if alignChrom == chrom and alignStrand == primerStrand and abs(alignLoc - loc5) < 5:
      if primer_ == primer:
         isIntendedSite = 1
         readPairCounts[NUM_PRIMER_AT_DESIGN_SITE] += 1
      else:
         isIntendedSite = 0
         readPairCounts[NUM_PRIMER_NOT_AT_DESIGN_SITE] += 1

      # get cigar strings for both ends
      cigar1 = "*" if read1.cigarstring == None else read1.cigarstring
      cigar2 = "*" if read2.cigarstring == None else read2.cigarstring
      
      # write output (NOTE: field position hard coded in Linux sort below!)
      outvec = (chrom, loc5, primerStrand, primer, umiSeq, isIntendedSite, alignChrom, alignStrand, alignLocRand, alignLoc, readId, read1.pos, read1.aend, cigar1, read2.pos, read2.aend, cigar2)
      fileout.write("|".join((str(x) for x in outvec)))
      fileout.write("\n")      
      
      # count read depths for debug of multiple primers priming nearby
      key = (alignChrom, alignStrand, alignLoc, primer)
      if key not in primingReadDepths:
         primingReadDepths[key] = 1
      else:
         primingReadDepths[key] += 1

   # done with bam input and the alignment output files
   bam.close()
   fileout.close()
   fileoutNoPrimer.close()
   
   # sort the priming site read depths by locus
   primingDepths = []
   for (key, numReads) in primingReadDepths.iteritems():
      (alignChrom, alignStrand, alignLoc, primer) = key
      primingDepths.append((alignChrom, alignStrand, alignLoc, primer, numReads))
   primingDepths.sort()
   
   # output warnings for multiple primers priming near each other
   fileout = open(filePrefixOut + ".warnings.txt", "w")
   for i in range(len(primingDepths)):
      (alignChrom, alignStrand, alignLoc, primer, numReads) = primingDepths[i]
      j = i - 1
      while True:
         if j < 0:
            break
         (alignChrom_, alignStrand_, alignLoc_, primer_, numReads_) = primingDepths[j]
         if alignChrom_ != alignChrom or alignStrand_ != alignStrand:
            break
         if abs(alignLoc_ - alignLoc) > 5:
            break
         if primer_ != primer:
            outvec = []
            outvec.extend(primingDepths[i])
            outvec.extend(primingDepths[j])
            outvec = (str(x) for x in outvec)
            fileout.write("|".join(outvec))
            fileout.write("\n")
         j -= 1
   fileout.close()
   
   # compute read pairs dropped, for short version of summary
   readPairsDropped = (readPairCounts[NUM_R1_R2_NOT_AT_SAME_LOCUS]
                    +  readPairCounts[NUM_R1_R2_SAME_ORIENTATION]
                    +  readPairCounts[NUM_SPLIT_ALIGNMENT]
                    +  readPairCounts[NUM_LOW_MAPQ]
                    +  readPairCounts[NUM_LT_25BP_ALIGNED]
                    +  readPairCounts[NUM_TOO_MUCH_SOFTCLIP]
                    +  readPairCounts[NUM_ENDOGENOUS_BP_ALIGNED_TOO_LOW]
                    +  readPairCounts[NUM_NO_PRIMER_BUT_NEAR_SPE_PRIMING_SITE]
                    +  readPairCounts[NUM_NO_PRIMER_NOT_NEAR_SPE_PRIMING_SITE])
                    
   # compute on-target SPE priming specificity (does NOT include universal priming, etc!)
   readPairsPrimerFound = readPairCounts[NUM_PRIMER_AT_DESIGN_SITE] \
                        + readPairCounts[NUM_PRIMER_NOT_AT_DESIGN_SITE]
   readPairsPrimerFoundOnTargetPct = 100.000 * readPairCounts[NUM_PRIMER_AT_DESIGN_SITE] / readPairsPrimerFound if readPairsPrimerFound > 0 else 0.00

   # compute on-target SPE priming specificity as a percent of all reads not dropped after universal trimming
   readPairsTotal = readPairCounts[NUM_PRIMER_SIDE_NOT_MAPPED] + readPairCounts[NUM_RANDOM_SIDE_NOT_MAPPED] + readPairsDropped +  + readPairsPrimerFound
   readPairsPrimerFoundOnTargetPctOfAll = 100.000 * readPairCounts[NUM_PRIMER_AT_DESIGN_SITE] / readPairsTotal
   
   # report read accounting, detail version
   fileout = open(filePrefixOut + ".detail.summary.txt", "w")
   fileout.write("{}\tread fragments dropped, primer side read not mapped\n".format(readPairCounts[NUM_PRIMER_SIDE_NOT_MAPPED]))
   fileout.write("{}\tread fragments dropped, random side read not mapped\n".format(readPairCounts[NUM_RANDOM_SIDE_NOT_MAPPED]))
   fileout.write("{}\tread fragments dropped, R1 and R2 not mapped to same locus\n".format(readPairCounts[NUM_R1_R2_NOT_AT_SAME_LOCUS]))
   fileout.write("{}\tread fragments dropped, FF or RR mapping orientation\n".format(readPairCounts[NUM_R1_R2_SAME_ORIENTATION]))
   fileout.write("{}\tread fragments dropped, split alignment\n".format(readPairCounts[NUM_SPLIT_ALIGNMENT]))
   fileout.write("{}\tread fragments dropped, low mapping quality MAPQ < 17\n".format(readPairCounts[NUM_LOW_MAPQ]))
   fileout.write("{}\tread fragments dropped, < 25 bp aligned to genome\n".format(readPairCounts[NUM_LT_25BP_ALIGNED]))
   fileout.write("{}\tread fragments dropped, > 3 bp soft clip on random side, or > 20 bp soft clip on primer side\n".format(readPairCounts[NUM_TOO_MUCH_SOFTCLIP]))
   fileout.write("{}\tread fragments dropped, < 15 bp endogenous after primer\n".format(readPairCounts[NUM_ENDOGENOUS_BP_ALIGNED_TOO_LOW]))
   fileout.write("{}\tread fragments dropped, no primer match, but near an SPE priming site\n".format(readPairCounts[NUM_NO_PRIMER_BUT_NEAR_SPE_PRIMING_SITE]))
   fileout.write("{}\tread fragments dropped, no primer match, not near an SPE priming site\n".format(readPairCounts[NUM_NO_PRIMER_NOT_NEAR_SPE_PRIMING_SITE]))
   fileout.write("{}\tread fragments dropped, off target\n".format(readPairCounts[NUM_PRIMER_NOT_AT_DESIGN_SITE]))
   fileout.write("{}\tread fragments with primer found, on-target\n".format(readPairCounts[NUM_PRIMER_AT_DESIGN_SITE]))
   fileout.write("{0:.2f}\tread fragments with primer found, on-target percent\n".format(readPairsPrimerFoundOnTargetPct))
   fileout.write("{0:.2f}\tread fragments with primer found, on-target percent of all\n".format(readPairsPrimerFoundOnTargetPctOfAll))
   #fileout.write("{}\tread fragments with primer found via edit distance to intended primer\n".format(readPairCounts[NUM_PRIMER_FOUND_VIA_EDIT_DIST]))
   #fileout.write("{}\tread fragments with primer found via Smith-Waterman shortcut\n".format(readPairCounts[NUM_PRIMER_FOUND_VIA_SW_SHORTCUT]))
   fileout.close()

   # report read accounting, brief version
   fileout = open(filePrefixOut + ".summary.txt", "w")
   fileout.write("{}\tread fragments dropped, R1 or R2 not mapped\n".format(readPairCounts[NUM_PRIMER_SIDE_NOT_MAPPED] + readPairCounts[NUM_RANDOM_SIDE_NOT_MAPPED]))
   fileout.write("{}\tread fragments dropped, not passing mapping filters (low mapq, split alignments, discordant pairs, etc.)\n".format(readPairsDropped))
   fileout.write("{}\tread fragments dropped, off target\n".format(readPairCounts[NUM_PRIMER_NOT_AT_DESIGN_SITE]))
   fileout.write("{}\tread fragments with primer found, on-target\n".format(readPairCounts[NUM_PRIMER_AT_DESIGN_SITE]))
   fileout.write("{0:.2f}\tread fragments with primer found, on-target percent\n".format(readPairsPrimerFoundOnTargetPct))
   fileout.close()
   
   # sort the "no-primer" alignment file by locus (helpful for debug)
   fileName = filePrefixOut + ".no-primer.txt"
   cmd = "sort -k1,1 -k2,2n -k3,3n -k4,4n -k5,5n -t\| -T./ --parallel={1} {0} > {0}.tmp".format(fileName,numCores)
   subprocess.check_call(cmd, shell=True)
   os.rename(fileName + ".tmp", fileName)
   
   # stop pipeline if no reads
   if readPairCounts[NUM_PRIMER_AT_DESIGN_SITE] == 0:
      raise UserWarning("Zero on-target reads found for read set: " + readSet)

   # report completion
   print("umi_filter: done")
