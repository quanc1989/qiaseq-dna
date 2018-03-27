import math
from collections import defaultdict

#----------------------------------------------------------------------------
# get read-per-barcode metrics
def getRpmtMetrics(readCounts):

   # handle zero MTs case
   if len(readCounts) == 0:
      outvec = tuple([0] * 6)
      return outvec
   
   # main metrics
   totalReads = sum(readCounts)
   totalMts = len(readCounts)
   rpmtMean = round(1.00 * sum(readCounts) / len(readCounts), 1)
   
   # sort read counts
   readCounts.sort()
   
   # get percentiles
   pctiles = []
   for pct in (25,50,75):
      k = pct * 0.01 * (len(readCounts)-1)
      f = math.floor(k)
      c = math.ceil(k)
      if f == c:
         pctile = readCounts[int(k)]
      else:
         d0 = readCounts[int(f)] * (c - k)
         d1 = readCounts[int(c)] * (k - f)
         pctile = d0 + d1
      pctiles.append(int(round(pctile)))
      
   # done
   outvec = [totalMts, totalReads, rpmtMean]
   outvec.extend(pctiles)
   outvec = tuple(outvec)
   return outvec


#----------------------------------------------------------------------------
# main function
def run(cfg):
   print("sum_primer_umis starting...")
   readSet = cfg.readSet

   # read primers, init data structure
   primers = {}
   for line in open(cfg.primerFile, "r"):
      (chrom, loc3, direction, primer) = line.strip().split("\t")
      strand = 0 if direction == "L" or direction == "0" else 1
      loc3 = int(loc3)
      loc5 = loc3 - len(primer) + 1 if strand == 0 else loc3 + len(primer) - 1
      primers[primer] = (strand, chrom, loc5, loc3, [])

   # get MTs and supporting read counts from disk
   for line in open(readSet + ".umi_mark.for.sum.primer.txt", "r"):
      (chrom, strand, umiLoc, umi, numReads, numAlignments, mtReadIdx, isResample, fragLen, primer, primerLoc5) = line.strip().split("|")
      primers[primer][-1].append(int(numReads))

   # define metric names
   metricNames = ("UMIs", "read fragments"
   ,"read fragments per UMI, mean"
   ,"read fragments per UMI, 25th percentile"
   ,"read fragments per UMI, 50th percentile"
   ,"read fragments per UMI, 75th percentile")

   # write column header for main table
   fileout = open(readSet + ".sum.primer.umis.txt", "w")
   fileout.write("|".join(("read set", "primer", "strand", "chrom", "loc5", "loc3")))
   fileout.write("|")
   fileout.write("|".join(metricNames))
   fileout.write("\n")

   # for each primer, get read per MT metrics
   mtListAllPrimers = []
   for (primer, primerVec) in primers.iteritems():
   
      # unpack primer vec
      (strand, chrom, loc5, loc3, mtList) = primerVec
      
      # get read-per-barcode metrics
      metricVals = getRpmtMetrics(mtList)
      
      # output
      outrow = [readSet, primer]
      outrow.extend(primerVec[0:-1])  # dropping the MT list now
      outrow.extend(metricVals)
      outrow = (str(x) for x in outrow)
      fileout.write("|".join(outrow))
      fileout.write("\n")
      
      # save all-primers MT list
      mtListAllPrimers.extend(mtList)
      
   fileout.close()
