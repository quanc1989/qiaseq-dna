import math

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
   for line in open(readSet + ".umi_mark.alignments.txt", "r"):
      (chrom, strand, mtLoc, mt, numReads, numAlignments, mtReadIdx, isResample, fragLen, primerLoc5, primer, umiSeq, alignLocR, alignLocP, readId, read1L, read1R, cigar1, read2L, read2R, cigar2) = line.strip().split("|")
      (strand, chrom, loc5, loc3, mtList) = primers[primer]
      mtList.append(int(numReads))

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
