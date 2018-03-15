import math

# constants
FRAG_BINS = 46

#---------------------------------------------------------------------   
# metrics function
#---------------------------------------------------------------------   
def getMetrics(vals):

   # handle zero case
   if len(vals) == 0:
      outvec = tuple([0] * 6)
      return outvec
   
   # main metrics
   valSum  = sum(vals)
   valNum  = len(vals)
   valMean = round(1.00 * valSum / valNum, 1)
   
   # start output vector
   outvec = [valNum, valSum, valMean]
   
   # sort values
   vals.sort()
   
   # get percentiles (quartiles here)
   for pct in (25,50,75):
      k = pct * 0.01 * (len(vals)-1)
      f = math.floor(k)
      c = math.ceil(k)
      if f == c:
         pctile = vals[int(k)]
      else:
         d0 = vals[int(f)] * (c - k)
         d1 = vals[int(c)] * (k - f)
         pctile = d0 + d1
      outvec.append(int(round(pctile)))
      
   # done
   return tuple(outvec)

#---------------------------------------------------------------------   
# main function
#---------------------------------------------------------------------   
def run(cfg):
   print("umi_frags starting...")
   # params
   readSet = cfg.readSet

   # init counters and vector for molecule metrics - this can get long - as many as input molecule count!
   moleculeData = []
   totMolecules = 0
   totMoleculesWithResampleRead = 0
   totReads = 0
   totReadsFromResample = 0
   
   # read unique molecule information from "umi" module
   for line in open(readSet + ".umi_mark.alignments.txt", "r"):
   
      # parse line
      (pChrom, pStrand, mtLoc, mt, mtReads, mtReads_, mtReadIdx, isResample, fragLen, pLoc5, primer, umiSeq, alignLocR, alignLocP, readId, read1L, read1R, cigar1, read2L, read2R, cigar2) = line.strip().split("|")
      mtReadIdx = int(mtReadIdx)
      numReads = int(mtReads_)
      isResampleRead = int(isResample) > 0
      
      # read accounting
      totReads += 1
      if isResampleRead:
         totReadsFromResample += 1
      
      # save one value per molecule in large list
      if mtReadIdx == 0:
         moleculeData.append((int(fragLen),numReads))
      
      # count unique molecules by whether it has a resampling read
      elif mtReadIdx == numReads - 1:
         totMolecules += 1
         if isResampleRead:
            totMoleculesWithResampleRead += 1

   # sort the molecule data vector by fragment length
   moleculeData.sort()
   
   # output read per molecule distribution for each fragment length bin (bin size 20 bp)
   fragLenBinCurrent = 0
   readCountBuffer = []
   linesOut = [None] * FRAG_BINS
   for fragLen, numReads in moleculeData:
   
      # get frag len bin
      fragLenBin = min(fragLen / 20, FRAG_BINS - 1)
      
      # flush data if next bin found
      if fragLenBin > fragLenBinCurrent:
         linesOut[fragLenBinCurrent] = getMetrics(readCountBuffer)
         del(readCountBuffer[:])
         fragLenBinCurrent = fragLenBin

      # save read count for later
      readCountBuffer.append(numReads)         
   
   # process last frag len bin
   linesOut[fragLenBinCurrent] = getMetrics(readCountBuffer)
   
   # open output file for fragment length distribution, write column headers
   fileout = open(readSet + ".umi_frags.len-distrib.txt", "w")
   outvec = ["read set", "frag len", "UMIs", "reads", "rpumiMean", "rpumi25", "rpumi50", "rpumi75"]
   fileout.write("|".join(outvec))
   fileout.write("\n")
   
   # write metrics by frag len
   for idxFragLen in range(0,FRAG_BINS):
   
      # first columm is frag len
      outvec = [readSet, idxFragLen * 20 + 10]
      
      # add metrics
      metrics = linesOut[idxFragLen]
      if metrics == None:
         metrics = [0 for x in range(6)]
      outvec.extend(metrics)
      
      # write to disk
      fileout.write("|".join([str(x) for x in outvec]))
      fileout.write("\n")
      
   # done
   fileout.close()
   
   # get overall metrics for whole read set
   metricsFragLens = getMetrics([x[0] for x in moleculeData])
   metricsReadNums = getMetrics([x[1] for x in moleculeData])

   # set up metric names and values for output
   metricNames = ("UMIs", "read fragments"
   ,"read fragments per UMI, mean"
   ,"read fragments per UMI, 25th percentile"
   ,"read fragments per UMI, 50th percentile"
   ,"read fragments per UMI, 75th percentile"
   ,"sample fragment length, 25th percentile (bp)"
   ,"sample fragment length, 50th percentile (bp)"
   ,"sample fragment length, 75th percentile (bp)"
   )
   metricVals = list(metricsReadNums)
   metricVals.extend(metricsFragLens[3:])

   # write metrics to summary file
   fileout = open(readSet + ".umi_frags.summary.txt", "w")
   for metricVal, metricName in zip(metricVals,metricNames):
      fileout.write("{}\t{}\n".format(metricVal,metricName))
   
   # write internal priming metrics
   fileout.write("{:4.1f}\t% of UMIs with >= one read from internal downstream priming\n".format(100.0 * totMoleculesWithResampleRead / totMolecules))
   fileout.write("{:4.1f}\t% of reads from internal downstream priming\n".format(100.0 * totReadsFromResample / totReads))
   fileout.close()
   
   # also save mean reads-per-molecule, for possible use in variant calling
   cfg.readsPerUmi = metricsReadNums[2]
   
#----------------------------------------------------------------------------------------------
# run from the command line
#----------------------------------------------------------------------------------------------
if __name__ == "__main__":
   cfg = lambda:0
   cfg.readSet = "SRR3493403"
   run(cfg)
