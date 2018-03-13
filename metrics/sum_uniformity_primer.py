#------------------------------------------------------------------------------------------
def getUniformityMetrics(depths):

   # get total cummulative depth
   depthTotal = sum(depths)
   
   # if zero depth at all primers/sites, return zeros
   if depthTotal == 0:
      outvec = [0] * 5
      return tuple(outvec)
   
   # sort
   depths.sort()
   
   # get mean depth
   depthMean = 1.00 * depthTotal / len(depths)
      
   # get % of mean metrics
   pctGtVec = []
   for pct in (5.0, 10.0, 20.0, 30.0):
      for idx in range(len(depths)):
         pctMean = 100.00 * depths[idx] / depthMean
         if pctMean >= pct:
            x = 100.00 - (100.00 * idx) / len(depths)
            pctGtVec.append(x)
            break

   # done
   outvec = [depthMean]
   outvec.extend(pctGtVec)
   return tuple(outvec)
   
#------------------------------------------------------------------------------------------
def run(cfg):
   print("sum_uniformity_primer starting...")
   readSet = cfg.readSet

   # initialize a depth vector for each depth type
   depthVecs = ([], [])  # MTs, reads
      
   # loop over primer read depths from previous step, compute raw and MT-corrected uniformity
   for line in open(readSet + ".sum.primer.umis.txt", "r"):
   
      # skip column header line
      if line.startswith("read set"):
         continue
      
      # parse line
      vals = line.strip().split("|")
      (readSet, primer, strand, chrom, loc5, loc3, numMts, numReads) = vals[0:8]
      
      # save depth count for later
      depthVecs[0].append(int(numMts))
      depthVecs[1].append(int(numReads))
      
   # open output file
   fileout = open(readSet + ".sum.uniformity.primer.summary.txt", "w")
   
   # write num of primers
   fileout.write("{}\t# of primers\n".format(len(depthVecs[0])))

   # define metric labels and get metric values for each depth type
   metricNames = ("mean primer", "% of primers >= 5% of mean", "% of primers >= 10% of mean", "% of primers >= 20% of mean", "% of primers >= 30% of mean")
   depthTypes = ("UMI", "read fragment")
   
   # write primer uniformity metrics to output file
   for idxDepthType in range(len(depthTypes)):
      depthType = depthTypes[idxDepthType]
      metricVals = getUniformityMetrics(depthVecs[idxDepthType])
      for idxMetric in range(len(metricNames)):
         metricName = metricNames[idxMetric]
         metricVal  = metricVals[idxMetric]
         fileout.write("{:.2f}\t{} {} depth\n".format(metricVal,metricName,depthType))
   
   # done
   fileout.close() 
