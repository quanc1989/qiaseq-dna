import os
import os.path
import subprocess

#------------------------------------------------------------------------
# bed merge
#------------------------------------------------------------------------
def bedMerge(bedIn):

   # done if input is an empty list
   if len(bedIn) == 0:
      return bedIn
      
   # sort regions (inefficient because already sorted sometimes)
   bedIn.sort()

   # do a 3-column merge
   bedOut = []
   bedOut.append(bedIn[0])
   for region2 in bedIn[1:]:
      (chrom2, locL2, locR2) = region2
      (chrom1, locL1, locR1) = bedOut[-1]
      if chrom2 != chrom1 or locL2 > locR1:
         bedOut.append((chrom2, locL2, locR2))
      elif locR2 > locR1:
         bedOut[-1] = (chrom1, locL1, locR2)
    
   # done
   return bedOut

#------------------------------------------------------------------------
# bed subtract
#------------------------------------------------------------------------
def bedSubtract(a,b):

   # A and B must both be merged
   a = bedMerge(a)
   b = bedMerge(b)

   # save region boundaries
   vec = []
   for aVec in a:
      (chrom, locL, locR) = aVec
      vec.append((chrom, locL,  1))
      vec.append((chrom, locR, -1))
   for bVec in b:
      (chrom, locL, locR) = bVec
      vec.append((chrom, locL, -1))
      vec.append((chrom, locR,  1))
   vec.sort()

   # subtract C = A - B
   c = []
   depthNet = 0
   locLast = None
   for (chrom, loc, depth) in vec:
      if depthNet == 1 and loc - locLast > 0:
         c.append((chrom, locLast, loc))
      locLast = loc
      depthNet += depth

   # merge bookends
   c = bedMerge(c)
   return c
   
#------------------------------------------------------------------------
# get target region bed   
#------------------------------------------------------------------------
def getTargetBed(cfg):

   # make target region
   bedTarget = []

   # read ROI from disk (assuming this is already merged and sorted!)
   if "roiBedFile" in cfg.__dict__:
      for line in open(cfg.roiBedFile, "r"):
         if line.startswith("track "):
            continue
         vals = line.strip().split("\t")
         chrom, locL, locR = vals[0:3]
         bedTarget.append((chrom, int(locL), int(locR)))
  
   # if no ROI file, use primer regions
   else:
      targetWindowSize = 150
      for line in open(cfg.primerFile, "r"):
         (chrom, loc3, direction, primer) = line.strip().split("\t")
         strand = 0 if direction == "L" or direction == "0" else 1
         loc3 = int(loc3)
         if strand == 0:
            loc5 = loc3 - len(primer) + 1
            bedRow = (chrom, loc3+1,loc5+targetWindowSize)
         else:
            loc5 = loc3 + len(primer) - 1
            bedRow = (chrom, loc5-targetWindowSize+1,loc3)
         bedTarget.append(bedRow)

      # merge ROI region (also causes a sort)
      bedTarget = bedMerge(bedTarget)

      # save ROI to disk, for possible use downstream
      trackName = cfg.readSet + ".roi"
      fileName = trackName + ".bed"
      cfg.roiBedFile = fileName
      fileout = open(fileName, "w")
      fileout.write("track name='{0}' description='{0}'\n".format(trackName))
      for row in bedTarget:
         fileout.write("\t".join((str(x) for x in row)))
         fileout.write("\n")
      fileout.close()

   bpTarget = sum((x[2] - x[1] for x in bedTarget))   
   # done
   return (bedTarget,bpTarget)

#------------------------------------------------------------------------
# get UMI depth bedgraph, in target region only
#------------------------------------------------------------------------
def getDepthsInRoi(bedTarget,depthFile):
   
   # read depth bedgraph file
   bedgraphDepthsAll = []
   for line in open(depthFile, "r"):
      if line.startswith("track "):
         continue
      (chrom, locL, locR, depth) = line.strip().split("\t")
      bedgraphDepthsAll.append((chrom, int(locL), int(locR), int(depth)))

   # add regions in ROI, but outside of bedgraph to the depth bedgraph
   bedTargetZeroDepth = bedSubtract(bedTarget, [x[0:3] for x in bedgraphDepthsAll])
   for (chrom, locL, locR) in bedTargetZeroDepth:
      bedgraphDepthsAll.append((chrom,locL,locR,0))
   bedgraphDepthsAll.sort()

   # make new bedgraph for only the target region
   bedgraphDepths = []
   idx = 0
   for (chrom, locL, locR) in bedTarget:

      # spin bedgraph up to this loc
      while True:
         (chrom_, locL_, locR_, depth) = bedgraphDepthsAll[idx]

         # region entirely before ROI, keep going
         if chrom_ != chrom or locR_ <= locL:
            pass
         
         # region entirely within the ROI
         elif locL <= locL_ and locR_ <= locR:
            bedgraphDepths.append((chrom, locL_, locR_, depth))
   
         # region entirely spanning the ROI
         elif locL_ <= locL and locR <= locR_:
            bedgraphDepths.append((chrom, locL , locR , depth))

         # region spanning the ROI left edge
         elif locL_ <= locL and locR_ <= locR:
            bedgraphDepths.append((chrom, locL , locR_, depth))
     
         # region spanning the ROI right edge
         elif locL <= locL_ and locR <= locR_:
            bedgraphDepths.append((chrom, locL_, locR , depth))

         # region beyond the ROI
         else:
            raise Exception("error generating bedgraph within ROI")

         # get next depth block
         if locR_ < locR or chrom_ != chrom:
            idx+= 1
            
         # set up next ROI region - spin bedgraph backward
         else:
            while locR_ >= locR and chrom_ == chrom and idx > 0:
               idx -= 1
               (chrom_, locL_, locR_, umiDepth) = bedgraphDepthsAll[idx]
            break

   # debug check
   bpTarget = sum((x[2]-x[1] for x in bedTarget))
   bpDepths = sum((x[2]-x[1] for x in bedgraphDepths))
   if bpDepths != bpTarget:
      raise Exception("depth bedgraph does not match target region!")

   # done
   return bedgraphDepths

#------------------------------------------------------------------------------------------
# enrichment depth uniformity metrics
#------------------------------------------------------------------------------------------
def getUniformityMetrics(bedgraphDepths,fileout,metricType):

   # make depth vector - as big as the target region
   depths = []
   for chrom, locL, locR, depth in bedgraphDepths:
      for loc in range(locL,locR):
         depths.append(depth)

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
   outvec = [depthMean]
   for pct in (5.0, 10.0, 20.0, 30.0):
      for idx in range(len(depths)):
         pctMean = 100.00 * depths[idx] / depthMean
         if pctMean >= pct:
            x = 100.00 - (100.00 * idx) / len(depths)
            outvec.append(x)
            break

   # write to summary file
   metricNames = ("mean"
   , "% of bases >= 5% of mean"
   , "% of bases >= 10% of mean"
   , "% of bases >= 20% of mean"
   , "% of bases >= 30% of mean")
   for metricVal, metricName in zip(outvec, metricNames):
      fileout.write("{:.2f}\t{} {} depth\n".format(metricVal,metricName,metricType))
   
#------------------------------------------------------------------------
# get LOD estimates
#------------------------------------------------------------------------
def getLodEstimates(cfg,fileoutSummary,bedgraphDepths):
   # get read set
   readSet = cfg.readSet

   # write bedgraph to disk for use by R script, without track name
   fileout = open(readSet + ".umi_depths.lod.txt", "w")
   for row in bedgraphDepths:
      fileout.write("|".join((str(x) for x in row)))
      fileout.write("\n")
   fileout.close()

   # compute mean MT depth
   mtSum = 0
   bpSum = 0
   for (chrom, locL, locR, mtDepth) in bedgraphDepths:
      bp = locR - locL
      bpSum += bp
      mtSum += bp * mtDepth
   umiDepthMean = int(round(1.00 * mtSum / bpSum))
   cfg.umiDepthMean = umiDepthMean
   print("umi_depths: mean UMI depth over target region:", umiDepthMean)

   # get path to this script - R script in same dir
   scriptName = os.path.abspath(__file__)
   (scriptName, ext) = os.path.splitext(scriptName)
   scriptName = scriptName + "_lod.R"

   # call R script for LOD estimate
   fileInName  = readSet + ".umi_depths.lod.txt"
   fileOutName = readSet + ".umi_depths.lod.temp.txt" 
   cmd = "Rscript {} {} {} {}".format(scriptName, umiDepthMean, fileInName, fileOutName)
   subprocess.check_call(cmd,shell=True)

   # add the bedgraph track line to the R output file
   fileout = open(readSet + ".umi_depths.lod.bedgraph", "w")
   fileout.write("track type=bedGraph name='" + readSet + ".umi_depths.lod'\n")
   for line in open(readSet + ".umi_depths.lod.temp.txt", "r"):
      fileout.write(line)
   fileout.close()
   
   # format the LOD percentiles and write to summary file
   for line in open(readSet + ".umi_depths.lod.temp.txt.quantiles.txt","r"):
      (metricName, metricVal) = line.strip().split("|")
      metricName = int(metricName.replace("%",""))
      metricVal = float(metricVal)
      thorst = "st" if metricName == 1 else "th"
      fileoutSummary.write("{:6.4f}\t{:2d}{} percentile estimated minimum detectible allele fraction (LOD)\n".format(metricVal, metricName,thorst))

   # remove the temporary files
   os.remove(readSet + ".umi_depths.lod.txt")
   os.remove(readSet + ".umi_depths.lod.temp.txt")
   os.remove(readSet + ".umi_depths.lod.temp.txt.quantiles.txt")

#------------------------------------------------------------------------
# output one locus for UMI depth bedgraph   
#------------------------------------------------------------------------
def handleOneLocus(bedIn,fileout):

   # done if nothing to do yet
   if len(bedIn) == 0:
      return

   # prep input
   vec = []
   for row in bedIn:
      (chrom, locL, locR) = row
      vec.append((chrom, locL,  1))
      vec.append((chrom, locR, -1))
   vec.sort()
   
   # count coverage depth
   depthVec = []
   depthNet = 0
   for (chrom, loc, depth) in vec:
      depthNet += depth
      depthVec.append((chrom, loc, depthNet))
      
   # convert to bedgraph format
   bedgraph = []
   (chrom_, loc_, depth_) = depthVec[0]
   for (chrom, loc, depth) in depthVec[1:]:
      if loc != loc_ or chrom != chrom_:
         bedgraph.append((chrom, loc_, loc, depth_))
      chrom_ = chrom
      loc_   = loc
      depth_ = depth
      
   # output
   for row in bedgraph:
      fileout.write("\t".join((str(x) for x in row)))
      fileout.write("\n")
   
#------------------------------------------------------------------------
# make a UMI depth bedgraph   
#------------------------------------------------------------------------
def makeUmiDepthBedgraph(cfg,readSupportMin,trackName):
   # params
   readSet = cfg.readSet
   
   # open output file
   trackName = readSet + "." + trackName
   fileout = open(trackName + ".bedgraph", "w")
   fileout.write("track type=bedGraph name='" + trackName + "'\n")
   
   # init genome locus buffer
   locusChrom = "foobar"
   locusLocR = -2001
   locusBedFrags = []
   
   # read unique molecule information from "umi" module (file must be sorted by random fragmentation position)
   for line in open(readSet + ".umi_mark.alignments.txt", "r"):
   
      # parse line
      (pChrom, pStrand, mtLoc, mt, mtReads, mtReads_, mtReadIdx, isResample, fragLen, pLoc5, primer, umiSeq, alignLocR, alignLocP, readId, read1L, read1R, cigar1, read2L, read2R, cigar2) = line.strip().split("|")
      mtReadIdx = int(mtReadIdx)
      
      # skip PCR replicates, unless requested
      if mtReadIdx > 0 and readSupportMin > 0:
         continue
         
      # filter molecules with insufficient read support
      if int(mtReads_) < readSupportMin:
         continue

      # covert to int
      pLoc5, read1L, read1R, read2L, read2R  = map(int,(pLoc5, read1L, read1R, read2L, read2R))
      
      # take off the primer region
      if int(pStrand) == 1:
         pLoc3 = pLoc5 - len(primer) + 1
         read2R = min(read2R, pLoc3)
         read1R = min(read1R, pLoc3)
         locs = (read2L,read2R,read1L,read1R)
      else:
         pLoc3 = pLoc5 + len(primer)
         read2L = max(read2L, pLoc3)
         read1L = max(read1L, pLoc3)
         locs = (read1L,read1R,read2L,read2R)

      # merge R1 and R2
      loc0,loc1,loc2,loc3 = locs
      if loc2 > loc1 or readSupportMin == 0:
         bedFrag = [(pChrom, loc0, loc1), (pChrom, loc2, loc3)]
      else:
         bedFrag = [(pChrom, loc0, loc3)]

      # flush previous data to disk, if new locus island found (note: lot of memory needed for long contigous single locus panel!)
      if pChrom != locusChrom or loc0 > locusLocR + 2000:
         handleOneLocus(locusBedFrags,fileout)
         del(locusBedFrags[:])
         locusChrom = pChrom
         locusLocR = loc3

      # save for later
      locusBedFrags.extend(bedFrag)
      locusLocR = max(locusLocR,loc3)

   # process final locus
   handleOneLocus(locusBedFrags,fileout)
   fileout.close()

#---------------------------------------------------------------------   
# main function
#---------------------------------------------------------------------   
def run(cfg):
   print("umi_depths starting...")
   readSet  = cfg.readSet
   numCores = cfg.numCores
   
   # sort read alignments by R2 random fragmentation location (to remove strand sort from previous step)
   fileNameIn = readSet + ".umi_mark.alignments.txt"
   cmd = "sort -k1,1 -k3,3n -t\| -T./ --parallel={1} {0} > {0}.temp.txt".format(fileNameIn, numCores)  # sort order is (alignChrom, mtLoc)
   subprocess.check_call(cmd, shell=True)
   os.rename(fileNameIn + ".temp.txt", fileNameIn)
   print("umi_depths: done sorting UMI read alignment file by random fragmentation position")
   
   # make depth bedgraph files, save to disk
   makeUmiDepthBedgraph(cfg,0, "umi_depths.raw_reads")  # raw read depth
   makeUmiDepthBedgraph(cfg,1, "umi_depths")            # UMI depth
   makeUmiDepthBedgraph(cfg,4, "umi_depths.ge4reads")   # UMI with >=4 supporting reads depth
   print("umi_depths: done making depth bedgraphs")
   
   # get target region from disk, or make it if not specified
   bedTarget,bpTarget = getTargetBed(cfg)
   
   # open summary file
   fileout = open(readSet + ".umi_depths.summary.txt", "w")
   fileout.write("{}\t# of target bases\n".format(bpTarget))
   # remove non-target regions from UMI depth bedgraph, compute uniformity metrics
   bedgraphDepths = getDepthsInRoi(bedTarget, readSet + ".umi_depths.bedgraph")
   getUniformityMetrics(bedgraphDepths, fileout, "UMI")
   print("umi_depths: done making UMI depth bedgraph over ROI only, and computing depth uniformity")
   
   # remove non-target regions from raw read depth bedgraph, compute uniformity metrics
   bedgraphDepthsRaw = getDepthsInRoi(bedTarget, readSet + ".umi_depths.raw_reads.bedgraph")
   getUniformityMetrics(bedgraphDepthsRaw, fileout, "read")
   os.remove(readSet + ".umi_depths.raw_reads.bedgraph")
   print("umi_depths: done making UMI depth bedgraph over ROI only, and computing depth uniformity")
   
   # get allele fraction limit-of-detection (LOD) estimate for each locus
   getLodEstimates(cfg,fileout,bedgraphDepths)
   print("umi_depths: done making computing allele fraction LOD")
   
   # close summary file
   fileout.close()
   
   # write bases of ROI that have MT depth below 20% of mean depth - can be concatenated
   umiDepthLow = 0.20 * cfg.umiDepthMean
   fileout = open(readSet + ".umi_depths.LT20PctOfMean.txt", "w")
   fileout.write("\t".join(("read set", "chrom", "loc", "UMI depth")))
   fileout.write("\n")
   for (chrom, locL, locR, umiDepth) in bedgraphDepths:
      if umiDepth < umiDepthLow:
         for loc in range(locL, locR):
            outvec = (readSet, chrom, loc, umiDepth)
            fileout.write("\t".join((str(x) for x in outvec)))
            fileout.write("\n")
   fileout.close()
   print("umi_depths: done writing regions below 20% of mean UMI depth")
   
#----------------------------------------------------------------------------------------------
# run from the command line
#----------------------------------------------------------------------------------------------
if __name__ == "__main__":
   cfg = lambda:0
   cfg.readSet = "NEB_S2"
   cfg.numCores = 64
   cfg.primerFile =  "/srv/qgen/example/DHS-101Z.primers.txt"
   cfg.roiBedFile =  "/srv/qgen/example/DHS-101Z.roi.bed"
   run(cfg)
