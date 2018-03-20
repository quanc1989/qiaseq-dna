#---------------------------------------------------------------------------------------
# merge regions (3 or 4 col output)
#---------------------------------------------------------------------------------------
def merge(regionsIn):

   # done if input is an empty list
   if len(regionsIn) == 0:
      return regionsIn
      
   # sort regions (inefficient because already sorted sometimes)
   regionsIn.sort()

   # allocate output list
   regionsOut = []
   
   # do a 3-col merge
   if len(regionsIn[0]) == 3:
      regionsOut.append(regionsIn[0])
      for region2 in regionsIn[1:]:
         (chrom2, locL2, locR2) = region2
         (chrom1, locL1, locR1) = regionsOut[-1]
         if chrom2 != chrom1 or locL2 > locR1:
            regionsOut.append((chrom2, locL2, locR2))
         elif locR2 > locR1:
            regionsOut[-1] = (chrom1, locL1, locR2)

   # do a 4-col merge
   else:
      regionsOut.append(regionsIn[0][0:4])
      for region2 in regionsIn[1:]:
         (chrom2, locL2, locR2, anno2) = region2[0:4]
         (chrom1, locL1, locR1, anno1) = regionsOut[-1]
         if chrom2 != chrom1 or locL2 > locR1:
            regionsOut.append((chrom2, locL2, locR2, anno2))
         elif anno1 == anno2:
            if locR2 > locR1:
               regionsOut[-1] = (chrom1, locL1, locR2, anno1)
         else:
            anno = set()
            anno.update(anno1.split(";"))
            anno.update(anno2.split(";"))
            anno = list(anno)
            anno.sort()
            anno = ";".join(anno)
            locR = max(locR2,locR1)
            regionsOut[-1] = (chrom1, locL1, locR, anno)

   # done        
   return regionsOut

#---------------------------------------------------------------------------------------
# subtract B regions from A region (returns 3-col bed)
#---------------------------------------------------------------------------------------
def subtract(a,b):

   # first remove unused column 4
   if len(a) > 0 and len(a[0]) > 3:
      a = [x[0:3] for x in a]
   if len(b) > 0 and len(b[0]) > 3:
      b = [x[0:3] for x in b]

   # A and B must both be merged
   a = merge(a)
   b = merge(b)

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
   c = merge(c)
   return c
   
#---------------------------------------------------------------------------------------------
# intersect A and B regions, return B elements that are in A region (all B cols are preserved)
#---------------------------------------------------------------------------------------------
def intersect(a,b):

   # a must be merged!!!
   a = merge(a)

   # save region boundaries
   vec = []
   for aVec in a:
      (chrom, locL, locR) = aVec[0:3]
      vec.append((chrom, locL,  1, -1))
      vec.append((chrom, locR, -1, -1))
   bIdx = 0
   for bVec in b:
      (chrom, locL, locR) = bVec[0:3]
      vec.append((chrom, locL,  1, bIdx))
      vec.append((chrom, locR, -1, bIdx))
      bIdx += 1
   vec.sort()

   # find B regions in A region (at least 1 bp intersection)
   bOut = set()
   bLocal = set()
   depthA = 0
   for (chrom, loc, depth, bIdx) in vec:
      if bIdx == -1:
         depthA += depth
         bOut.update(bLocal)
      else:
         if depth == 1:
            bLocal.add(bIdx)
         else:
            bLocal.remove(bIdx)
         if depthA == 1:
            bOut.update(bLocal)
   
   # return B regions
   out = []
   for bIdx in bOut:
      out.append(b[bIdx])
   out.sort()
   return out

#---------------------------------------------------------------------------------------
# annotate A region with column 4 from overlapping B regions (4-col output)
#   B regions must be merged, or at least have unique col4
#---------------------------------------------------------------------------------------
def annotate(a,b):

   # save region boundaries
   vec = []
   aIdx = 0
   for aVec in a:
      (chrom, locL, locR) = aVec[0:3]
      vec.append((chrom, locL,  1, 1, aIdx))
      vec.append((chrom, locR, -1, 1, aIdx))
      aIdx += 1
   for bVec in b:
      (chrom, locL, locR, col4) = bVec[0:4]
      vec.append((chrom, locL,  1, 0, col4))
      vec.append((chrom, locR, -1, 0, col4))
   vec.sort()

   # allocate A annotation sets
   aAnno = [set() for i in range(0,len(a))]
   
   # annotate A region with column 4 from overlapping B regions
   aLocal = set()
   bLocal = set()
   for (chrom, loc, direction, type, payload) in vec:
      # process B record (source annotation)
      if type == 0:
         if direction == -1:
            bLocal.remove(payload)
         else:
            bLocal.add(payload)
            for aIdx in aLocal:
               aAnno[aIdx].add(payload)
      # process A record (destination annotation)
      else:
         if direction == -1:
            aLocal.remove(payload)
         else:
            aLocal.add(payload)
            for col4 in bLocal:
               aAnno[payload].add(col4)

   # concatenate the COL4 annotation
   aIdx = 0
   aOut = []
   for aVec in a:
      (chrom, locL, locR) = aVec[0:3]
      anno = list(aAnno[aIdx])
      anno.sort()
      anno = [str(x) for x in anno]
      if len(anno) > 0 and len(anno) < 10:
         anno = ";".join(anno)
      else:
         anno = chrom + ":" + str(round(locL / 1000000.0,1)) + "Mb"
      aOut.append((chrom, locL, locR, anno))
      aIdx += 1

   # done
   return aOut

#---------------------------------------------------------------------------------------
# trim regions
#---------------------------------------------------------------------------------------
def trim(regionsIn, bpToTrim):

   # trim some bases off input region
   regionsOut = []
   for row in regionsIn:
      (chrom, locL, locR) = row[0:3]
      locL += bpToTrim
      locR -= bpToTrim
      if locR > locL:
         out = [chrom, locL, locR]
         out.extend(row[3:])
         out = tuple(out)
         regionsOut.append(out)
   
   # if trim was negative (adding bases), then need to merge new overlaps
   if bpToTrim < 0:
      regionsOut = merge(regionsOut)

   # done
   return regionsOut
   
#----------------------------------------------------------------------------------------
# write bed to disk (actually does not need to be bed)
#----------------------------------------------------------------------------------------
def write(regions, fileName):
   fileout = open(fileName, "w")
   fileout.write("track name='{0}' description='{0}'".format(fileName.replace(".bed","")))
   fileout.write("\n")
   for row in regions:
      row = (str(x) for x in row)
      fileout.write("\t".join(row))
      fileout.write("\n")
   fileout.close()

#------------------------------------------------------------------
# get bedgraph of depth (e.g. amplicons)
#------------------------------------------------------------------
def getCoverageDepth(bedIn):
   # prep input
   vec = []
   for row in bedIn:
      (chrom, locL, locR) = row[0:3]
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
      
   # done
   return bedgraph

   