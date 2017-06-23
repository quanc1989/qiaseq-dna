import editdist
import sys
MAX_BC_LEN = 50
mtx = [[0] *(MAX_BC_LEN+1) for idx in xrange(MAX_BC_LEN+1)]

#---------------------------------------------------------------------------------------
# fast check to determine if a barcodeB can be merged with barcodeA
#---------------------------------------------------------------------------------------
def isSimilar(bcA, bcB, bcLen, allowNs = True):
   aLen = len(bcA)
   bLen = len(bcB)
   maxLen = max(aLen, bLen)
   prefixMatchLen = 0
   misMatchPos = -1
   while bcA[prefixMatchLen] == bcB[prefixMatchLen] or (allowNs == True and bcB[prefixMatchLen] == 'N'):
      prefixMatchLen += 1
      if prefixMatchLen >= aLen or prefixMatchLen >= bLen:
         break
   if bLen > prefixMatchLen or aLen > prefixMatchLen:
      misMatchPos = prefixMatchLen + 1
   else:
      misMatchPos = prefixMatchLen
   if misMatchPos > bcLen:
      misMatchPos = bcLen
   else:
      misMatchPos -= 1 # covert to 0-based index
   
   bcA = bcA[::-1]
   bcB = bcB[::-1]
   aLen = aLen - prefixMatchLen
   bLen = bLen - prefixMatchLen
   suffixMatchLen = 0
   while bcA[suffixMatchLen] == bcB[suffixMatchLen] or (allowNs == True and bcB[suffixMatchLen] == 'N'):
      suffixMatchLen += 1
      if suffixMatchLen >= aLen or suffixMatchLen >= bLen:
         break
   if prefixMatchLen + suffixMatchLen >= maxLen-1:
      return (True, misMatchPos)
   else:
      return (False, misMatchPos)

#---------------------------------------------------------------------------------------
#function to cluster the barcodes for one amplicon
# inPairs: list of (readID, barcode) tuples
# bcLen: barcode length
# minRealNum: minimum number of reads with a barcode for that barcode to be considered as real in the initial scan
# minRealFrac: the barcode should a minimum of (minRealFrac * (no. of reads in the most frequent barcode) for it to be considered real in teh initial scan
# minMergeFactor: reads(barcodeA) must be >= minMergeFactor * reads(barcodeB) for barcodeB to be merged with barcodeA
#---------------------------------------------------------------------------------------
def cluster(inPairs, bcLen, minRealNum = 3, minRealFrac = 0.1, minMergeFactor = 6):
   uniqueIDs = {}
   totalCnt = 0
   prefixLen = bcLen/2
   misMatchCnt = [0] * (bcLen+1)
   distCnt = [0,0,0,0]

   # garbage collector barcode
   allNBarcode = "N" * bcLen

   # group reads by barcode
   for (readID, barcode) in inPairs:
      if barcode not in uniqueIDs:
         uniqueIDs[barcode] = set()
      uniqueIDs[barcode].add(readID)
      totalCnt += 1

	# count how many times each unique barcode occurs, and get the read count for the most frequent barcode
   largestUniq = 0
   uniqBCCnts = {}
   barcodeParent = {}
   childBarcodes = {}
   for barcode in uniqueIDs:
      uniqBCCnts[barcode] = len(uniqueIDs[barcode])
      barcodeParent[barcode] = "_UNKNOWN_" #ambiguous barcode
      childBarcodes[barcode] = []
      if uniqBCCnts[barcode] > largestUniq:
         largestUniq = uniqBCCnts[barcode]

   if allNBarcode not in uniqBCCnts:
      uniqBCCnts[allNBarcode] = 0
      barcodeParent[allNBarcode] = "_SELF_"
      childBarcodes[allNBarcode] = []

   # iteration 1: mark barcodes as real or merge them with other barcodes if they are within  1 bp of a real barcode
   prefixHash = {}
   suffixHash = {}
   sortedBarcodeList = sorted(uniqBCCnts.iteritems(), key=lambda x: x[1] , reverse=True)
   for (bcA, bcACnt) in sortedBarcodeList:
      prefix = bcA[:prefixLen]
      suffix = bcA[0-prefixLen:]
      if prefix not in prefixHash:
         prefixHash[prefix] = []
      if suffix not in suffixHash:
         suffixHash[suffix] = []
      prefixHash[prefix].append(bcA)
      suffixHash[suffix].append(bcA)
      if bcA.find("N") == -1 and uniqBCCnts[bcA] > minRealFrac * largestUniq and uniqBCCnts[bcA] >= minRealNum:
         barcodeParent[bcA] = "_SELF_"  # this is a real barcode
         continue
      for realBCList in (prefixHash[prefix], suffixHash[suffix]):
          for bcB in realBCList:
             if barcodeParent[bcB] != "_SELF_":
                continue
             (similar, misMatchPos) = isSimilar(bcB, bcA, bcLen)
             if similar:
                barcodeParent[bcA] = bcB
                childBarcodes[bcB].append(bcA)
                misMatchCnt[misMatchPos] += bcACnt
                distCnt[1] += bcACnt
                break
          if barcodeParent[bcA] != "_UNKNOWN_":
             break	# already assigned a parent in prefix list

   # iteration 2: mark barcodes as real or merge them with other barcodes if they are within 1 bp of a another real or merged barcode 
   level2Parent = {}
   level3Parent = {}
   for (bcA, bcACnt) in sortedBarcodeList:
      prefix = bcA[:prefixLen]
      suffix = bcA[0-prefixLen:]
      if barcodeParent[bcA] != "_UNKNOWN_":
         continue
      for realBCList in (prefixHash[prefix], suffixHash[suffix]):
         for bcB in realBCList:
            if barcodeParent[bcB] == "_SELF_":
               (similar, misMatchPos) = isSimilar(bcB, bcA, bcLen)
               if similar:
                  barcodeParent[bcA] = bcB
                  childBarcodes[bcB].append(bcA)
                  misMatchCnt[misMatchPos] += bcACnt
                  distCnt[1] += bcACnt
                  break
               continue
            elif barcodeParent[bcB] != "_UNKNOWN_": 
               (similar, misMatchPos) = isSimilar(bcB, bcA, bcLen)
               if similar: # checking if bcA is within 1 bp of bcB
                  if bcA not in level2Parent:
                     level2Parent[bcA] = set()
                     distCnt[2] += bcACnt
                  level2Parent[bcA].add(barcodeParent[bcB])	# do not make this parent yet: doing so will cause >2 mismatch links to parents

         # already assigned a parent: no need to check any further
         if barcodeParent[bcA] != "_UNKNOWN_" or bcA in level2Parent:	
            break

      # do a complete global alignment if we do not still find similarity 
      if barcodeParent[bcA] == "_UNKNOWN_" and bcA not in level2Parent:
         for (bcB, bcBCnt) in sortedBarcodeList:
            if barcodeParent[bcB] != "_SELF_":
               continue
            if bcBCnt < minMergeFactor * bcACnt:
               break
            editDistance = editdist.distance(bcB, bcA)
            if editDistance <= 2:
               if bcA not in level2Parent:
                  level2Parent[bcA] = set()
                  distCnt[2] += bcACnt
               level2Parent[bcA].add(bcB)
               break
            elif len(bcA) == bcLen-3:
               if bcA == bcB[3:] or bcA == bcB[0:bcLen-3]:
                  if bcA not in level3Parent:
                     level3Parent[bcA] = set()
                     distCnt[3] += bcACnt
                  level3Parent[bcA].add(bcB)
                  break
      if barcodeParent[bcA] == "_UNKNOWN_" and bcA not in level2Parent and bcA not in level3Parent:
         barcodeParent[bcA] = "_SELF_" # not within 1-bp of any child of any real barcode: this must be real as well

   # clean up and make level2 parent as the full parent
   for bcA in uniqBCCnts: 
      if bcA in level2Parent:
         bcAParent = list(level2Parent[bcA])[0]		# arbitrarily pick the first one if multiple level2 parents
         barcodeParent[bcA] = bcAParent
         childBarcodes[bcAParent].append(bcA)
         distCnt[2] += uniqBCCnts[bcA]
      elif bcA in level3Parent:
         bcAParent = list(level3Parent[bcA])[0]		# arbitrarily pick the first one if muliple level3 parents
         barcodeParent[bcA] = bcAParent
         childBarcodes[bcAParent].append(bcA)
         distCnt[3] += uniqBCCnts[bcA]

   #DEBUG
   #for (bcA, bcACnt) in sortedBarcodeList:
   #   if barcodeParent[bcA] == "_SELF_":
   #      distCnt[0] += bcACnt
   #      print ("\t".join((bcA,str(bcACnt))))
   #      for bcB in childBarcodes[bcA]:
   #         print("\t\t" + "\t".join((bcB, str(uniqBCCnts[bcB]))))
   #clusterInfo = []
   #clusterInfo.append(totalCnt)
   #clusterInfo.extend(distCnt)
   #clusterInfo.extend(misMatchCnt)
   #print("cluster stats:\t" + "\t".join((str(x) for x in clusterInfo)))
   
   # output 
   readDict = {}
   for (bcA, bcACnt) in sortedBarcodeList:
      if bcA == allNBarcode:
         continue
      if barcodeParent[bcA] == "_SELF_":
         readDict[bcA] = []
         readDict[bcA].extend(uniqueIDs[bcA])
         for bcB in childBarcodes[bcA]:
            readDict[bcA].extend(uniqueIDs[bcB])

   # reformat ouput, and hack around a bug in the algorithm above
   mts = {}
   for (mt, readIds) in readDict.iteritems():
      numReads = len(readIds)
      for readId in readIds:
         if readId in mts:  # ERROR! - read assigned to more than one MT!!!!
            #print("umi_cluster: bug in barcode clustering - read assigned to more than one MT centroid, readId: {}, MT1: {}, MT2: {}".format(readId,mts[readId],mt))
            continue
         mts[readId] = (mt,numReads)

   # done
   return mts
