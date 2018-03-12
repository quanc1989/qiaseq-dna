def run(cfg):
   print("sum_specificity starting...")
   readSet = cfg.readSet

   # open output
   fileout = open(readSet + ".sum.specificity.txt", "w")

   # write column headers
   fileout.write("|".join(("read set","chrom","loc5","loc3","strand","primer","on target","off target","offT 1","offT 2","offT 3","offT 4","offT 5","offT 6","offT 7","offT 8","offT 9","offT 10","offT 1","offT 2","offT 3","offT 4","offT 5","offT 6","offT 7","offT 8","offT 9","offT 10")))
   fileout.write("\n")
   
   # get read depth at all primer sites (WARNING: this is memory unbounded! fix later.)
   siteDepthsOnT  = {}
   siteDepthsOffT = {}
   for line in open(readSet + ".umi_filter.alignments.txt", "r"):
      (pChrom, pLoc5, pStrand, primer, barcode, isIntendedSite, alignChrom, alignStrand, alignLoc1, alignLoc2, readId, read1L, read1R, cigar1, read2L, read2R, cigar2) = line.strip().split("|")

      (pLoc5, pStrand, alignLoc, alignStrand, isIntendedSite) = (int(x) for x in (pLoc5, pStrand, alignLoc2, alignStrand, isIntendedSite))
      
      # on target sites
      if isIntendedSite == 1:
         key = (primer, pChrom, pStrand, pLoc5)
         if key in siteDepthsOnT:
            siteDepthsOnT[key] += 1
         else:
            siteDepthsOnT[key]  = 1
            
      # off-target sites
      else:
         key = (primer, alignChrom, alignStrand, alignLoc)
         if key in siteDepthsOffT:
            siteDepthsOffT[key] += 1
         else:
            siteDepthsOffT[key]  = 1
            
   # re-organize off-target site depths by primer
   stacksOffT = {}
   for (key, numReads) in siteDepthsOffT.iteritems():
      (primer, alignChrom, alignStrand, alignLoc) = key
      if primer in stacksOffT:
         vec = stacksOffT[primer]
      else:
         vec = []
         stacksOffT[primer] = vec
      vec.append((numReads, alignChrom, alignLoc, alignStrand))
         
   # descending read depth sort off-target sites for each primer; also extend vector to minimum of 10
   for primer in stacksOffT.iterkeys():
      vec = stacksOffT[primer]
      vec.sort(reverse=True)
      if len(vec) <= 10:
         for i in range(len(vec), 10):
            vec.append((0, "_", "_", "_"))
      else:
         numReads10up = sum([x[0] for x in vec[9:]])
         vec[9] = (numReads10up, "_", "_", "_")
         stacksOffT[primer] = vec[0:10]

   # read primers, output read depths
   for line in open(cfg.primerFile, "r"):
      (chrom, loc3, direction, primer) = line.strip().split("\t")
      strand = 0 if direction == "L" or direction == "0" else 1
      loc3 = int(loc3)
      loc5 = loc3 - len(primer) + 1 if strand == 0 else loc3 + len(primer) - 1
      
      # start outvec
      outvec = [readSet, chrom, loc5, loc3, strand, primer]
      
      # add on-target read count
      key = (primer, chrom, strand, loc5)
      readsOnT = 0 if key not in siteDepthsOnT else siteDepthsOnT[key]
      outvec.append(readsOnT)
      
      # add off-target top 10 site read count
      readsOffT = [0] * 10 if primer not in stacksOffT else [x[0] for x in stacksOffT[primer]]
      outvec.append(sum(readsOffT))
      outvec.extend(readsOffT)
      
      # add off-target top 10 site locations
      sitesOffT = ["_"] * 10 if primer not in stacksOffT else [" ".join([str(y) for y in x[1:]]) for x in stacksOffT[primer]]
      outvec.extend(sitesOffT)
      
      # write to disk
      fileout.write("|".join((str(x) for x in outvec)))
      fileout.write("\n")

   # done
   fileout.close()      
