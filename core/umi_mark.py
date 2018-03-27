import operator
import os
import subprocess

# our modules
import umi_cluster

# globals
WINDOW_SIZE = 6
WINDOW_OFFSET = 3

#------------------------------------------------------------------------------------------------------------------------------------
# handle reads for one putative original input molecule
#------------------------------------------------------------------------------------------------------------------------------------
def handleOneMolecule(alignments, fileout1, fileout2):

   # nothing to do yet
   if len(alignments) == 0:
      return
   
   # write number of read alignments for this molecule to all alignment records
   mtReads_ = len(alignments)

   # add read index with this molecule and output, and tag innner-priming resampling (putative)
   pLoc5Original = None
   mtReadIdx = 0

   maxFragLen = 0
   fragPrimer = None
   
   for vec in alignments:

      # unpack input
      (pChrom, pStrand, mtLoc, mt, mtReads, fragLen, pLoc5, primer, umiSeq, alignLocR, alignLocP, readId, read1L, read1R, cigar1, read2L, read2R, cigar2) = vec

     
      # mark resampling internal-priming reads
      if mtReadIdx == 0:
         pLoc5Original = pLoc5
      isResample = 1 if pLoc5 != pLoc5Original else 0
      
      # output
      vec = (pChrom, pStrand, mtLoc, mt, mtReads, mtReads_, mtReadIdx, isResample, fragLen, pLoc5, primer, umiSeq, alignLocR, alignLocP, readId, read1L, read1R, cigar1, read2L, read2R, cigar2)
      fileout1.write("|".join((str(x) for x in vec)))
      fileout1.write("\n")
      mtReadIdx += 1

      # identify primer causing longest fragment
      if fragLen > maxFragLen:
         maxFragLen = fragLen
         fragPrimer = (primer,pLoc5)
         

   # write file used by sum.primer.umis.py to disk
   (primer, primerLoc5)  = fragPrimer
   outvec = (pChrom, pStrand, mtLoc, mt, mtReads, mtReads_, mtReadIdx, isResample, maxFragLen, primer, primerLoc5)
   outvec = (str(x) for x in outvec)
   fileout2.write("|".join(outvec))
   fileout2.write("\n")
      

   # debug check - some conflicting evidence for molecule call - caused by (1) imperfection in this algorithm or (2) normal expected random collisions at high depths
   #if mtReads_ != mtReads:
   #   print("WARNING: POSSIBLE INCORRECT MOLECULE CALL", pChrom, pStrand, mtLoc, mt, mtReads, mtReads_)

#------------------------------------------------------------------------------------------------------------------------------------
# handle reads with randomly fragmented side at mapped to approximately same (+/- 3 bp) genome locus
#------------------------------------------------------------------------------------------------------------------------------------
def handleOneLocus(buffers, fileout1, fileout2):
   # debug print
   #print("********************", buffers.end1, len(buffers.mts0),len(buffers.mts1), len(buffers.buffer0), len(buffers.buffer1))
   
   # unpack buffer objects
   buffer1 = buffers.buffer1
   buffer0 = buffers.buffer0
   mts1 = buffers.mts1
   mts0 = buffers.mts0
   
   # if something in buffer1, do UMI clustering
   if len(buffer1) > 0:
   
      # buffer1: make UMI list including readId
      umiSeqs = []
      for vec in buffer1:
         (pChrom, pLoc5, pStrand, primer, umiSeq, alignLocR, alignLocP, readId, read1L, read1R, cigar1, read2L, read2R, cigar2) = vec
         umiSeqs.append((readId, umiSeq))

      # buffer1: cluster the UMI read seqs, save UMI centroid for each read
      mts1 = umi_cluster.cluster(umiSeqs, 12)
         
   # buffer0 - determine which window has the better UMI cluster
   alignmentsOut = []
   mtLoc = str(buffers.end1 - WINDOW_SIZE - WINDOW_OFFSET)
   for vec in buffer0:
      (pChrom, pLoc5, pStrand, primer, umiSeq, alignLocR, alignLocP, readId, read1L, read1R, cigar1, read2L, read2R, cigar2) = vec 
      
      # skip read if previously written out (but not physically removed from buffer0 list object, for speed)
      if readId not in mts0:
         continue
      
      # assign the UMI centroid, using best-supported (most reads) of the two overlapped UMI sequence clusterings
      (mt, mtReads) = mts0[readId]
      if readId in mts1:
         (mt_, mtReads_) = mts1[readId]
         if mtReads_ > mtReads:
            continue
         del(mts1[readId])

      # save new alignment record from buffer0
      fragLen = abs(pLoc5 - alignLocR) + 1
      vec = (pChrom, pStrand, mtLoc, mt, mtReads, fragLen, pLoc5, primer, umiSeq, alignLocR, alignLocP, readId, read1L, read1R, cigar1, read2L, read2R, cigar2)
      alignmentsOut.append(vec)
      
   # buffer0 - sort reads by UMI centroid, frag len, primer loc, umi read
   alignmentsOut.sort(key=operator.itemgetter(3,5,6,8), reverse=True)

   # buffer0 - write read alignments to disk, by UMI in order to check for "bad cluster"
   mtLast = "foobar"
   alignments = []
   for vec in alignmentsOut:
      (pChrom, pStrand, mtLoc, mt, mtReads, fragLen, pLoc5, primer, umiSeq, alignLocR, alignLocP, readId, read1L, read1R, cigar1, read2L, read2R, cigar2) = vec

      # new UMI centroid
      if mt != mtLast:
      
         # process UMI worth of reads
         handleOneMolecule(alignments, fileout1, fileout2)
         
         # prepare for next UMI
         mtLast = mt
         del(alignments[:])
         
      # save read in buffer
      alignments.append(vec)

   # process the final UMI
   handleOneMolecule(alignments, fileout1, fileout2)
      
   # buffer0: clear for next iteration
   del(buffer0[:])
   mts0.clear()

   # slide window forward
   temp = mts0
   buffers.mts0 = mts1
   buffers.mts1 = temp
   temp = buffer0
   buffers.buffer0 = buffer1
   buffers.buffer1 = temp
   buffers.end1 += WINDOW_OFFSET
   
   # copy read alignment overlapping half of the two windows reads from the old buffer1 (now buffer0) into the new buffer1 
   buffer1 = buffers.buffer1
   buffer1Start = buffers.end1 - WINDOW_SIZE
   for vec in buffers.buffer0:
      (pChrom, pLoc5, pStrand, primer, umiSeq, alignLocR, alignLocP, readId, read1L, read1R, cigar1, read2L, read2R, cigar2) = vec
      if alignLocR >= buffer1Start:
         buffer1.append(vec)

#----------------------------------------------------------------------------------------------------------------------------------------
# main - mark putative input molecules using BOTH the partially-error-corrected UMI sequence AND the random fragmentation genome position
#----------------------------------------------------------------------------------------------------------------------------------------
def run(cfg):
   print("umi_mark starting...")

   # get params
   readSet  = cfg.readSet
   numCores = cfg.numCores

   # sort read alignments by R2 random fragmentation location (used to identify putative input molecules, in addition to UMI tag)
   fileNameIn = readSet + ".umi_filter.alignments.txt"
   cmd = "sort -k7,7 -k8,8n -k9,9n -k2,2n -t\| -T./ --parallel={1} {0} > {0}.tmp.txt".format(fileNameIn, numCores)  # sort order is (alignChrom, alignStrand, alignLocRand, primerLoc5)
   subprocess.check_call(cmd, shell=True)
   os.rename(fileNameIn + ".tmp.txt", fileNameIn)
   print("umi_mark: done sorting read alignments by random fragmentation position")
   
   # open output file
   fileout1 = open(readSet + ".umi_mark.alignments.txt"  , "w")
   fileout2 = open(readSet + ".umi_mark.for.sum.primer.txt","w")

   # read on-target read pair alignment positions (already sorted by random fragment locus), buffer by locus
   buffers = lambda:0
   buffers.buffer0 = []
   buffers.buffer1 = []
   buffers.mts0 = {}
   buffers.mts1 = {}
   buffers.end1 = None
   chromStrandLast = None
   for line in open(fileNameIn, "r"):
      # unpack read line
      (pChrom, pLoc5, pStrand, primer, umiSeq, isIntendedSite, alignChrom, alignStrand, alignLocR, alignLocP, readId, read1L, read1R, cigar1, read2L, read2R, cigar2) = line.strip().split("|")

      # skip over read alignments from off-target priming
      if isIntendedSite == "0":
         continue
      
      # parse line
      pLoc5 = int(pLoc5)          # primer-side genome position (from design)
      alignLocR = int(alignLocR)  # random-side genome position
      chromStrand = pChrom + "-" + pStrand  # chrom and strand combined
      
      # start new buffer if this is the first read pair alignment record
      if buffers.end1 == None:
         buffers.end1 = alignLocR + WINDOW_SIZE

      # process buffer 1 if this read is past end of accumulation window or chrom/strand change
      if alignLocR >= buffers.end1 or chromStrand != chromStrandLast:
         handleOneLocus(buffers, fileout1, fileout2)

         # do again if current read is beyond end of new window or chrom/strand change (this processes buffer 0, now called buffer 1 after previous handleOneLocus() call)
         if alignLocR >= buffers.end1 or chromStrand != chromStrandLast:
            handleOneLocus(buffers, fileout1, fileout2)
            del(buffers.buffer0[:])
            del(buffers.buffer1[:])
            buffers.end1 = alignLocR + WINDOW_SIZE

            # update chrom/strand if changed
            if chromStrand != chromStrandLast:
               chromStrandLast = chromStrand
      
      # removed redunant fields (alignChrom, alignStrand match the primer design), save the current read
      vec = (pChrom, pLoc5, pStrand, primer, umiSeq, alignLocR, alignLocP, readId, read1L, read1R, cigar1, read2L, read2R, cigar2)
      buffers.buffer1.append(vec)

   # process last two windows
   if buffers.end1 != None:
      handleOneLocus(buffers, fileout1, fileout2)
      handleOneLocus(buffers, fileout1, fileout2)
      
   # done
   fileout1.close()
   fileout2.close()
   print("umi_mark: done")
