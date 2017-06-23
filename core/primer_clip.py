import os
import os.path

# 3rd party
import pysam

#-----------------------------------------------------------------------
# soft clip one read
#-----------------------------------------------------------------------
def clipRead(read, aligmentReoriented, cigar, basesGenomeNeeded):

   # init return flag for 100% softclip
   dropFlag = False

   # use CIGAR to get number of read bases in the primer binding region that are now aligned to genome
   basesGenome = 0
   basesRead = 0
   for idx in range(0,len(cigar)):
      (op, bases) = cigar[idx]
      if op == 4:    # soft clip - only the read - ignore to get bases of the read aligned
         continue
      elif op == 1:  # insertion to reference - only the read
         basesRead += bases
      elif op == 0: # genome and read in tandem
         bases = min(bases, basesGenomeNeeded - basesGenome) # don't go past the primer 3' binding design site
         basesRead   += bases
         basesGenome += bases
      elif op == 2:  # deletion from reference - only the genome
         basesGenome += bases
      else:
         raise Exception("unexpected CIGAR code")
      if basesGenome >= basesGenomeNeeded:
         break
         
   # debug check
   if basesRead < 0:
      raise Exception("bad CIGAR accounting!")
      
   # if last op was a deletion, move the read clipping one base forward - this could be a misaligned primer base - also forces wipeout of deletion in S walk below
   if op == 2 and (idx + 1) < len(cigar):
      (op, bases) = cigar[idx+1]
      if op == 0:
         basesGenome += 1
         basesRead   += 1
      elif op == 1:
         basesRead   += 1
   
   # if read bases to clip, modify the cigar
   if basesRead > 0:
   
      # increase the soft clip amount, if already some soft clipping
      (op, bases) = cigar[0]
      if op == 4:
         cigar[0] = (op,basesRead+bases)
      
      # otherwise push a new code onto the cigar
      else:
         cigar.insert(0,(4, basesRead))
       
      # walk the cigar, subtract the M,D,I bases overwritten by the new S bases
      basesToSubtract = basesRead
      for idx in range(1,len(cigar)):
         (op, bases) = cigar[idx]
         if op == 2:  # deletion from reference does not help with cigar accounting
            cigar[idx] = (op,0)
         elif bases <= basesToSubtract:
            cigar[idx] = (op,0)
            basesToSubtract -= bases
         else:
            basesNew = bases - basesToSubtract
            cigar[idx] = (op,basesNew)
            basesToSubtract = 0
         if basesToSubtract == 0:
            break

      # debug check
      if basesToSubtract > 0:
         raise Exception("ERROR: entire read soft clipped!", read.qname)
         
      # remove any zero-base ops
      cigarNew = []
      for val in cigar:
         (op, bases) = val
         if bases > 0:
            cigarNew.append(val)
            
      # collapse same adjacent codes (to handle example such as 24S105S)
      cigarTmp = cigarNew
      cigarNew = []
      cigarNew.append(cigarTmp[0])
      for idx in range(1,len(cigarTmp)):
         val = cigarTmp[idx]
         (op,  bases ) = val
         (op_, bases_) = cigarNew[-1]
         if op != op_:
            cigarNew.append(val)
         else:
            bases = bases + bases_
            cigarNew[-1] = (op, bases)
       
      # fix cases in which new alignment starts with an insertion, right after the primer soft clip 
      #  (saw case of G-T priming on sample T insertion variant, causing soft clip too short)
      #  (not sure about the specificity/sensitivity trade off for deletions - this should usually not be a problem)
      if len(cigarNew) > 1:
         (op0, bases0) = cigarNew[0]
         (op1, bases1) = cigarNew[1]
         if op0 != 4:
            raise Exception()
         if op1 == 1:
            cigarNew[0] = (4, bases0 + bases1)  # move insertion bases to soft-clip
            del(cigarNew[1])

      # identify 100% softclipped reads - need to drop these (and PE mate read)
      if len(cigarNew) == 1:
         (op, bases) = cigarNew[0]
         if op == 4:
            dropFlag = True
     
      # reverse the new cigar if negative strand
      if aligmentReoriented:
         cigarNew.reverse()
      
      # update the cigar
      read.cigar = cigarNew
      
      # update the read alignment start position, if it has changed from the cigar edit (aend must be auto recalculated by pysam)
      if not aligmentReoriented:
         read.reference_start += basesGenome

      # change the tlen
      if abs(read.tlen) < basesGenome:
         read.tlen = 0
      elif read.tlen > 0 :
         read.tlen -= basesGenome
      else:
         read.tlen += basesGenome
         
      # debug check on valid cigar string
      numCigarBases = sum((x[1] for x in read.cigar if x[0] != 2))
      if numCigarBases != read.query_length:
         print(read.qname,read.cigarstring, cigar, read.query_length, aligmentReoriented, read.is_read1, basesGenomeNeeded, basesGenome)
         raise Exception("new cigar does not match read length!")
         
      # done
      return dropFlag

#-----------------------------------------------------------------------
# main
#-----------------------------------------------------------------------
def run(cfg, bamFileIn, bamFileOut, resampleOnly):
   print("primer_clip starting...")
   deleteLocalFiles = cfg.deleteLocalFiles
   
   # open files
   bamIn  = pysam.AlignmentFile(bamFileIn , "rb")
   bamOut = pysam.AlignmentFile(bamFileOut, "wb", template=bamIn)
   
   # loop over input BAM - must be sorted by read id
   numReadPairs = 0
   numReadPairsTrimmed = 0
   numReadPairsDropped = 0
   for read1 in bamIn:

      # crash if read is not paired
      if not read1.is_paired:
         raise Exception("read not paired! " + read1.qname)
         
      # get mate, assuming mate is the next record in the BAM file - alignments assume to be filtered upstream
      read2 = bamIn.next()
      if read1.qname != read2.qname:
         print(read1.qname, read2.qname)
         raise Exception("read mate is not next in BAM record order!")
         
      # switch variable names to match actual R1 and R2
      if read1.is_read2:
         temp = read2
         read2 = read1
         read1 = temp
      
      # report progress
      numReadPairs += 1
      if numReadPairs % 1000000 == 0:
         print("# read pairs processed: " + "{0:,}".format(numReadPairs))

      # check for internal-resample-only request
      if resampleOnly and read1.get_tag("re") == 0:
         bamOut.write(read1)
         bamOut.write(read2)
         continue
         
      # flag for pair with one or both reads 100% soft clipped
      dropReadPair = False
         
      # loop over R1 then R2
      for read in (read1, read2):
      
         # get tag containing primer design location
         designSite = read.get_tag("pr")
         
         # parse primer design site
         (chrom, strand, loc5, primerLen) = designSite.split("-")
         primerLen = int(primerLen)
         strand = int(strand)
         loc5 = int(loc5)
         loc3 = loc5 + primerLen - 1 if strand == 0 else loc5 - primerLen + 1
     
         # save the read alignment orientation
         aligmentReoriented = (read.is_reverse and read.is_read1) or (not read.is_reverse and read.is_read2)
         
         # copy the cigar list, reverse the cigar when necessary to put primer at start
         if aligmentReoriented:
            cigar = read.cigar[::-1]
            basesGenomeNeeded = read.aend - loc3
         else:
            cigar = list(read.cigar)
            basesGenomeNeeded = loc3 - read.pos + 1
         
         # soft clip the primer bases
         if basesGenomeNeeded > 0:
            dropFlag = clipRead(read, aligmentReoriented, cigar, basesGenomeNeeded)
            if dropFlag:
               dropReadPair = True
               continue
               
            # count reads with R1 trimmed
            if read.is_read1:
               numReadPairsTrimmed += 1

         # not sure whether or not variant callers use the MD tag, but just delete it for now (can be filled via samtools fillmd)
         if read.has_tag("MD"):
            read.set_tag("MD",None)

      #---------------------------------------------------------------------------------------------------
      # clip any umi/common-11-mer oligo from 3' end of R1 that aligns beyond the R2 alignment start
      #---------------------------------------------------------------------------------------------------
      
      # save the read alignment orientation - opposite of primer trim
      aligmentReoriented = not aligmentReoriented
      
      # copy the cigar list, reverse the cigar when necessary to put umi/common-11-mer at start
      if aligmentReoriented:
         cigar = read1.cigar[::-1]
         basesGenomeNeeded = read1.aend - read2.aend
      else:
         cigar = list(read1.cigar)
         basesGenomeNeeded = read2.pos - read1.pos
   
      # soft clip the umi/common-11-mer oligo bases from R1 3' end
      if basesGenomeNeeded > 0:
         dropFlag = clipRead(read1, aligmentReoriented, cigar, basesGenomeNeeded)
         if dropFlag:
            dropReadPair = True
            continue
         
      # write the modified read pair to the output BAM
      if dropReadPair:
         numReadPairsDropped += 1
      else:
         bamOut.write(read1)
         bamOut.write(read2)

   # done         
   bamIn.close()
   bamOut.close()

   # report debug count
   print("# read pairs total", numReadPairs)
   print("# read pairs dropped because 100% soft clipped:", numReadPairsDropped)
   print("# read pairs with R1 trimmed: {}".format(numReadPairsTrimmed))
   
   # delete input file if local
   if deleteLocalFiles and len(os.path.dirname(bamFileIn)) == 0:
      os.remove(bamFileIn)
