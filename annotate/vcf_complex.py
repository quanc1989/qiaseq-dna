# Chang Xu. 07MAR2016
import pysam

#----------------------------------------------------------------------------------------------
# Reconstruct complex variants from primitives
#----------------------------------------------------------------------------------------------
def recon(cluster, combined, genomeFile):
   refseq = pysam.FastaFile(genomeFile)
   lastR = -1
   refStr = ''
   refAlt = ''
   for j in combined:
      (tmpChr,tmpPos,tmpID,tmpRef,tmpAlt,tmpPI,tmpFltr,tmpInfo,tmpFmt,tmpGeno,tmpType) = cluster[j]
      if lastR == -1:
         # reconstruct ref and alt sequences
         refStr += tmpRef
         refAlt += tmpAlt
         finalPos = tmpPos # temporary starting position of the complex variant/mnp
         finalPI = int(tmpPI) # temporary prediction index of the complex variant/mnp
         finalID = tmpID 
         finalInfo = tmpInfo # use the info of first variant as the cluster info
         finalFmt  = tmpFmt
         finalGeno = tmpGeno
      elif lastR <= tmpPos:
         refStr = refStr + refseq.fetch(reference=tmpChr, start=lastR, end=tmpPos-1).upper() + tmpRef   #NOTE reference intervals are (a,b]
         refAlt = refAlt + refseq.fetch(reference=tmpChr, start=lastR, end=tmpPos-1).upper() + tmpAlt
      else:
         print('Warning: Last variant right index is greater than the current position')
         break

      # update last variant's right index
      if tmpType in ('SNP', 'Ins'):
         lastR = tmpPos
      elif tmpType == 'Del':
         lastR = tmpPos + len(tmpRef) - 1

      # keep the lowest prediction index as PI for the cluster
      if tmpPI <= finalPI:
         finalPI = tmpPI

      # complex variant type
      finalInfo = finalInfo.replace('SNP', 'COMPLEX').replace('INDEL', 'COMPLEX')

   # remove the leftmost base if necessary
   while len(refStr) > 1 and len(refAlt) > 1:
      if refStr[0] == refAlt[0]:
         refStr = refStr[1:]
         refAlt = refAlt[1:]
         finalPos += 1
      else:
         break

   # output the complex variant/mnp. 
   outline = (tmpChr, finalPos, finalID, refStr, refAlt, finalPI, 'PASS', finalInfo, finalFmt, finalGeno)
   return outline

#----------------------------------------------------------------------------------------------
# Verify complex variants and MNPs
#----------------------------------------------------------------------------------------------
def verifyCluster(bam, cluster, genomeFile):
   # debug check
   if len(cluster) == 1:
      raise Exception("Need at least two variants")

   # open input data files
   refseq = pysam.FastaFile(genomeFile)
   samfile = pysam.AlignmentFile(bam, 'rb')
   
   # check if the adjacent two variants are on the same reads
   outlines = []
   lastVarIsComp = False
   combined = [0]
   for d in range(1,len(cluster)):
      cluster0 = cluster[d-1]
      cluster1 = cluster[d]
      chrom = cluster0[0]
      pos   = (cluster0[ 1], cluster1[ 1])
      #ids  = (cluster0[ 2], cluster1[ 2])
      ref   = (cluster0[ 3], cluster1[ 3])
      alt   = (cluster0[ 4], cluster1[ 4])
      #pIs  = (cluster0[ 5], cluster1[ 5])
      #fltrs= (cluster0[ 6], cluster1[ 6])
      #infos= (cluster0[ 7], cluster1[ 7])
      #fmts = (cluster0[ 8], cluster1[ 8])
      #genos= (cluster0[ 9], cluster1[ 9])
      types = (cluster0[10], cluster1[10])

      # get read alignments supporting each of the two primitive variants
      qname = [set(), set()]
      for i in (0,1):
         locus = "{0}:{1}:{1}".format(chrom, pos[i])
         
         if types[i] == 'SNP':
            endRefPos = pos[i]
            for read in samfile.pileup(region = locus, truncate=True, max_depth=1000000, stepper='nofilter'):
               for pileupRead in read.pileups:
                  if pileupRead.indel == 0 and not pileupRead.is_del:
                     right = pileupRead.alignment.query_sequence[pileupRead.query_position]
                     start = pileupRead.alignment.reference_start + 1
                     end = pileupRead.alignment.reference_end
                     if start != None and end != None and start <= pos[i] and end >= endRefPos and right == alt[i]:
                        if pileupRead.alignment.is_read1:
                           pairOrder = 'R1'
                        if pileupRead.alignment.is_read2:
                           pairOrder = 'R2'
                        qname[i].add(pileupRead.alignment.query_name + ':' + pairOrder)

         elif types[i] == 'Ins':
            endRefPos = pos[i] + 1
            for read in samfile.pileup(region = locus, truncate=True, max_depth=1000000, stepper='nofilter'):
               for pileupRead in read.pileups:
                  if pileupRead.indel > 0:
                     left = pileupRead.alignment.query_sequence[pileupRead.query_position]
                     inserted = pileupRead.alignment.query_sequence[(pileupRead.query_position + 1) : (pileupRead.query_position + 1 +  pileupRead.indel)]
                     right = left + inserted
                     start = pileupRead.alignment.reference_start + 1
                     end = pileupRead.alignment.reference_end
                     if start != None and end != None and start <= pos[i] and end >= endRefPos and right == alt[i]:
                        if pileupRead.alignment.is_read1:
                           pairOrder = 'R1'
                        if pileupRead.alignment.is_read2:
                           pairOrder = 'R2'
                        qname[i].add(pileupRead.alignment.query_name + ':' + pairOrder)

         elif types[i] == 'Del':
            endRefPos = pos[i] + len(ref[i]) - 1  # NOTE: should remove "-1" if we want a flanking base aligned, but this is sometimes soft-clipped primer base, e.g. 78M3D30S
            for read in samfile.pileup(region = locus, truncate=True, max_depth=1000000, stepper='nofilter'):
               for pileupRead in read.pileups:
                  if pileupRead.indel < 0:
                     right = pileupRead.alignment.query_sequence[pileupRead.query_position]
                     deleted = refseq.fetch(reference=chrom, start=pos[i], end=pos[i]+abs(pileupRead.indel))
                     deleted = deleted.upper()
                     left = right + deleted
                     start = pileupRead.alignment.reference_start + 1
                     end = pileupRead.alignment.reference_end
                     if start != None and end != None and start <= pos[i] and end >= endRefPos and right == alt[i]:
                        if pileupRead.alignment.is_read1:
                           pairOrder = 'R1'
                        if pileupRead.alignment.is_read2:
                           pairOrder = 'R2'
                        qname[i].add(pileupRead.alignment.query_name + ':' + pairOrder)

      # get common read set
      commonRead = qname[0].intersection(qname[1])

      # reconstruct complex variants and mnp if adjacent two variants have 95% reads in common
      clusterCond = 1.0 * len(commonRead) / len(qname[0]) >= 0.95 and 1.0 * len(commonRead) / len(qname[1]) >= 0.95

      # if the current variant links with the last variant
      if clusterCond:
         # add variant to current potential combine group
         combined.append(d)

         # if the current variant is the last in the cluster, reconstruct the complex variant and output
         if d == len(cluster) - 1:  
            outlines.append(recon(cluster, combined, genomeFile)) 

      # no link to previous variant
      else:
         # potential combine group has only 1 variant, so output the primitive
         if len(combined) == 1:
            outlines.append(cluster0[0:10])

            # if the current variant is the last in the cluster, also output it as primitive
            if d == len(cluster) - 1:
               outlines.append(cluster1[0:10])

         # end of a potential complex variant
         else:
            outlines.append(recon(cluster, combined, genomeFile))

         # reset combined, starting from the current variant
         combined = [d]

   return outlines

#----------------------------------------------------------------------------------------------
# main function
#----------------------------------------------------------------------------------------------
def run(cfg, bamFileIn, vcfFileIn, vcfFileOut):
   # log start for pipeline use
   print("vcf_complex: start conversion from primitive to complex representation...")

   # get parameters
   gap = int(cfg.vcfComplexGapMax)
   genomeFile = cfg.genomeFile

   # output vcf with mnp and complex variants
   outVcf = open(vcfFileOut, 'w')

   # init variant variables
   lastVariant = None
   cluster = None
   outlines = []
   cnt = 0

   # read VCF input file
   for line in  open(vcfFileIn, 'r'):

      # echo header lines
      if line.startswith('#'):
         outVcf.write(line)
         continue

      # parse input line
      vals = line.strip().split()
      (chrom, pos, id0, ref, alt, pI, fltr, info, fmt, geno) = vals[0:10]
      pos = int(pos)

      # echo filter-failed primitives, and bi-allelic primitives
      if fltr != 'PASS' or alt.find(",") > 0:
         outlines.append((chrom, pos, id0, ref, alt, pI, fltr, info, fmt, geno))
         continue

      # get variant type
      if len(ref) == 1 and len(alt) == 1:
         varType = 'SNP'
      elif len(ref) == 1 and len(alt) > 1:
         varType = 'Ins'
      elif len(ref) > 1 and len(alt) == 1:
         varType = 'Del'
      else:
         raise Exception()

      # save full vector of current variant
      currentVariant = (chrom, pos, id0, ref, alt, pI, fltr, info, fmt, geno, varType)

      # handle start case
      if lastVariant == None:
         cluster = [currentVariant]
         lastVariant = currentVariant
         continue

      # get previous PASS variant info
      (lastChrom, lastPos, lastID, lastRef, lastAlt, lastPI, lastFltr, lastInfo, lastFmt, lastGeno, lastType) = lastVariant

      # conditions on which variants may be clustered. default gap = 3
      cond1 = (lastType in ('SNP', 'Ins') and varType in ('Ins', 'Del') and chrom == lastChrom and abs(pos - lastPos) <= gap)
      cond2 = (lastType == 'SNP' and varType == 'SNP' and chrom == lastChrom and abs(pos - lastPos) <= 1)  # two SNPs have to be next to each other
      cond3 = (lastType == 'Ins' and varType == 'SNP' and chrom == lastChrom and abs(pos - lastPos) - 1 <= gap)
      cond4 = (lastType == 'Del' and varType in ('Ins', 'Del') and chrom == lastChrom and abs(pos - (lastPos + len(lastRef) - 1)) <= gap)
      cond5 = (lastType == 'Del' and varType == 'SNP' and chrom == lastChrom and abs(pos - 1 - (lastPos + len(lastRef) - 1)) <= gap)

      # a candidate cluster
      if cond1 or cond2 or cond3 or cond4 or cond5:

         # add variant to cluster
         cluster.append(currentVariant)

      # else cluster has ended
      else:
         # flush pre-existing cluster to output
         if len(cluster) == 1:
            outlines.append(cluster[0][0:10])
         else:
            outInfo = verifyCluster(bamFileIn, cluster, genomeFile)
            outlines.extend(outInfo)
            cnt += 1

         # init next cluster
         cluster = [currentVariant]

      # save current variant for next iteration
      lastVariant = currentVariant

   # flush last cluster to output
   if cluster == None:
      pass
   elif len(cluster) == 1:
      outlines.append(cluster[0][0:10])
   else:
      outInfo = verifyCluster(bamFileIn, cluster, genomeFile)
      outlines.extend(outInfo)
      cnt += 1

   # report number of clusters checked
   print('vcf_complex: processed ' + str(cnt) + ' variant clusters')

   # sort rows
   outlines.sort()

   # write output VCF
   for row in outlines:
      if row[0] == "chrM":   # hack needed for proper snpSift annotation
         row = list(row)
         row[0] = "chrMT"
      outVcf.write("\t".join((str(x) for x in row)))
      outVcf.write("\n")
   outVcf.close()

   # done if not running as part of longer pipeline
   if "readSet" not in cfg.__dict__:
      return

   # get number of PASS variants called
   numVariantsCalled = 0
   for (chrom, pos, id0, ref, alt, pI, fltr, info, fmt, geno) in outlines:
      if fltr == "PASS":
         numVariantsCalled += 1

   # write count of variants output (used for pipeline integration)
   fileout = open(cfg.readSet + ".vcf_complex.summary.txt","w")
   fileout.write("{}\tvariants called by smCounter\n".format(numVariantsCalled))
   fileout.close()

#----------------------------------------------------------------------------------------------
# pythonism to run from the command line
#----------------------------------------------------------------------------------------------
if __name__ == "__main__":
   import sys
   bamFileIn = sys.argv[1]
   vcfFileIn = sys.argv[2]
   vcfFileOut = sys.argv[3]
   cfg = lambda:0
   cfg.genomeFile = sys.argv[4]
   cfg.vcfComplexGapMax = sys.argv[5]
   run(cfg, bamFileIn, vcfFileIn, vcfFileOut)
