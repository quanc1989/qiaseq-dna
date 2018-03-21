import subprocess
import string
import os
import array
import shutil

# third party modules
import pysam

# our modules
import bed

#----------------------------------------------------------------------------------------------------------   
# main TVC function - prepare primer-trimmed BAM and call TVC
#----------------------------------------------------------------------------------------------------------   
def run(cfg):
   # get params
   print("tvc: start...")
   readSet = cfg.readSet
   uBam = cfg.uBam
   numCpus = cfg.numCores
   samtoolsMem = cfg.samtoolsMem
   samtoolsDir = cfg.samtoolsDir
   vcflibDir = cfg.vcflibDir
   flowTagsNeeded = set(("ZP","ZA","ZG","ZB","ZC","ZM","ZF","RG"))
   
   # open raw read file *.basecaller.bam
   bam = pysam.AlignmentFile(uBam, "rb", check_sq = False)
   
   # save BAM header in a string for later
   bamHeaderRaw = bam.text
   
   # open output file, and cutadapt 3' trim file
   fileout = open(readSet + ".tvc.flowtags.txt","w")
   fileIn3 = open(readSet + ".cutadapt.3.R1.txt","r")

   # create a dict of read id -> tag
   umi_dict = {}
   tag_name = "mi"
   with open(readSet +  ".umi.tag.txt","r") as IN:
      for line in IN:
         read_id, umi, umi_qual = line.strip('\n').split('\t')
         umi_dict[read_id] = umi

   print "\nDone creating readID -> UMI  dict\n"
   
   # merge 5' trim, 3' trim, and raw read flow tags
   for line in open(readSet + ".cutadapt.5.R1.txt","r"):
      vals5 = line.strip().split("\t")
      readId = vals5[0]      
      umi = umi_dict["@"+readId]
      # spin the 3' trim file, and the raw read file forward
      while True:
         line = fileIn3.readline()
         vals3 = line.strip().split("\t")
         readId3 = vals3[0]
         read = bam.next()
         if read.query_name != readId3:
            raise Exception("3' cutadapt trim info file not in same order as raw read file")
         if readId3 == readId:
            break

      # debug check            
      if readId3 != readId:
         raise Exception("3' and 5' cutadapt trim files out of sync")

      # skip to next read if 5' adapter not found
      if int(vals5[1]) == -1:
         continue
      
      # debug check
      if int(vals3[1]) == -1:
         raise Exception("trim of both 5' and 3' adapters expected")
         
      # fix line strip
      while len(vals3) < 11:
         vals3.append("")
      while len(vals5) < 11:
         vals5.append("")

      # avoid hassles of pysam type conversions for optional tags, by using SAM text format (might be slow)
      samVals = read.tostring(bam).split("\t")
      
      # get barcode region quality from raw read
      umiC = samVals[ 9][0:12]
      umiQ = samVals[10][0:12]
      if umiC != umi:
         raise Exception("unexpected barcode sync")

      # get trimmed sequences
      seq3  = vals3[5] + vals3[6]             # 6 should be empty string
      seq5  = umi  + vals5[4] + vals5[5]  # 4 should be empty string
      
      # get trimmed qual vals
      qual3 = vals3[9] + vals3[10]           # 10 should be empty string
      qual5 = umiQ + vals5[8] + vals5[9] #  8 should be empty string
      
      # debug check on read length
      readLenRaw = len(samVals[9])
      readLenTrm = len(seq5) + len(vals5[6]) + len(seq3)
      if readLenTrm != readLenRaw:
         raise Exception("Length of trimmed read not equal to raw read!")
      
      # get raw read flow signal tags
      outvec = [readId, seq5, qual5, seq3, qual3]
      for tag in samVals[11:]:
         if tag[0:2] in flowTagsNeeded:
            outvec.append(tag)

      # write to disk
      fileout.write("|".join(outvec))
      fileout.write("\n")
      
   # done
   fileout.close()
   fileIn3.close()
   bam.close()
   
   # sort the trimmed seq / flow tag file by read id
   cmd = "sort -k1,1 -t\| --parallel={0} {1}.tvc.flowtags.txt > {1}.tvc.flowtags.sorted.txt".format(numCpus,readSet)
   subprocess.check_call(cmd, shell=True)
   os.remove("{}.tvc.flowtags.txt".format(readSet))
   
   # sort oligoClip file by read id
   cmd = samtoolsDir + "samtools sort -n -m " + samtoolsMem + " -@" + numCpus \
   + " -T " + readSet \
   + " -o " + readSet + ".tvc.temp.bam " \
            + readSet + ".bam " \
   + " > "  + readSet + ".tvc.sort.log 2>&1 "
   subprocess.check_call(cmd, shell=True)

   # set up reverse comlement
   dnaComplementTranslation = string.maketrans("ATGC", "TACG")
   
   # open readId-sorted main BAM file, build header for TVC output bam
   bamIn = pysam.AlignmentFile(readSet + ".tvc.temp.bam", "rb")
   
   # dump bam header to sam file, because we are using older version of pysam that cannot directly take text lines for header
   headerTagsNeeded = set(["CO", "RG", "PG"])
   fileOut = open(readSet + ".tvc.header.sam", "w")
   for line in bamIn.text.split("\n"): # init with TMAP tags
      if len(line) > 3 and line[1:3] != "RG":
         fileOut.write(line)
         fileOut.write("\n")
   for line in bamHeaderRaw.split("\n"): # add Ion BaseCaller flow tags
      if line[1:3] in headerTagsNeeded:
         fileOut.write(line)
         fileOut.write("\n")
   fileOut.close()
   samHeaderOnly = pysam.AlignmentFile(readSet + ".tvc.header.sam", "r")
   
   # open output BAM file
   bamOut = pysam.AlignmentFile(readSet + ".tvc.bam", "wb", template=samHeaderOnly)
   samHeaderOnly.close()
   os.remove(readSet + ".tvc.header.sam")
   
   # make TVC input file - add hard clipped regions back as soft-clipped alignments
   fileIn = open(readSet + ".tvc.flowtags.sorted.txt","r")
   for read in bamIn:
   
      # drop fake R2 (primer side)
      if not read.is_read2:
         continue
         
      # change back to single end
      read.is_read1 = False
      read.is_read2 = False
      read.is_paired = False
      read.mate_is_reverse = False
      read.mate_is_unmapped = False
      
      # parse readId
      vals = read.query_name.split(":")
      #readIdBam = ":".join(vals[0:-2])
      readIdBam = read.query_name
      
      # spin the flowtag file forward (not all reads in the bam)
      readId = None
      while True:
         line = fileIn.readline()
         vals = line.strip().split("|")
         (readId, seq5, qual5, seq3, qual3) = vals[0:5]
         if readId == readIdBam:
            break
      
      # debug check
      if readId == None:
         raise Exception("missing read id in TVC flowtag merge")
      
      # handle negative strand alignment
      if read.is_reverse:
         tmp  = seq5
         seq5 = seq3
         seq3 = tmp
         seq5 = seq5[::-1]
         seq5 = seq5.translate(dnaComplementTranslation)
         seq3 = seq3[::-1]
         seq3 = seq3.translate(dnaComplementTranslation)
         tmp   = qual5[::-1]
         qual5 = qual3[::-1]
         qual3 = tmp
      
      # copy the cigar      
      cigar = list(read.cigar)
      
      # add 5' trim back on
      (op, bases) = cigar[0]
      if op == 4:
         cigar[0] = (op, bases + len(seq5))
      else:
         cigar.insert(0,(4, len(seq5)))
   
      # add 3' trim back on
      (op, bases) = cigar[-1]
      if op == 4:
         cigar[-1] = (op, bases + len(seq3))
      else:
         cigar.append((4, len(seq3)))
      
      # save cigar edits
      read.cigar = cigar
      
      # pysam requires saving qual values first
      qual = read.qual
      
      # fix up the seq 
      read.query_sequence = seq5 + read.query_sequence + seq3

      # fix up the quality
      read.qual = qual5 + qual + qual3
      
      # add flow quality tags
      for tag in vals[5:]:
         (tagName,tagType,tagVal) = tag.split(":")
         if tagType == "Z":
            pass
         elif tagType == "i":
            tagVal = int(tagVal)
         elif tagType == "B":
            if tagVal.startswith("f,"):
               tagVal = array.array("f",[float(x) for x in tagVal[2:].split(",")])
            elif tagVal.startswith("i,"):
               tagVal = array.array("i",[int(x)   for x in tagVal[2:].split(",")])
            elif tagVal.startswith("s,"):
               tagVal = array.array("h",[int(x)   for x in tagVal[2:].split(",")])
            else:
               raise Exception()
         else:
            raise Exception()
         read.set_tag(tagName,tagVal)
         
      # output modified read
      bamOut.write(read)
      
   # done
   fileIn.close()
   bamIn.close()
   bamOut.close()
   os.remove(readSet + ".tvc.temp.bam")
   
   # sort final TVC input bam
   cmd = samtoolsDir + "samtools sort -m " + samtoolsMem + " -@" + numCpus \
   + " -T " + readSet \
   + " -o " + readSet + ".tvc.sorted.bam " \
            + readSet + ".tvc.bam " \
   + " > "  + readSet + ".tvc.sort.log 2>&1 "
   subprocess.check_call(cmd, shell=True)
   
   # index final TVC input bam
   cmd = samtoolsDir + "samtools index " + readSet + ".tvc.sorted.bam"
   subprocess.check_call(cmd, shell=True)
   
   # run TVC
   roiBedFile = cfg.roiBedFile
   torrentBinDir     = cfg.torrentBinDir
   torrentGenomeFile = cfg.genomeFile
   torrentVcfFile = readSet + ".tvc.vcf"
   cmd = os.path.join(torrentBinDir , "tvc") + " --output-dir _TVC_ " \
    + " -n " + numCpus \
    + " -b " + readSet + ".tvc.sorted.bam" \
    + " -t " + roiBedFile \
    + " -r " + torrentGenomeFile \
    + " -o " + torrentVcfFile \
    + " --snp-min-allele-freq 0.005" \
    + " --snp-min-cov-each-strand 0 " \
    + " --snp-min-coverage 3" \
    + " --snp-min-var-coverage 2" \
    + " --snp-min-variant-score 6" \
    + " --snp-strand-bias 1" \
    + " --snp-strand-bias-pval 0" \
    + " --mnp-min-allele-freq 0.005" \
    + " --mnp-min-cov-each-strand 0" \
    + " --mnp-min-coverage 3" \
    + " --mnp-min-var-coverage 2" \
    + " --mnp-min-variant-score 6" \
    + " --mnp-strand-bias 1" \
    + " --mnp-strand-bias-pval 0" \
    + " --indel-min-allele-freq 0.05" \
    + " --indel-min-cov-each-strand 0" \
    + " --indel-min-coverage 3" \
    + " --indel-min-var-coverage 2" \
    + " --indel-min-variant-score 10" \
    + " --indel-strand-bias 1" \
    + " --indel-strand-bias-pval 0" \
    + " > " + readSet + ".tvc.log 2>&1"
   print("tvc: command line is " + cmd) 
   subprocess.check_call(cmd, shell=True)
   print("tvc: done running TVC")
   
   # move TVC VCF to current directory
   os.rename("_TVC_/" + torrentVcfFile, torrentVcfFile)
   
   # call up vcflib command to split multi-allelic, remove GT tag, get primitives
   cmd = "{0}vcfbreakmulti {1}.tvc.vcf | " \
       + "{0}vcfkeepgeno - DP AF AD VF | " \
       + "{0}vcfallelicprimitives --tag-parsed AP > {1}.tvc.primitives.vcf 2> {1}.tvc.vcflib.log"
   cmd = cmd.format(vcflibDir,readSet)
   subprocess.check_call(cmd,shell=True)
   
   # (1) drop TVC primitive variants that have allele fraction below 0.05
   # (2) make BED file for smCounter - regions +/- 10 bp from a TVC primitive variant
   bedTvc = []
   fileout   = open(readSet + ".tvc.primitives.temp.vcf", "w")
   for line in open(readSet + ".tvc.primitives.vcf", "r"):

      # echo VCF header
      if line.startswith("#"):
         fileout.write(line)
         continue
         
      # parse line
      chrom, pos, id, ref, alt, qual,filter,info,format,sampleId = line.strip().split("\t")

      # make sure data is as expected      
      if alt.find(",") >= 0:
         raise Exception("tvc: not expecting multi-allelic variant in primitives file")
      
      # get left location, zero-based
      locL = int(pos) - 1
      
      # look for indel, include right flanking base
      altLen = len(alt)
      refLen = len(ref)
      if altLen == refLen:  # SNP or MNP
         locR = locL + refLen
         isIndel = False
      else:  # INDEL
         locR = locL + refLen + 1
         isIndel = True
      
      # get AF tag - assume always present - exception will be thrown if None remains
      alleleFraction = None
      for tag in info.split(";"):
         if tag.find("=") > 0:
            tagName, tagVal = tag.split("=")
            if tagName == "AF":
               if tagVal.find(",") >= 0:
                  raise Exception("tvc: not expecting TVC primitives to be multi-allelic")
               alleleFraction = float(tagVal)
                  
      # drop TVC primitive variants with low INDEL allele fraction
      if (isIndel and alleleFraction < 0.05) or alleleFraction < 0.005:
         continue
         
      # echo line to new TVC VCF primitives file
      fileout.write(line)
      
      # save region, with 10 bp flanking 
      locL = max(0,locL - 10)
      locR += 10
      if chrom == "chrM" and locR > 16569:  # horrific hack for chrM NC_012920 reference
         locR = 16569
      if locL < locR:
         bedTvc.append((chrom,locL, locR))
      
   # close filtered TVC VCF primitives file, rename for later use
   fileout.close()
   os.rename(readSet + ".tvc.primitives.temp.vcf", readSet + ".tvc.primitives.vcf")
   
   # merge BED and write to disk
   bedTvc = bed.merge(bedTvc)
   bed.write(bedTvc, readSet + ".tvc_roi.bed")
   print("tvc: done running TVC and making smCounter ROI bed")



#----------------------------------------------------------------------------------------------------------   
# use TVC variant list to filter smCounter variant list
#----------------------------------------------------------------------------------------------------------   
def smCounterFilter(cfg,vc):

   # get params
   readSet = cfg.readSet
   print("tvc: start smCounterFilter")
   
   # store TVC primitives in hash - all have PASS flag
   tvcPrimitives = set()
   for line in open(readSet + ".tvc.primitives.vcf", "r"):
      if not line.startswith("#"):
         chrom, pos, id, ref, alt, qual,filter,info,format,sampleId = line.strip().split("\t")
         key = (chrom, pos, ref, alt)
         tvcPrimitives.add(key)
   print("tvc: TVC primitive variants: {}".format(len(tvcPrimitives)))
   
   # loop over smCounter variants, for each check if all primitives are in the TVC primitives hash
   lowQSuffix = ".smCounter.GT12PI" if vc.lower() == "v1" else ".smCounter.lowQ"

   shutil.copyfile(readSet + lowQSuffix + ".txt", readSet + lowQSuffix + ".tmp.txt")
   fileIn2  = open(readSet + ".smCounter.cut.txt", "r")
   fileOut1 = open(readSet + ".smCounter.cut.tmp.vcf", "w")
   fileOut2 = open(readSet + ".smCounter.cut.tmp.txt", "w")
   fileOut3 = open(readSet + lowQSuffix + ".tmp.txt", "a")
   numVariantsRetained = 0
   numVariantsFiltered = 0
   for line in open(readSet + ".smCounter.cut.vcf", "r"):
      if line.startswith("##"):
         fileOut1.write(line)
      elif line.startswith("#CHROM"):
         fileOut1.write(line)
         fileOut2.write(fileIn2.readline())
      else:
         # read the detail file in parallel with the VCF file
         line2 = fileIn2.readline()
         
         # parse smCounter VCF line, and flat file line
         chrom,pos,id,ref,alt,qual,filter,info,format,sampleId = line.strip().split("\t")


         if vc.lower() == "v1":
            CHROM,POS,REF,ALT,TYPE,DP,MT,UMT,PI,THR,VMT,VMF,VSM,FILTER = line2.strip().split("\t")
         else:
            CHROM,POS,REF,ALT,TYPE,DP,VDP,VAF,UMT,VMT,VMF,QUAL,FILTER = line2.strip().split("\t")
         
         # debug check
         if chrom != CHROM or pos != POS or ref != REF or alt != ALT or filter != FILTER:
            raise Exception("tvc: smCounter VCF and TXT out of sync!!!")
            
         # prevent possible coding mistakes below
         CHROM = None
         POS = None
         REF = None
         ALT = None
         FILTER = None

         # handle not multi-allelic variant (which is primitive, not MNP in smCounter code)
         if alt.find(",") == -1:
            key = (chrom, pos, ref, alt)
            keepVariant = key in tvcPrimitives
         
         #check for multi-allelic case (
         #NOTE: smCounter only outputs bi-allelic heterogygous germline case, and only gives detail for highest PI variant - only ALT, TYPE, GT, and AD have multi-allelic info - FILTER always the same for both)

         #NOTE2 : smCounter-v2 does not have restrictions on the zygocity, 

         else:
            alt1, alt2 = alt.split(",")  # smCounter only gives 2 ALTs
            key1 = (chrom, pos, ref, alt1)
            key2 = (chrom, pos, ref, alt2)
            
            # drop both alleles - neither are in TVC primitives
            if key1 not in tvcPrimitives and key2 not in tvcPrimitives:
               keepVariant = False
            
            # keep both alleles - both are in TVC primitives
            elif key1 in tvcPrimitives and key2 in tvcPrimitives:
               keepVariant = True
               
            # need to UNDO the smCounter 2-allele hacks - only in the ALT, TYPE, GT, and AD fields - this code must stay in sync with smCounter code!!!!!
            else:
               # keeping the smCounter variant, but removing one of the two alleles
               keepVariant = True               
              
               # fix ALT
               alt = alt1 if key1 in tvcPrimitives else alt2               

               if vc.lower == 'v1':
                  # fix SAMPLE - note that VMF was only ever output for the first ALT
                  GT, AD, VMF_ = sampleId.split(":")
                  if not AD.endswith(",1") or GT != "1/2" or VMF_ != VMF:
                     raise Exception("tvc: unexpected bi-allelic formats in smCounter VCF")
                  AD = AD[:-2]
                  GT = "1" if chrom == "chrY" or chrom == "chrM" else "0/1"
                  sampleId = ":".join((GT,AD,VMF))
                  
                  # fix INFO (in VCF) and TYPE (in flat file) - must keep in sync with smCounter code
                  type1, type2 = TYPE.split(",")
                  TYPE = type1 if key1 in tvcPrimitives else type2
                  info = ';'.join(('TYPE='+TYPE, 'DP='+DP, 'MT='+MT, 'UMT='+UMT, 'PI='+PI, 'THR='+THR, 'VMT='+VMT, 'VMF='+VMF, 'VSM='+VSM))

               else: ## SmCounter-v2 has different metrics
                  if key1 in tvcPrimitives:
                     index = 0
                  else:
                     index = 1

                  # fix SAMPLE
                  GT, AD, VMF_ = sampleId.split(":")
                  if not AD.endswith(",1") or GT != "1/2":
                     raise Exception("tvc: unexpected bi-allelic formats in smCounter-v2 VCF") ## Even if > 2 alt alleles, AD is still reported as 1/2 in smCounter-v2
                  AD = AD[:-2]
                  GT = "1" if chrom == "chrY" or chrom == "chrM" else "0/1"
                  sampleId = ":".join((GT,AD,VMF_.split(",")[index]))

                  ## Fix info
                  repeat = info.split(";")[1]
                  umt = info.split(';')[3]
                  assert umt == sUMT, "tvc: could not sync UMT values for smCounter-v2 vcf and variants file.\nVcf:{v1} VariantsFile:{v2}".format(v1=umt,v2=sUMT)
                  
                  ## split and use the appropriate metric corresponding to the allele being output
                  TYPE = TYPE.split(",")[index]
                  DP = DP.split(",")[index]
                  VMT = VMT.split(",")[index]
                  VMF = VMF.split(",")[index]
                  
                  info = ';'.join(["TYPE="+TYPE,repeat,"DP="+DP,"UMT="+UMT,"VMT="+VMT,"VMF="+VMF])
                  
         # remove smCounter filter if TVC called the variant, and VMF >= 0.30
         if keepVariant and filter != "PASS" and float(VMF) >= 0.30:
            filter = "PASS"
            
         # add smCounter FILTER if we are dropping this variant
         elif not keepVariant:
            if filter == "PASS":
               filter = "NO_TVC_CALL"
            else:
               filter = filter + ";NO_TVC_CALL"
               
         # reconstruct the possibly edited output lines
         line  = "\t".join((chrom,pos,id,ref,alt,qual,filter,info,format,sampleId)) + "\n"
         if vc.lower() == "v1":
            line2 = "\t".join((chrom,pos,ref,alt,TYPE,DP,MT,UMT,PI,THR,VMT,VMF,VSM,filter)) + "\n"
         else:
            line2 = "\t".join((chrom,pos,ref,alt,TYPE,DP,VDP,VAF,UMT,VMT,VMF,QUAL,filter)) + "\n"
         
         # smCounter primitives match TVC - keep this variant
         if keepVariant:
            fileOut1.write(line)
            fileOut2.write(line2)
            numVariantsRetained += 1
            
         # variant not in TVC primitives - move to GT12PI/lowQ supplement text file
         else:
            fileOut3.write(line2)
            numVariantsFiltered += 1
            print("tvc: filtering smCounter variant: " + line.strip())
            
   # close files
   fileIn2.close()
   fileOut1.close()
   fileOut2.close()
   fileOut3.close()

   # rename the files for downstream use
   for suffix in (".smCounter.cut.tmp.vcf",".smCounter.cut.tmp.txt",lowQSuffix+".tmp.txt"):
      fileNameSrc = readSet +  suffix
      fileNameDes = fileNameSrc.replace(".tmp.", ".")
      fileNameBak = fileNameSrc.replace(".tmp.", ".old.")
      shutil.copyfile(fileNameDes,fileNameBak)
      shutil.move(fileNameSrc, fileNameDes)
   
   # done
   print("tvc: # variants filtered: {}".format(numVariantsFiltered))
   print("tvc: # variants retained: {}".format(numVariantsRetained))
   print("tvc: done smCounterFilter")
   return numVariantsRetained

#----------------------------------------------------------------------------------------------------------   
# command line test
#----------------------------------------------------------------------------------------------------------   
if __name__ == '__main__':
   cfg = lambda:0
   cfg.readSet = "S5-0333-3.subset"
   cfg.vcflibDir = "/srv/qgen/bin/vcflib/bin/"
   smCounterFilter(cfg)
