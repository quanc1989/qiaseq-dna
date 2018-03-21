import os
import subprocess
import string

def runCopyNumberEstimates(cfg):
   ''' Run CNV analysis using quandico
   '''
   if cfg.runCNV.lower() == "false":
      return
   referenceUmiFiles = cfg.refUmiFiles.split(",")
   assert len(referenceUmiFiles) >= 1, "No reference UMI Files supplied !"
    
   # read reference UMI counts, normalize by read set, and get median for each primer across read sets (not necessary for only one reference read set)
   umiCountsAll = {}
   for fileName in referenceUmiFiles:
      # read MT counts from disk file
      umiCounts = []
      umiCountsTotal = 0
      for line in open(fileName,'r'):
          if not line.startswith("read set|"):
             vals = line.strip().split("|")
             (readSet, primer, strand, chrom, loc5, loc3, umiCount) = vals[0:7]
             (strand, loc5, umiCount) = map(int,(strand, loc5, umiCount))
             key = (chrom,strand,loc5,primer)
             umiCounts.append((key,umiCount))
             umiCountsTotal += umiCount

      # normalize to 1,000 mean MT depth and save in multi-readset hash
      meanUmiDepth = float(umiCountsTotal) / len(umiCounts)
      for key , umiCount in umiCounts:
          if key in umiCountsAll:
             vec = umiCountsAll[key]
          else:
             vec = []
             umiCountsAll[key] = vec
          vec.append(int(round(1000.00 * umiCount / meanUmiDepth)))

   # need to complement first base of negative strand primers
   dnaComplementTranslation = string.maketrans("ATGC", "TACG")
   
   # open output file for quandico input
   readSet = cfg.readSet
   umiFileReference = readSet + ".copy-number.reference.txt"
   fileout = open(umiFileReference, "w")
      
   # for each primer, get median MT depth across all reference read sets and write to disk
   for key, vec in umiCountsAll.iteritems():
      # get median MT depth for this primer
      idx = len(vec) / 2
      if len(vec) % 2 == 1: # odd length
         umiCount = vec[idx]
      else:
         umiCount = int(round((vec[idx-1] + vec[idx]) / 2.00))
              
      # unpack primer info
      (chrom, strand, loc5, primer) = key
           
      # write output in format needed by quandico
      refGenomeBase = primer[0]
      if strand == 1:
         refGenomeBase = refGenomeBase.translate(dnaComplementTranslation)
      outvec = (chrom, loc5+1, strand, refGenomeBase, "foobar", umiCount)
      fileout.write("\t".join((str(x) for x in outvec)))
      fileout.write("\n")
           
   # done writing input reference file
   fileout.close()

   # convert sample readset MT counts to format needed by quandico (NOTE: not normalized - quandico will do that)
   umiFileSample = readSet + ".copy-number.sample.txt"
   fileout = open(umiFileSample,"w")
   fileName = readSet + ".sum.primer.umis.txt"
   for line in open(fileName,"r"):
      if line.startswith("read set|"):
           continue
      vals = line.strip().split("|")
      (readSet_, primer, strand, chrom, loc5, loc3, mtCount) = vals[0:7]
      (strand, loc5, mtCount) = map(int,(strand, loc5, mtCount))
      refGenomeBase = primer[0]
      if strand == 1:
         refGenomeBase = refGenomeBase.translate(dnaComplementTranslation)
      outvec = (chrom, loc5+1, strand, refGenomeBase, "foobar", mtCount)
      fileout.write("\t".join((str(x) for x in outvec)))
      fileout.write("\n")
   fileout.close()

   # make work directory
   if not os.path.exists("_quandico_work_"):
      os.mkdir("_quandico_work_")
   
   # get code dir for quandico
   codeDirQuandico = cfg.quandicoDir
   
   # get reference genome path
   genomeFile = cfg.genomeFile + ".fai"

   # call quandico (note: hack in quandico.pl looks for ".copy-number" in the -b parameter)
   filePrefix = "{}.copy-number".format(readSet)
   cmd = "perl {}quandico.pl ".format(codeDirQuandico) \
       + "-E {}cluster.pl ".format(codeDirQuandico) \
       + "-y {}R/ ".format(codeDirQuandico) \
       + "-t _quandico_work_ " \
       + "-s data={} -s x=2 -s y=0 ".format(umiFileSample)  \
       + "-r data={} -r x=2 -r y=0 ".format(umiFileReference)  \
       + "-G data={} -G name=GRCh37 ".format(genomeFile) \
       + "-b {} ".format(filePrefix) \
       + " > {}.log 2>&1 ".format(filePrefix)
   print(cmd)
   subprocess.check_call(cmd, shell=True)                

def removeNormalVariants(cfg):
   ''' Remove normal variants from tumor vcf
   '''
   readSetNormal = cfg.readSetMatchedNormal
   readSetTumor = cfg.readSet
    
   # do nothing if zero variants from tumor read set
   if not os.path.isfile(readSetTumor + ".smCounter.anno.txt"):
       return
    
   # save normal variants from flat detail file
   normalVariants = set()
   fileName = readSetNormal + ".smCounter.anno.txt"
   if os.path.isfile(fileName):
       with open(fileName,'r') as IN:
           for line in open(fileName,'r'):
               if not line.startswith("CHROM"):
                   vals = line.strip().split("\t")
                   normalVariants.add(tuple(vals[0:5]))
                    
   # filter tumor VCF and corresponding tumor flat file
   for fileSuffix in (".smCounter.anno.txt",".smCounter.anno.vcf"):
       with open(readSetTumor + fileSuffix + ".temp",'w') as OUT:
           with open(readSetTumor + fileSuffix,'r') as IN:
               for line in IN:
                   if line.startswith("#"):
                       OUT.write(line)
                   else:
                       vals = line.strip().split('\t')
                       key = tuple(vals[0:5])
                       if key not in normalVariants:
                           OUT.write(line)
    # replace tumor vcf and flat file with the updated temp files
   for fileSuffix in (".smCounter.anno.txt",".smCounter.anno.vcf"):
       os.system("mv {temp} {f}".format(temp = readSetTumor + fileSuffix + ".temp",
                                        f    = readSetTumor + fileSuffix))

