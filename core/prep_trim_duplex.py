import os
import subprocess


#-------------------------------------------------------------------------------------------------------
# trimming using recognition of universal sequence on both sides, and extract duplex tag TT or CC
#-------------------------------------------------------------------------------------------------------
def trimDuplex(cfg):
   print("prep: starting read prep duplex - trimming ends and UMI extraction")
   # get parameters
   readSet = cfg.readSet
   readFile1 = cfg.readFile1
   readFile2 = cfg.readFile2
   cutadaptDir = cfg.cutadaptDir

   filePrefixOut = readSet + ".prep"
   
   # validate existence of read files
   for fileName in (readFile1, readFile2):
      if not os.path.isfile(fileName):
         raise UserWarning("input read file does not exist: " + fileName)

   # trim R2 reads (primer side), 3' end (also gets 12 bp barcode and ILMN adapter not cut by ILMN software)
   cmd  = cutadaptDir + "cutadapt -e 0.18 -O 3 " \
        + "-a AGGACTCCAAT -n 1 " \
        + "-o " + filePrefixOut + ".temp3.R2.fastq " \
        + readFile2 \
        + " > " + filePrefixOut + ".cutadapt.3.R2.log 2>&1 " 
   subprocess.check_call(cmd, shell=True)
   #deleteLocalFastq(readFile2)

   # trim R1 reads (barcode side), 3' end (gets ILMN adapter not cut by ILMN software - customer sequencing primer settings usually wrong)
   cmd  = cutadaptDir + "cutadapt -e 0.18 -O 3 " \
        +  "-a CAAAACGCAATACTGTACATT -n 1 " \
        + "-o " + filePrefixOut + ".temp3.R1.fastq " \
        + readFile1 \
        + " > " + filePrefixOut + ".cutadapt.3.R1.log 2>&1 "
   subprocess.check_call(cmd, shell=True)
   
   # trim both 5' ends using paired-end-mode, dropping reads that do not have both univerals
   cmd = cutadaptDir + "cutadapt -e 0.18 -O 18 --discard-untrimmed --minimum-length 35 " \
       + "-u 12 -g ^TTCTGAGCGAYYATAGGAGTCCT -G ^AATGTACAGTATTGCGTTTTG " \
       + "-o {0}.temp5.R1.fastq -p {0}.temp5.R2.fastq {0}.temp3.R1.fastq {0}.temp3.R2.fastq > {0}.cutadapt.5.log 2>&1 ".format(filePrefixOut)
   subprocess.check_call(cmd, shell=True)
   
   # delete unneeded fastq file
   os.remove(filePrefixOut + ".temp3.R2.fastq")
   
   # init counters
   numReadPairsTotal = 0
   numReadPairsDropped = 0
   numReadPairsTT = 0
   numReadPairsCC = 0
   numReadPairsNN = 0
   
   # open input and output files
   fileout1 = open(readSet + ".prep.R1.fastq","w")
   fileout2 = open(readSet + ".prep.R2.fastq","w")
   filein1  = open(filePrefixOut + ".temp5.R1.fastq","r")
   filein2  = open(filePrefixOut + ".temp5.R2.fastq","r")
   filein3  = open(filePrefixOut + ".temp3.R1.fastq","r")  # needed to pull barcode
   
   # loop over all R1 reads
   lineIdx = 0
   for line1 in filein1:
      line2 = filein2.readline()
      
      # echo lines 2,3,4
      if lineIdx > 0:
         fileout1.write(line1)
         fileout2.write(line2)
         
      # add read id and duplex tag to line 1 readId         
      else:
         idx = line1.find(" ")
         readId1 = line1[:idx]
         readId2 = line2[:idx]
         assert readId2 == readId1,"R1 and R2 do not have the same name"
         
         # get the read sequence from temp3.R1.fastq file, which contains all reads, all 5' bases
         line3 = None
         while True:
            line3 = filein3.readline()
            numReadPairsTotal += 1
            i = line3.find(" ")
            readId3 = line3[:i]
            line3 = filein3.readline()
            filein3.readline() # discard line 3
            filein3.readline() # discard line 4
            if readId3 == readId1:
               break
            else:
               numReadPairsDropped += 1
         if line3 == None:
            raise Exception()

         # extract barcode and duplex tag (rough initial cut)
         umi = line3[0:12]
         duplexTagRegion = line3[21:25]
         if duplexTagRegion.find("CC") >= 0:
            duplexTag = "CC"
            numReadPairsCC += 1
         elif duplexTagRegion.find("TT") >= 0:
            duplexTag = "TT"
            numReadPairsTT += 1
         else:
            duplexTag = "NN"
            numReadPairsNN += 1
            
         # add the duplex tag to read id and the umi tag  after the space so that bwa can create a seperate tag for it
         # note : bwa mem -C does not seem to create 2 tags
         umiTag = cfg.tagNameUmiSeq + ":Z:" + umi
         duplexTag = "DU:Z:" + duplexTag
         
         line1 = readId1 + ":" +  duplexTag + " " + umiTag
         line2 = readId1 + ":" +  duplexTag + " " + umiTag
         fileout1.write(line1+"\n")
         fileout2.write(line2+"\n")
      
      # prepare next iteration
      lineIdx += 1
      if lineIdx == 4:
         lineIdx = 0
         
   # close files
   fileout1.close()
   fileout2.close()
   filein1.close()
   filein2.close()
   filein3.close()
   
   # write summary file 
   fileout = open(readSet + ".prep.summary.txt", "w")
   fileout.write("{}\tread fragments total\n".format(numReadPairsTotal))
   fileout.write("{}\tread fragments dropped, less than 40 bp\n".format(numReadPairsDropped))
   fileout.write("{}\tread fragments with duplex tag CC\n".format(numReadPairsCC))
   fileout.write("{}\tread fragments with duplex tag TT\n".format(numReadPairsTT))
   fileout.write("{}\tread fragments with duplex tag NN\n".format(numReadPairsNN))
   fileout.close()   
      
   # delete unneeded fastq files
   os.remove(filePrefixOut + ".temp3.R1.fastq")
   os.remove(filePrefixOut + ".temp5.R1.fastq")
   os.remove(filePrefixOut + ".temp5.R2.fastq")
   
