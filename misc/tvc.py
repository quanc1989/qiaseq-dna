import os
import os.path
import subprocess
import socket
import string
import copy

import pysam # only really needed for Ion runs



#-------------------------------------------------------------------------------------------------------
# Ion single-end reads - trim barcode and adapters, drop reads not full length
#-------------------------------------------------------------------------------------------------------
def trimIon(cfg):
   # get params
   readSet = cfg.readSet
   readFile1 = cfg.readFile1
   cutadaptDir = cfg.cutadaptDir
      
   # look for universal seq at 3' end of read, drop read if not found (these are long fragments)
   cmd = cutadaptDir + "cutadapt -e 0.18 -O 17 --discard-untrimmed " \
       + "-a CAAAACGCAATACTGTACATT -n 1 " \
       + "--info-file " + readSet + ".cutadapt.3.R1.txt " \
       + "-o " + readSet + ".temp1.R1.fastq " \
       + readFile1 \
       + " > " + readSet + ".cutadapt.3.R1.log 2>&1 "
   subprocess.check_call(cmd, shell=True)
      
   # read original read count from cutadapt log
   numReadsTotal = None
   for line in open(readSet + ".cutadapt.3.R1.log","r"):
      if line.startswith("Total reads processed:"):
         (key,val) = line.strip().split(":")
         numReadsTotal = int(val.strip().replace(",",""))
         break

   # pull first 12 bp off, move to header
   fileout = open(readSet + ".temp0.R1.fastq","w")
   numReadsWithAdapter3 = 0
   numReadsDroppedTooShort = 0
   lines = []
   lineIdx = 0
   for line in open(readSet + ".temp1.R1.fastq","r"):
      lines.append(line.strip())
      lineIdx += 1

      # all four lines of R1 in memory
      if lineIdx == 4:
         lineIdx = 0
         numReadsWithAdapter3 += 1
      
         # skip reads too short
         if len(lines[1]) < 40:
            numReadsDroppedTooShort += 1
            del(lines[:])
            continue
         
         # get barcode, trim R1 - 12 bp barcode
         barcode  = lines[1][0:12]
         barcodeQ = lines[3][0:12]
         lines[1] = lines[1][12:]
         lines[3] = lines[3][12:]

         # add barcode to R1 id line
         line = lines[0]
         idx = line.find(" ")
         if idx == -1:
            lines[0] = line        + ":" + barcode
         else:
            lines[0] = line[0:idx] + ":" + barcode + line[idx:]
         
         # write output fastq
         for i in range(4):
            fileout.write(lines[i] + "\n")
            
         # clear for next read           
         del(lines[:])
   fileout.close()
   numReadsDroppedNoAdapter3 = numReadsTotal - numReadsWithAdapter3

   # done with temp1 file
   os.remove(readSet + ".temp1.R1.fastq")

   # trim 11-mer universal adapter from 5' end
   cmd = cutadaptDir + "cutadapt -e 0.18 -O 9 --discard-untrimmed " \
       + "-g ^ATTGGAGTCCT -n 1 " \
       + "--info-file " + readSet + ".cutadapt.5.R1.txt " \
       + "-o " + readSet + ".trimmed.R1.fastq " \
       + readSet + ".temp0.R1.fastq" \
       + " > " + readSet + ".cutadapt.5.R1.log 2>&1 "
   subprocess.check_call(cmd, shell=True)
         
   # done with temp0 file
   os.remove(readSet + ".temp0.R1.fastq")

   # get read count from cutadapt log
   numReadsFinal = 0
   for line in open(readSet + ".cutadapt.5.R1.log","r"):
      if line.startswith("Reads with adapters:"):
         (key,val) = line.strip().split(":")
         (val,foo) = val.split("(")
         numReadsFinal = int(val.strip().replace(",",""))
         break
   numReadsDroppedNoAdapter5 = numReadsTotal - numReadsDroppedNoAdapter3 - numReadsDroppedTooShort - numReadsFinal

   # write summary file 
   fileout = open(readSet + ".align.summary.txt", "w")
   fileout.write("{}\tread fragments total\n".format(numReadsTotal))
   fileout.write("{}\tread fragments dropped, not full-length\n".format(numReadsDroppedNoAdapter3))
   fileout.write("{}\tread fragments dropped, less than 40 bp\n".format(numReadsDroppedTooShort))
   fileout.write("{}\tread fragments dropped, universal 11-mer not found\n".format(numReadsDroppedNoAdapter5))
   fileout.close()   

   # set up reverse comlement
   dnaComplementTranslation = string.maketrans("ATGC", "TACG")
   
   # make fake R2 (primer side) file
   fileout = open(readSet + ".trimmed.R2.fastq", "w")
   lineIdx = 0
   for line in open(readSet + ".trimmed.R1.fastq", "r"):
      line = line.strip()
      if lineIdx == 0 or lineIdx == 2:
         fileout.write(line)
      elif lineIdx == 3:
         fileout.write(line[::-1])
      else:
         seq = line[::-1]
         seq = seq.translate(dnaComplementTranslation)
         fileout.write(seq)
      fileout.write("\n")
      lineIdx += 1
      if lineIdx == 4:
         lineIdx = 0
   fileout.close()

#-------------------------------------------------------------------------------------------------------
# Ion single-end reads - align single-end FASTQ using TMAP
#-------------------------------------------------------------------------------------------------------
def alignToGenomeIon(cfg):
   # get some parameters from config
   readSet = cfg.readSet
   numCpus = cfg.numCpus
   samtoolsDir = cfg.samtoolsDir
   samtoolsMem = cfg.samtoolsMem
   torrentBinDir = cfg.torrentBinDir
   torrentGenomeFile = cfg.torrentGenomeFile
   
   # align full-length reads to reference genome using TMAP
   cmd = "{}tmap mapall -n {} -r {}.trimmed.R1.fastq -f {} -v -Y -u --prefix-exclude 5 -o 2 stage1 map4 > ".format(torrentBinDir,numCpus,readSet,torrentGenomeFile) \
   + readSet + ".temp.bam 2> " \
   + readSet + ".align.tmap.log "
   subprocess.check_call(cmd, shell=True)

   # delete unneeded fastq file
   os.remove(readSet + ".trimmed.R1.fastq")
   
   # add a fake reverse compliment read alignment (i.e. simulate paired-end primer-side read) for use in downstream code
   bamIn  = pysam.Samfile(readSet + ".temp.bam", "rb")
   bamOut = pysam.AlignmentFile(readSet + ".align.bam", "wb", template=bamIn)
   for read1 in bamIn:
      read1.is_paired = True
      read1.is_read1 = True
      read1.is_read2 = False
      read2 = copy.deepcopy(read1)
      read2.is_read2 = True
      read2.is_read1 = False
      if not read1.is_unmapped:
         read2.is_reverse = not read1.is_reverse
         read1.mate_is_reverse = read2.is_reverse
         read2.mate_is_reverse = read1.is_reverse
         read1.mate_is_unmapped = False
         read2.mate_is_unmapped = False
      bamOut.write(read1)
      bamOut.write(read2)
   bamIn.close()
   bamOut.close()

   # delete unneeded bam file
   os.remove(readSet + ".temp.bam")
   
   # sort by locus for IGV viewing, and for mtMerge with hashing by chromosome
   cmd = samtoolsDir + "samtools sort -m " + samtoolsMem + " -@" + numCpus \
   + " -T " + readSet \
   + " -o " + readSet + ".align.sorted.bam " \
            + readSet + ".align.bam " \
   + "> "   + readSet + ".align.sort.log 2>&1 "
   subprocess.check_call(cmd, shell=True)
   
   # make BAM index for IGV 
   cmd = samtoolsDir + "samtools index " + readSet + ".align.sorted.bam "
   subprocess.check_call(cmd, shell=True)
