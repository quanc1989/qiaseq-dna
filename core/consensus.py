import os
import os.path
import subprocess

import pysam

#-------------------------------------------------------------------------
# make unique molecule consensus reads FASTQ
#-------------------------------------------------------------------------
def run(cfg,bamFileIn,primerFileIn):
   print("consensus: starting...")

   # get params
   readSet          = cfg.readSet
   javaExe          = cfg.javaExe
   fgbioJar         = cfg.fgbioJar
   samtoolsDir      = cfg.samtoolsDir
   samtoolsMem      = cfg.samtoolsMem
   numCores         = cfg.numCores
   genomeFile       = cfg.genomeFile
   deleteLocalFiles = cfg.deleteLocalFiles
   tagNameUmi       = cfg.tagNameUmi

   # sort the BAM file by column 12 (the "Mi" unique molecule identifier tag)
   cmd = samtoolsDir + "samtools sort"  \
   + " -m " + samtoolsMem \
   + " -@ " + numCores \
   + " -T " + readSet \
   + " -t " + tagNameUmi \
   + " -o " + readSet + ".consensus.temp0.bam" \
   + "    " + bamFileIn \
   + " >> " + readSet + ".consensus.shell.log 2>&1"
   subprocess.check_call(cmd, shell=True)
   print("consensus: done sorting BAM by 'Mi' tag")
   
   # remove input SAM if it is a local file
   if deleteLocalFiles and len(os.path.dirname(bamFileIn)) == 0:
      os.remove(bamFileIn)
   
   # call Fulcrum Genomics consensus reads using fgbio CallMolecularConsensusReads
   cmd = javaExe + " -jar " \
   + " -Dsamjdk.use_async_io_write_samtools=true" \
   + " -Dsamjdk.use_async_io_read_samtools=true" \
   + " -Dsamjdk.compression_level=1 " \
   + fgbioJar + " CallMolecularConsensusReads" \
   + " -i " + readSet + ".consensus.temp0.bam" \
   + " -o " + readSet + ".consensus.temp1.bam" \
   + " --min-reads=2 --min-input-base-quality=25 -1 40 -2 40 -D -S unsorted " \
   + " --tag=" + tagNameUmi \
   + " > " + readSet + ".consensus.CallMolecularConsensusReads.log 2>&1"
   subprocess.check_call(cmd, shell=True)
   print("consensus: done with fgbio CallMolecularConsensusReads")
   os.remove(readSet + ".consensus.temp0.bam")
   
   # add @RG SM tag to header (bug in CallMolecularConsensusReads!!!)
   cmd = samtoolsDir + "samtools view -H " \
             + readSet + ".consensus.temp1.bam" \
   + " 2>> " + readSet + ".consensus.shell.log" \
   + " | awk '{if ($0 ~/^@RG/) print $0\"\tSM:foobar\"; else print $0;}'" \
   + " > "   + readSet + ".consensus.temp1.header.sam" \
   + " 2>> " + readSet + ".consensus.shell.log"
   subprocess.check_call(cmd, shell=True)
   
   # reheader the unaligned BAM using the corrected header file
   cmd = samtoolsDir + "samtools reheader -P " \
             + readSet + ".consensus.temp1.header.sam" \
   + " "     + readSet + ".consensus.temp1.bam" \
   + " > "   + readSet + ".consensus.temp2.bam" \
   + " 2>> " + readSet + ".consensus.shell.log" 
   subprocess.check_call(cmd, shell=True)
   os.remove(readSet + ".consensus.temp1.bam")

   # filter consensus reads using Fulcrum Genomics fgbio FilterConsensusReads
   cmd = javaExe + " -jar " \
   + fgbioJar + " FilterConsensusReads" \
   + " -i " + readSet + ".consensus.temp2.bam" \
   + " -o " + readSet + ".consensus.temp3.bam" \
   + " -r " + genomeFile  \
   + " --min-reads=2" \
   + " --max-base-error-rate=0.2" \
   + " --min-base-quality=0" \
   + " > " + readSet + ".consensus.FilterConsensusReads.log 2>&1"
   subprocess.check_call(cmd, shell=True)
   print("consensus: done with fgbio FilterConsensusReads")
   os.remove(readSet + ".consensus.temp2.bam")

   # copy header of unaligned BAM from fbbio FilterConsensusReads
   cmd = samtoolsDir + "samtools view -H " \
            + readSet + ".consensus.temp3.bam" \
   + " 1> " + readSet + ".consensus.temp4.sam" \
   + " 2>>" + readSet + ".consensus.shell.log" 
   subprocess.check_call(cmd, shell=True)
   
   # sort unaligned BAM by readId (which is also the molecule tag), using Linux sort, not samtools sort -n
   cmd = samtoolsDir + "samtools view " \
   + readSet + ".consensus.temp3.bam" \
   + " | sort -k1,1 -T./ --parallel=" + numCores \
   + " 1>> " + readSet + ".consensus.temp4.sam" \
   + " 2>> " + readSet + ".consensus.shell.log"
   subprocess.check_call(cmd, shell=True)
   os.remove(readSet + ".consensus.temp3.bam")
   print("consensus: done sorting unaligned BAM from FilterConsensusReads by readId")
   
   # convert unaligned SAM to FASTQ format
   cmd = samtoolsDir + "samtools fastq -n " \
            + readSet + ".consensus.temp4.sam" \
   + " -1 " + readSet + ".consensus.temp5.R1.fastq" \
   + " -2 " + readSet + ".consensus.temp5.R2.fastq" \
   + " >> " + readSet + ".consensus.shell.log 2>&1" 
   subprocess.check_call(cmd, shell=True)
   print("consensus: done converting unaligned BAM to FASTQ format")
   os.remove(readSet + ".consensus.temp4.sam")
   
   # sort primer info by molecule tag (which is also the read id)
   cmd = "sort -k1,1 -T./ --parallel=" + numCores \
   + "     " + primerFileIn \
   + " 1>  " + readSet + ".consensus.primers.txt" \
   + " 2>> " + readSet + ".consensus.shell.log"
   subprocess.check_call(cmd, shell=True)
   print("consensus: done sorting primer file by moleculeId/readId")

   # merge primer tag onto comment section of FASTQ read id line (to enable use of BWA MEM -C to propagate to aligned BAM file)
   for readEnd in ("R1","R2"):
   
      # open files
      fileInP = open(readSet + ".consensus.primers.txt","r")
      fileInR = open(readSet + ".consensus.temp5.{}.fastq".format(readEnd),"r")
      fileOut = open(readSet + ".consensus.{}.fastq".format(readEnd),"w")

      # loop over read file
      moleculeId_ = "foobar"
      lineNum = -1
      for line in fileInR:
      
         # fastq 4-line accounting
         lineNum += 1
         if lineNum == 4:
            lineNum = 0
            
         # echo other fastq lines
         if lineNum != 0:
            fileOut.write(line)
            continue
      
         # get molecule id
         readId = line.strip()
         readId = readId[1:]
         moleculeId = readId[1:]  
         
         # spin the primer file forward if necessary (because it might have more molecules than the reads file does)
         while moleculeId_ != moleculeId:
            line = fileInP.readline()
            if len(line) == 0:
               break
            moleculeId_, primerTag = line.strip().split("|")
         
         # crash if not found
         if moleculeId_ != moleculeId:
            raise Exception("consensus: could not find primer info for molecule: " + moleculeId)
         
         # add the primer info tag and write output
         fileOut.write("@{} pr:Z:{}\n".format(readId, primerTag))

      # done
      fileInP.close()
      fileInR.close()
      fileOut.close()
      os.remove(readSet + ".consensus.temp5.{}.fastq".format(readEnd))

#-------------------------------------------------------------------------
# filter genome alignments of consensus reads
#-------------------------------------------------------------------------
def filter(cfg,bamFileIn,bamFileOut):
   print("consensus filter: starting...")
   
   # get params
   deleteLocalFiles = cfg.deleteLocalFiles

   # constants for read pair accounting
   NUM_PRIMER_SIDE_NOT_MAPPED = 0
   NUM_RANDOM_SIDE_NOT_MAPPED = 1
   NUM_R1_R2_NOT_AT_SAME_LOCUS = 2
   NUM_R1_R2_SAME_ORIENTATION = 3
   NUM_SPLIT_ALIGNMENT = 4
   NUM_LOW_MAPQ = 5
   NUM_LT_25BP_ALIGNED = 6
   NUM_WRITTEN_OUT = 7
   NUM_METRICS_TOTAL = 8
   
   # open BAM read alignment files
   bamIn  = pysam.Samfile(bamFileIn , "rb")
   bamOut = pysam.Samfile(bamFileOut, "wb", template=bamIn)

   # loop over read alignments
   readPairCounts = [0] * NUM_METRICS_TOTAL
   for read in bamIn:
   
      # this is dangerous, but drop these for now
      if read.is_secondary or read.is_supplementary:
         continue

      # crash if read is not paired
      if not read.is_paired:
         print(read.qname)
         raise Exception("read not paired!")
      
      # this should be R1      
      read1 = read
      
      # get mate, assuming mate is the next record in the BAM file
      while True:
         read = bamIn.next()
         if not read.is_secondary and not read.is_supplementary:
            break
            
      # this should be R2
      read2 = read
      
      # debug check
      if read1.qname != read2.qname:
         print(read1.qname, read2.qname)
         raise Exception("read mate is not next in BAM record order!")
         
      # debug check 
      if not read1.is_read1 or not read2.is_read2:
         raise Exception("R1/R2 mixed up!")
     
      # skip but count unmapped R1 reads, even if R2 mapped.  Need to look at these later...
      if read1.is_unmapped:
         readPairCounts[NUM_PRIMER_SIDE_NOT_MAPPED] += 1
         continue
      
      # skip but count unmapped R2 reads
      if read2.is_unmapped:
         readPairCounts[NUM_RANDOM_SIDE_NOT_MAPPED] += 1
         continue
         
      # skip reads not mapped to same chrom
      chrom1 = bamIn.getrname(read1.tid)
      chrom2 = bamIn.getrname(read2.tid)
      if chrom1 != chrom2:
         readPairCounts[NUM_R1_R2_NOT_AT_SAME_LOCUS] += 1
         continue
         
      # skip reads not mapped to same locus
      locRead1 = int(read1.aend) - 1 if read1.is_reverse else read1.pos
      locRead2 = int(read2.aend) - 1 if read2.is_reverse else read2.pos
      if abs(locRead1 - locRead2) > 2000:
         readPairCounts[NUM_R1_R2_NOT_AT_SAME_LOCUS] += 1
         continue

      # skip pairs with odd alignment orientation
      if read1.is_reverse == read2.is_reverse:
         readPairCounts[NUM_R1_R2_SAME_ORIENTATION] += 1
         continue
      
      # drop read pair if either end has a supplementary split alignment
      if read1.has_tag("SA") or read2.has_tag("SA"):
         readPairCounts[NUM_SPLIT_ALIGNMENT] += 1
         continue
            
      # drop read pair if R1 or R2 read has low mapq 
      if read2.mapq < 17 or read1.mapq < 17:
         readPairCounts[NUM_LOW_MAPQ] += 1
         continue
         
      # require some significant alignment to genome
      if read2.aend - read2.pos < 25 or read1.aend - read1.pos < 25:
         readPairCounts[NUM_LT_25BP_ALIGNED] += 1
         continue
      
      # output
      bamOut.write(read1)
      bamOut.write(read2)
      readPairCounts[NUM_WRITTEN_OUT] += 1
      
   # done
   bamOut.close()
   bamIn.close()
   
   # delete input BAM file if local
   if deleteLocalFiles and len(os.path.dirname(bamFileIn)) == 0:
      os.remove(bamFileIn)
   
   # report drop totals
   print("{} read fragments dropped, primer side read not mapped".format(readPairCounts[NUM_PRIMER_SIDE_NOT_MAPPED]))
   print("{} read fragments dropped, random side read not mapped".format(readPairCounts[NUM_RANDOM_SIDE_NOT_MAPPED]))
   print("{} read fragments dropped, R1 and R2 not mapped to same locus".format(readPairCounts[NUM_R1_R2_NOT_AT_SAME_LOCUS]))
   print("{} read fragments dropped, FF or RR mapping orientation".format(readPairCounts[NUM_R1_R2_SAME_ORIENTATION]))
   print("{} read fragments dropped, split alignment".format(readPairCounts[NUM_SPLIT_ALIGNMENT]))
   print("{} read fragments dropped, low mapping quality MAPQ < 17".format(readPairCounts[NUM_LOW_MAPQ]))
   print("{} read fragments dropped, less than 25 bp aligned to genome".format(readPairCounts[NUM_LT_25BP_ALIGNED]))
   print("{} read fragments written".format(readPairCounts[NUM_WRITTEN_OUT]))
