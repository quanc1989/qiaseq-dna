import glob
import gzip
import os
import os.path
import math
import subprocess
from multiprocessing.dummy import Pool as ThreadPool

#-------------------------------------------------------------------------------------
def worker(cmd):
   try:
      subprocess.check_call(cmd, shell=True)
      return True
   except:
      return False

#-------------------------------------------------------------------------------------
def splitReadFile(readFile,filePrefixOut,readSide,numBatchesMax,deleteLocalFiles):

   # validate existence of input read file
   if not os.path.isfile(readFile):
      raise UserWarning("input read file does not exist: " + readFile)
         
   # get number of reads in the fastq file
   isGzFile = readFile.endswith(".fastq.gz")
   if isGzFile:
      line = subprocess.check_output("zcat {} | wc -l".format(readFile), shell=True)
      numLines = int(line)
   else:
      line = subprocess.check_output("wc -l " + readFile, shell=True)
      vals = line.split(" ")
      numLines = int(vals[0])
   if numLines % 4 != 0 or numLines == 0:
      raise UserWarning("truncated or empty fastq read file {}".format(readFile))
   numReads = numLines / 4
   
   # set batch size
   batchSize = 4 * int(math.ceil(1.00 * numReads / numBatchesMax))
   
   # open input fastq read file
   if isGzFile:
      fileIn = gzip.open(readFile,"rb")
   else:
      fileIn = open(readFile,"r")

   # write fastq batches to disk - one batch to be done on each CPU core
   numLinesOut = 0
   batchNum = 0
   fileOut = None
   for line in fileIn:
      # open new file if needed
      if numLinesOut == 0:
         fileOut = open("{}.{:04d}.{}.fastq".format(filePrefixOut,batchNum,readSide), "w")
         
      # write fastq line to disk
      fileOut.write(line)
      numLinesOut += 1
         
      # close output file if done with batch
      if numLinesOut == batchSize:
         fileOut.close()
         numLinesOut = 0
         batchNum += 1

   # close out last batch
   if numLinesOut != 0:
      fileOut.close()
      batchNum += 1
   numBatches = batchNum
   
   # delete local input file if no longer needed
   if deleteLocalFiles and len(os.path.dirname(readFile)) == 0:
      os.remove(readFile)
   
   # done
   return numReads, numBatches

#-------------------------------------------------------------------------------------
def run(cfg):
   # report start
   print("prep: starting read prep - trimming ends and UMI extraction")
   
   # get params
   readSet          = cfg.readSet
   readFile1        = cfg.readFile1
   readFile2        = cfg.readFile2
   trimScript       = cfg.trimScript
   cutadaptDir      = cfg.cutadaptDir
   deleteLocalFiles = cfg.deleteLocalFiles
   numCores         = int(cfg.numCores)
   tagNameUmiSeq    = cfg.tagNameUmiSeq
   
   # set output file prefix
   filePrefixOut = readSet + ".prep"
   
   # check for adapter on primer side - customer forgot to use the custom sequencing primer, or are a small part of a large HiSeq run
   if readFile1.endswith(".fastq"):  
      cmd = "head -n 4000 " + readFile1 + " > " 
   else:
      cmd = "zcat " + readFile1 + " | head -n 4000 > " 
   cmd += filePrefixOut + ".temp1000.R1.fastq"
   subprocess.check_call(cmd, shell=True)
   cmd = cutadaptDir + "cutadapt -e 0.18 -O 18" \
       + " -g ^AATGTACAGTATTGCGTTTTG -n 1" \
       + " -o /dev/null " \
       +         filePrefixOut + ".temp1000.R1.fastq" \
       + " > " + filePrefixOut + ".cutadapt.5.R1.log 2>&1"
   subprocess.check_call(cmd, shell=True)
   os.remove(filePrefixOut + ".temp1000.R1.fastq")
   pctWithAdapter = None
   for line in open(filePrefixOut + ".cutadapt.5.R1.log", "r"):
      if line.startswith("Reads with adapters:"):
         idx1 = line.find("(")
         idx2 = line.find("%)")
         pctWithAdapter = float(line[idx1+1:idx2])
         break
   print("prep: check for wrong sequencing primer using first 1,000 reads: % wrong is {}".format(pctWithAdapter))
   if pctWithAdapter > 95.0:
      print("WARNING: R1 reads start with PCR adapter - custom sequencing primer was not used on primer side!")
      cmd = cutadaptDir + "cutadapt -e 0.18 -O 18" \
          + " -g ^AATGTACAGTATTGCGTTTTG -n 1" \
          + " -o " + filePrefixOut + ".fixed.R1.fastq " \
                   + readFile1  \
          + " > "  + filePrefixOut + ".cutadapt.5.R1.log 2>&1"
      subprocess.check_call(cmd, shell=True)
      readFile1 = filePrefixOut + ".fixed.R1.fastq"
   
   # split both read files into chunks to be processed in parallel
   numReads1, numBatches = splitReadFile(readFile1,filePrefixOut,"R1",numCores,deleteLocalFiles)
   numReads2, numBatches = splitReadFile(readFile2,filePrefixOut,"R2",numCores,deleteLocalFiles)
   
   # debug check
   if numReads2 != numReads1:
      raise UserWarning("prep: R1 and R2 read count not equal!")
      
   # debug check
   if numReads1 == 0:
      raise UserWarning("prep: input read files are empty!")

   # set up trimming work to be run in parallel sub-processes, using another python script
   workIn = []
   for batchNum in range(numBatches):
      filePrefixBatch = "{}.{:04d}".format(filePrefixOut,batchNum)
      cmd = "python {0} {1} {2} {3} > {3}.log 2>&1 ".format(trimScript,cutadaptDir,tagNameUmiSeq,filePrefixBatch)
      workIn.append(cmd)
      
   # run cutadapt and UMI extraction in parallel sub-processes
   print("prep: starting parallel trimming batches")
   pool = ThreadPool(min(numCores,len(workIn)))
   workOut = pool.map(worker, workIn)
   pool.close()
   pool.join()
   print("prep: completed parallel trimming batches")
   
   # make sure all batches of work completed successfully
   for batchNum in range(len(workOut)):
      if not workOut[batchNum]:
         raise Exception("read trimming failed for batch: {:04d}".format(batchNum))

   # concatenate the read files back into one file
   for readEnd in ("R1","R2"):
   
      # delete output file if it exists
      readFileOut = "{}.{}.fastq".format(filePrefixOut,readEnd)
      if os.path.isfile(readFileOut):
         os.remove(readFileOut)
   
      # concatenate read file and delete (Linux cat probaby faster than Python line reads)
      for batchNum in range(numBatches):
         readFileIn = "{}.{:04d}.{}.fastq".format(filePrefixOut, batchNum, readEnd)
         cmd = "cat {} >> {} ".format(readFileIn,readFileOut)
         subprocess.check_call(cmd,shell=True)
         os.remove(readFileIn)

   # concatenate the log files - note dangerous wildcards here!
   logFileOut = open(filePrefixOut + ".log","w")
   for logFileIn in glob.glob(filePrefixOut + ".0*.log"):
      for line in open(logFileIn,"r"):
        logFileOut.write(line)
      os.remove(logFileIn)
      
   # aggregate summary read count files - for some trim scripts these files contain important read count metrics
   output = []
   firstFile = True
   for sumFileIn in glob.glob(filePrefixOut + ".0*.summary.txt"):
      iLine = 0
      for line in open(sumFileIn,"r"):
         metricVal,metricName = line.strip().split("\t")
         metricVal = int(metricVal)
         if firstFile:
            output.append([metricVal,metricName])
         else:
            output[iLine][0] += metricVal
         iLine += 1
      firstFile = False
      os.remove(sumFileIn)
   sumFileOut = open(filePrefixOut + ".summary.txt","w")
   for row in output:
      sumFileOut.write("\t".join((str(x) for x in row)))
      sumFileOut.write("\n")
   sumFileOut.close()
   
   # report completion
   print("prep: done")
   
#-------------------------------------------------------------------------------------
if __name__ == "__main__":
   cfg = lambda:0
   cfg.readSet = "NEB_S2"
   cfg.readFile1 = "/mnt/webserver/datadisk/resources/jdicarlo/NEB_S2_L001_R1_001.fastq.gz"
   cfg.readFile2 = "/mnt/webserver/datadisk/resources/jdicarlo/NEB_S2_L001_R2_001.fastq.gz"
   cfg.cutadaptDir = "/srv/qgen/bin/cutadapt-1.10/"
   cfg.trimScript  = "prep_trim.py"
   cfg.tagNameUmiSeq = "mi"
   cfg.numCores  = "32"
   cfg.deleteLocalFiles = True
   run(cfg)
