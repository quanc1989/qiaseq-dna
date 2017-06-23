import subprocess
import os
import os.path

#-------------------------------------------------------------------------
# align reads to refernce genome using BWA MEM
#-------------------------------------------------------------------------
def run(cfg, readFile1, readFile2, bamFileOut):
   # report start
   print("align: starting read alignment")

   # get some parameters from config
   readSet          = cfg.readSet
   genomeFile       = cfg.genomeFile
   bwaDir           = cfg.bwaDir
   samtoolsDir      = cfg.samtoolsDir
   deleteLocalFiles = cfg.deleteLocalFiles
   numCores         = cfg.numCores
   
   # make file names for local BWA and samtools log files
   logFileBase     = os.path.basename(bamFileOut)
   logFileBwa      = logFileBase.replace(".bam",".bwa.log")
   logFileSamtools = logFileBase.replace(".bam",".samtools.log")

   # align reads to reference genome using BWA-MEM, and convert to BAM format
   cmd = bwaDir + "bwa mem -C -t " +  numCores \
   + " " + genomeFile   \
   + " " + readFile1    \
   + " " + readFile2    \
   + " 2>" + logFileBwa \
   + " | " + samtoolsDir + "samtools view -1 -@ " + numCores \
   + " - " \
   + " 1> " + bamFileOut \
   + " 2> " + logFileSamtools
   subprocess.check_call(cmd, shell=True)

   # delete local fastq files if not needed anymore, to conserve local disk space
   if deleteLocalFiles:
      for readFile in (readFile1, readFile2):
         if len(os.path.dirname(readFile)) == 0:
            os.remove(readFile)
   
   # report completion
   print("align: done creating bam file " + bamFileOut)
   
#-------------------------------------------------------------------------------------
if __name__ == "__main__":
   cfg = lambda:0
   cfg.readSet = "NEB_S2"
   cfg.genomeFile  = "/srv/qgen/data/genome/ucsc.hg19.fa"
   cfg.bwaDir      = "/srv/qgen/bin/bwa-0.7.15/"
   cfg.samtoolsDir = "/srv/qgen/bin/samtools-1.3.1/"
   cfg.deleteLocalFiles = False
   cfg.numCores = "32"
   run(cfg, "prep", "align")
