import sys
sys.path.append("/srv/qgen/code/qiaseq-dna/")

# our modules
import run_log
import run_config
import prep
import align

#--------------------------------------------------------------------------------------
# trim common regions, align to genome, remove PCR duplicates
#--------------------------------------------------------------------------------------
def run(readSet, paramFile):

   # initialize logger
   run_log.init(readSet)

   # read run configuration file to memory
   cfg = run_config.run(readSet,paramFile)

   # trim 3' ends of both reads, and extract UMI sequence
   prep.run(cfg)
   
   # align trimmed reads to genome using BWA MEM
   readFileIn1 = readSet + ".prep.R1.fastq"
   readFileIn2 = readSet + ".prep.R2.fastq"
   bamFileOut  = readSet + ".align.bam"
   align.run(cfg, readFileIn1, readFileIn2, bamFileOut)
   
   # close log file
   run_log.close()
   
#-------------------------------------------------------------------------------------
# main program for running from shell 
#-------------------------------------------------------------------------------------
if __name__ == "__main__":
   readSet   = sys.argv[1]
   paramFile = sys.argv[2] if len(sys.argv) == 3 else "run-params.txt"
   run(readSet, paramFile)
