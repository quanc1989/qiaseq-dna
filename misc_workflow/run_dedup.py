import sys

# our modules
sys.path.append(os.path.abspath(os.path.dirname(os.path.dirname(__file__))))
import core.run_log
import core.run_config
import core.prep
import core.align
import core.dedup

#--------------------------------------------------------------------------------------
# trim common regions, align to genome, remove PCR duplicates
#--------------------------------------------------------------------------------------
def run(readSet, paramFile):

   # initialize logger
   core.run_log.init(readSet)

   # read run configuration file to memory
   cfg = core.run_config.run(readSet,paramFile)

   # trim 3' ends of both reads, and extract UMI sequence
   core.prep.run(cfg)
   
   # align trimmed reads to genome using BWA MEM
   readFileIn1 = readSet + ".prep.R1.fastq"
   readFileIn2 = readSet + ".prep.R2.fastq"
   bamFileOut  = readSet + ".align.bam"
   core.align.run(cfg, readFileIn1, readFileIn2, bamFileOut)
   
   # run Picard MarkDuplicates
   bamFileIn  = readSet + ".align.bam"
   bamFileOut = readSet + ".dedup.bam"
   core.dedup.run(cfg,bamFileIn,bamFileOut)
   
   # close log file
   core.run_log.close()
   
#-------------------------------------------------------------------------------------
# main program for running from shell 
#-------------------------------------------------------------------------------------
if __name__ == "__main__":
   readSet   = sys.argv[1]
   paramFile = sys.argv[2] if len(sys.argv) == 3 else "run-params.txt"
   run(readSet, paramFile)
