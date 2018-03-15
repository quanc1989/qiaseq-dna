import sys

# our modules
sys.path.append(os.path.abspath(os.path.dirname(os.path.dirname(__file__))))
import core.run_log
import core.run_config
import core.prep
import core.align
import core.umi_filter
import core.umi_mark
import core.umi_merge
import core.primer_clip
import core.consensus
import core.samtools

#--------------------------------------------------------------------------------------
# call input molecules, build consenus reads, align to genome, trim primer region
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
   
   # call putative unique input molecules using BOTH UMI seq AND genome alignment position on random fragmentation side
   bamFileIn  = readSet + ".align.bam"
   core.umi_filter.run(cfg, bamFileIn)
   core.umi_mark.run(cfg)
   core.umi_merge.run(cfg, bamFileIn)

   # soft clip primer regions, but only on reads from internal resampling priming
   bamFileIn  = readSet + ".umi_merge.bam"
   bamFileOut = readSet + ".primer_clip.bam"
   core.primer_clip.run(cfg,bamFileIn,bamFileOut,True)

   # make consensus reads using Fulcrum Genomics fgbio tools, and merge primer info to fastq header comment
   bamFileIn    = readSet + ".primer_clip.bam"
   primerFileIn = readSet + ".umi_merge.primers.txt"
   core.consensus.run(cfg,bamFileIn,primerFileIn)
   
   # align consensus reads to genome using BWA MEM
   readFileIn1 = readSet + ".consensus.R1.fastq"
   readFileIn2 = readSet + ".consensus.R2.fastq"
   bamFileOut  = readSet + ".consensus.align.bam"
   core.align.run(cfg, readFileIn1, readFileIn2, bamFileOut)
   
   # filter consensus read alignments that might cause trouble for SNP/indel calling and/or primer region clipping
   bamFileIn  = readSet + ".consensus.align.bam"
   bamFileOut = readSet + ".consensus.filter.bam"
   core.consensus.filter(cfg,bamFileIn,bamFileOut)
   
   # soft clip primer regions from consensus reads
   bamFileIn  = readSet + ".consensus.filter.bam"
   bamFileOut = readSet + ".consensus.primer_clip.bam"
   core.primer_clip.run(cfg,bamFileIn,bamFileOut,False)

   # sort the final BAM file, to prepare for downstream variant calling
   bamFileIn  = readSet + ".consensus.primer_clip.bam"
   bamFileOut = readSet + ".consensus.sorted.bam"
   core.samtools.sort(cfg,bamFileIn,bamFileOut)

   # close log file
   core.run_log.close()
   
#-------------------------------------------------------------------------------------
# main program for running from shell 
#-------------------------------------------------------------------------------------
if __name__ == "__main__":
   readSet   = sys.argv[1]
   paramFile = sys.argv[2] if len(sys.argv) == 3 else "run-params.txt"
   run(readSet, paramFile)
