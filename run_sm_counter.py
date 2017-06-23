import sys
sys.path.append("/srv/qgen/code/qiaseq-dna/")

# our modules
import core.run_log
import core.run_config
import core.prep
import core.align
import core.umi_filter
import core.umi_mark
import core.umi_merge
import core.primer_clip
import varcall.sm_counter_wrapper
import varcall.vcf_complex
import varcall.vcf_annotate

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
   
   # soft clip primer regions from read alignments
   bamFileIn  = readSet + ".umi_merge.bam"
   bamFileOut = readSet + ".primer_clip.bam"
   core.primer_clip.run(cfg,bamFileIn,bamFileOut,False)

   # sort the final BAM file, to prepare for downstream variant calling
   bamFileIn  = readSet + ".primer_clip.bam"
   bamFileOut = readSet + ".bam"
   core.samtools.sort(cfg,bamFileIn,bamFileOut)

   # run smCounter variant calling
   numVariants = varcall.sm_counter_wrapper.run(cfg, paramFile)
   
   # create complex variants, and annotate using snpEff
   if numVariants > 0:
   
      # convert nearby primitive variants to complex variants
      bamFileIn  = readSet + ".bam"
      vcfFileIn  = readSet + ".smCounter.cut.vcf"
      vcfFileOut = readSet + ".smCounter.cplx.vcf"
      varcall.vcf_complex.run(cfg, bamFileIn, vcfFileIn, vcfFileOut)
         
      # annotate variants in the VCF file
      vcfFileIn  = readSet + ".smCounter.cplx.vcf"
      vcfFileOut = readSet + ".smCounter.anno.vcf"
      varcall.vcf_annotate.run(cfg, vcfFileIn, vcfFileOut)
         
   # close log file
   core.run_log.close()
   
#-------------------------------------------------------------------------------------
# main program for running from shell 
#-------------------------------------------------------------------------------------
if __name__ == "__main__":
   readSet   = sys.argv[1]
   paramFile = sys.argv[2] if len(sys.argv) == 3 else "run-params.txt"
   run(readSet, paramFile)
