import ConfigParser
import sys
import multiprocessing
# our modules
import core.run_log
import core.run_config
import core.prep
import core.align
import core.umi_filter
import core.umi_mark
import core.umi_merge
import core.primer_clip
import core.samtools
import core.tumor_normal
import core.sm_counter_wrapper
import metrics.sum_specificity
import metrics.sum_uniformity_primer
import metrics.sum_primer_umis
import metrics.sum_all
import metrics.umi_frags
import metrics.umi_depths
import misc.process_ion
import misc.tvc
import annotate.vcf_complex
import annotate.vcf_annotate

#--------------------------------------------------------------------------------------
# call input molecules, build consenus reads, align to genome, trim primer region
#--------------------------------------------------------------------------------------
def run(args):
   readSet, paramFile, vc = args
   # initialize logger
   core.run_log.init(readSet)

   # read run configuration file to memory
   cfg = core.run_config.run(readSet,paramFile)

   if cfg.platform.lower() == "illumina":

      if cfg.duplex.lower() == "true": ## Duplex sequencing run
         pass

      # trim 3' ends of both reads, and extract UMI sequence
      core.prep.run(cfg)
   
      # align trimmed reads to genome using BWA MEM
      readFileIn1 = readSet + ".prep.R1.fastq"
      readFileIn2 = readSet + ".prep.R2.fastq"
      bamFileOut  = readSet + ".align.bam"
      core.align.run(cfg, readFileIn1, readFileIn2, bamFileOut)
   else:
      misc.process_ion.trimIon(cfg)
      misc.process_ion.alignToGenomeIon(cfg)

   # call putative unique input molecules using BOTH UMI seq AND genome alignment position on random fragmentation side
   bamFileIn  = readSet + ".align.bam"
   core.umi_filter.run(cfg, bamFileIn)
   core.umi_mark.run(cfg)   
   metrics.umi_frags.run(cfg)   
   metrics.umi_depths.run(cfg)   
   core.umi_merge.run(cfg, bamFileIn)
   
   # soft clip primer regions from read alignments
   bamFileIn  = readSet + ".umi_merge.bam"
   bamFileOut = readSet + ".primer_clip.bam"
   core.primer_clip.run(cfg,bamFileIn,bamFileOut,False)

   # additional metrics to generate
   metrics.sum_primer_umis.run(cfg) # primer-level umi and read metrics
   metrics.sum_specificity.run(cfg) # priming specificity
   metrics.sum_uniformity_primer.run(cfg) # primer-level uniformity

   # sort the final BAM file, to prepare for downstream variant calling
   bamFileIn  = readSet + ".primer_clip.bam"
   bamFileOut = readSet + ".bam"
   core.samtools.sort(cfg,bamFileIn,bamFileOut)   
  
   if cfg.duplex.lower() == "false": # do not run smCounter for duplex reads

      if cfg.platform.lower() != "illumina": # ion reads
         misc.tvc.run(cfg)

      # run smCounter variant calling
      numVariants = core.sm_counter_wrapper.run(cfg, paramFile, vc)
      
      if cfg.platform.lower() != "illumina":
         numVariants = misc.tvc.smCounterFilter(cfg,vc)

      # create complex variants, and annotate using snpEff
      if numVariants > 0:
         # convert nearby primitive variants to complex variants
         bamFileIn  = readSet + ".bam"
         vcfFileIn  = readSet + ".smCounter.cut.vcf"
         vcfFileOut = readSet + ".smCounter.cplx.vcf"
         annotate.vcf_complex.run(cfg, bamFileIn, vcfFileIn, vcfFileOut)
         
         # annotate variants in the VCF file
         vcfFileIn  = readSet + ".smCounter.cplx.vcf"
         vcfFileOut = readSet + ".smCounter.anno.vcf"
         annotate.vcf_annotate.run(cfg, vcfFileIn, vcfFileOut,vc)        
   
   # aggregate all summary metrics
   metrics.sum_all.run(cfg)
   
   # close log file
   core.run_log.close()


def run_tumor_normal(readSet,paramFile,vc):
   ''' Wrapper around run() for tumor-normal analysis
   '''
   # 2 read set names which are space delimited
   readSets = filter(None,readSet.split(" "))
   assert len(readSets) == 2, "Tumor-Normal Analysis requires exactly 2 read sets !"

   # read parameter file
   parser = ConfigParser.SafeConfigParser()
   parser.optionxform = str
   parser.read(paramFile)

   tumor = None
   normal = None
   for section in parser.sections():
      if section not in ['general','smCounter']:
         for (paramName, paramVal) in parser.items(section):
            if paramName == 'sampleType' and paramVal.lower() == 'normal':
               normal = section
            elif paramName == 'sampleType' and paramVal.lower() == 'tumor':
               tumor = section

   assert tumor!=None and normal!=None, "Could not sync read set names supplied with config file !"
   
   p = multiprocessing.Pool(2)
   p.map(run,[(normal,paramFile,vc),(tumor,paramFile,vc)])
   p.close()
   p.join()
   ## Additional analysis steps
   cfg = core.run_config.run(tumor,paramFile)
   core.tumor_normal.removeNormalVariants(cfg)
   core.tumor_normal.runCopyNumberEstimates(cfg)

#-------------------------------------------------------------------------------------
# main program for running from shell 
#-------------------------------------------------------------------------------------
if __name__ == "__main__":

   if len(sys.argv) != 5 :
      print "\nRun as : python run_qiaseq_dna.py <param_file> <v1/v2> <single/tumor-normal> <readSet>\n"
      sys.exit(-1)

   paramFile = sys.argv[1]
   vc = sys.argv[2]
   analysis = sys.argv[3]
   readSet   = " ".join(sys.argv[4:]) # 2 readSets in case of tumor-normal

   if analysis.lower() == "tumor-normal":      
      run_tumor_normal(readSet,paramFile,vc)
   else: # Single sample, might still need to run quandico
      run((readSet,paramFile,vc))
      cfg = core.run_config.run(readSet,paramFile)
      tumor_normal.runCopyNumberEstimates(cfg)
