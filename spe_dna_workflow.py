import luigi
import sys
import os

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
import metrics.umi_frags
import metrics.umi_depths
import varcall.sm_counter_wrapper
import varcall.vcf_complex
import varcall.vcf_annotate

class config(luigi.Config):
    ''' Map config file variables into class attributes
    '''
    # tools
    cutadaptDir = luigi.Parameter()
    bwaDir = luigi.Parameter()
    samtoolsDir = luigi.Parameter()
    javaExe = luigi.Parameter()
    sswPyFile = luigi.Parameter()
    # general params
    numCores = luigi.Parameter()
    deleteLocalFiles = luigi.BoolParameter()
    samtoolsMem = luigi.Parameter()
    outputDetail = luigi.BoolParameter()
    genomeFile = luigi.Parameter()
    # umi module
    endogenousLenMin = luigi.Parameter()
    # SAM tag names
    tagNameUmiSeq = luigi.Parameter()
    tagNameUmi = luigi.Parameter()
    tagNamePrimer = luigi.Parameter()
    tagNameResample = luigi.Parameter()
    # varaint primitive to complex conversion
    vcfComplexGapMax = luigi.IntParameter()
    # variant annotation
    snpEffPath = luigi.Parameter()
    snpEffConfig = luigi.Parameter()
    dbSnpFile = luigi.Parameter()
    cosmicFile   = luigi.Parameter()
    clinVarFile  = luigi.Parameter()

class smcounter(luigi.Config):
    # smCounter params
    minBQ = luigi.IntParameter()
    minMQ = luigi.IntParameter()
    hpLen = luigi.IntParameter()
    mismatchThr = luigi.FloatParameter()
    mtDrop = luigi.IntParameter()
    maxMT = luigi.IntParameter()
    primerDist = luigi.IntParameter()
    threshold = luigi.IntParameter()
    bedtoolsPath = luigi.Parameter()
    bedTandemRepeats = luigi.Parameter()
    bedRepeatMaskerSubset = luigi.Parameter()

#--------------------------------------------------------------------------------------
# Create an aggregate object to store config values as used in the previous pipeline
#--------------------------------------------------------------------------------------
cfg = config()
for param,val in smcounter().__dict__['param_kwargs'].items():
    setattr(cfg,param,val)
    
class TrimReads(luigi.Task):
    ''' Trim 3' ends of both reads, and extract UMI sequence
    '''

class AlignReads(luigi.Task):
    ''' Align trimmed reads to genome using bwa mem
    '''

class AssignUMI(luigi.Task):
    ''' Call putative unique input molecules using BOTH UMI seq AND genome alignment position on random fragmentation side
    '''

class SoftClipBam(luigi.Task):
    ''' Soft clip primer regions from read alignments
    '''

class SortBam(luigi.Task):
    ''' Sort the final BAM file, to prepare for downstream variant calling
    '''

class RunSmCounter(luigi.Task):
    ''' Run smCounter variant calling
    '''
    
class AnnotateVCF(luigi.Task):
    ''' Annotate a VCF file
    '''
    def __init__(self,*args,**kwargs):
        '''
        '''
        super(AnnotateVCF,self).__init__(*args,**kwargs)
        ## Add readSet info to global config
        global cfg
        parser = ConfigParser.ConfigParser()
        parser.read(self.samples_cfg)        
        ## Handle only 1 readSet at a time for now
        found=False
        for section in parser.sections():
            if section == self.readSet:
                setattr(cfg,"readSet",section)
                setattr(cfg,"readFile1",parser.get(section,"readFile1"))
                setattr(cfg,"readFile2",parser.get(section,"readFile2"))
                setattr(cfg,"primerFile",parser.get(section,"primerFile"))
                setattr(cfg,"roiBedFile",parser.get(section,"roiBedFile"))
                found=True

        assert found==True,"Pipeline Failire !. Could not find readSet name in your config file !"

    def requires(self):
        '''
        '''       
        return self.clone(RunSmCounter)
            
    def run(self):
        '''
        '''
        # convert nearby primitive variants to complex variants
        bamFileIn = cfg.readSet + ".bam"
        vcfFileIn = cfg.readSet + ".smCounter.cut.vcf"
        vcfFileOut = cfg.readSet + ".smCounter.cplx.vcf"
        varcall.vcf_complex.run(cfg,bamFileIn,vcfFileIn,vcfFileOut)
        # annotate variants in the VCF file
        vcfFileIn = cfg.readSet + ".smCounter.cplx.vcf"
        vcfFileOut = cfg.readSet + ".smCounter.anno.vcf"
        varcall.vcf_annotate.run(cfg,vcfFileIn,vcfFileOut)
        # write the verification file
        with open(self.verification_file,'w') as OUT:
            OUT.write("done\n")

    def output(self):
        '''
        '''
        return luigi.LocalTarget(self.verification_file)
        
def main():
    '''
    '''   

    
if __name__ == '__main__':
    main()
