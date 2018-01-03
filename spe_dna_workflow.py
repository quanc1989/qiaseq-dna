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


## To Do:
## 1. Add param for output directory

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

class MyExtTask(luigi.ExternalTask):
    ''' Checks whether the file specified exists on disk
    '''
    file_loc = luigi.Parameter()
    def output(self):
        return luigi.LocalTarget(self.file_loc)


class TrimReads(luigi.Task):
    ''' Trim 3' ends of both reads, and extract UMI sequence
    '''
    # Parameters
    samplesCfg = luigi.Parameter(description="Config file containing information on readSet to run")
    readSet = luigi.Parameter(description="The readSet to analyze from the config file above")
    
    def __init__(self,*args,**kwargs):
        '''
        '''
        super(self.__class__.__name__,self).__init__(*args,**kwargs)
        currentDir = os.path.dirname(os.path.realpath(__file__))
        verificationDir = os.path.join(currentDir,"verification")
        self.verificationFile = os.path.join(verificationDir,
                                             self.__class__.__name__+'.verification.txt')
        
    def requires(self):
        '''
        '''
        return MyExtTask(cfg.readFile1)
    
    def run(self):
        '''
        '''
        core.prep.run(cfg)
        with open(self.verificationFile,'w') as OUT:
            OUT.write("done\n")
            
    def output(self):
        '''
        '''
        return luigi.LocalTarget(self.verificationFile)
        
class AlignReads(luigi.Task):
    ''' Align trimmed reads to genome using bwa mem
    '''
    # Parameters
    samplesCfg = luigi.Parameter(description="Config file containing information on readSet to run")
    readSet = luigi.Parameter(description="The readSet to analyze from the config file above")
    
    def __init__(self,*args,**kwargs):
        '''
        '''
        super(self.__class__.__name__,self).__init__(*args,**kwargs)
        currentDir = os.path.dirname(os.path.realpath(__file__))
        verificationDir = os.path.join(currentDir,"verification")
        self.verificationFile = os.path.join(verificationDir,
                                             self.__class__.__name__+'.verification.txt')
        
    def requires(self):
        '''
        '''
        return self.clone(TrimReads)
    
    def run(self):
        '''
        '''
        readFileIn1 = cfg.readSet + ".prep.R1.fastq"
        readFileIn2 = cfg.readSet + ".prep.R2.fastq"
        bamFileOut = cfg.readSet + ".align.bam"
        core.align.run(cfg,readFileIn1,readFileIn2,bamFileOut)
        with open(self.verificationFile,'w') as OUT:
            OUT.write("done\n")
        
    def output(self):
        '''
        '''
        return luigi.LocalTarget(self.verificationFile)

class AssignUMI(luigi.Task):
    ''' Call putative unique input molecules using BOTH UMI seq AND genome alignment position on random fragmentation side
    '''
    # Parameters
    samplesCfg = luigi.Parameter(description="Config file containing information on readSet to run")
    readSet = luigi.Parameter(description="The readSet to analyze from the config file above")    

    def __init__(self,*args,**kwargs):
        '''
        '''
        super(self.__class__.__name__,self).__init__(*args,**kwargs)
        currentDir = os.path.dirname(os.path.realpath(__file__))
        verificationDir = os.path.join(currentDir,"verification")
        self.verificationFile = os.path.join(verificationDir,
                                             self.__class__.__name__+'.verification.txt')

    def requires(self):
        '''
        '''
        return self.clone(AlignReads)
    
    def run(self):
        '''
        '''
        bamFileIn = cfg.readSet + ".align.bam"
        core.umi_filter.run(cfg, bamFileIn)
        core.umi_mark.run(cfg)
        metrics.umi_frags.run(cfg)
        metrics.umi_depths.run(cfg)
        core.umi_merge.run(cfg, bamFileIn)
        with open(self.verificationFile,'w') as OUT:
            OUT.write("done\n")

    def output(self):
        '''
        '''
        return luigi.LocalTarget(self.verificationFile)        
    
class SoftClipBam(luigi.Task):
    ''' Soft clip primer regions from read alignments
    '''
    # Parameters
    samplesCfg = luigi.Parameter(description="Config file containing information on readSet to run")
    readSet = luigi.Parameter(description="The readSet to analyze from the config file above")    

    def __init__(self,*args,**kwargs):
        '''
        '''
        super(self.__class__.__name__,self).__init__(*args,**kwargs)
        currentDir = os.path.dirname(os.path.realpath(__file__))
        verificationDir = os.path.join(currentDir,"verification")
        self.verificationFile = os.path.join(verificationDir,
                                             self.__class__.__name__+'.verification.txt')

    def requires(self):
        '''
        '''
        return self.clone(AssignUMI)

    def run(self):
        '''
        '''
        bamFileIn  = readSet + ".umi_merge.bam"
        bamFileOut = readSet + ".primer_clip.bam"
        core.primer_clip.run(cfg,bamFileIn,bamFileOut,False)
        with open(self.verificationFile,'w') as OUT:
            OUT.write("done\n")

    def output(self):
        '''
        '''
        return luigi.LocalTarget(self.verificationFile)
    
class SortBam(luigi.Task):
    ''' Sort the final BAM file, to prepare for downstream variant calling
    '''
    # Parameters
    samplesCfg = luigi.Parameter(description="Config file containing information on readSet to run")
    readSet = luigi.Parameter(description="The readSet to analyze from the config file above")    

    def __init__(self,*args,**kwargs):
        '''
        '''
        super(self.__class__.__name__,self).__init__(*args,**kwargs)
        currentDir = os.path.dirname(os.path.realpath(__file__))
        verificationDir = os.path.join(currentDir,"verification")
        self.verificationFile = os.path.join(verificationDir,
                                             self.__class__.__name__+'.verification.txt')

    def requires(self):
        '''
        '''
        return self.clone(SoftClipBam)
        
    def run(self):
        '''
        '''
        bamFileIn  = readSet + ".primer_clip.bam"
        bamFileOut = readSet + ".bam"
        core.samtools.sort(cfg,bamFileIn,bamFileOut)
        with open(self.verificationFile,'w') as OUT:
            OUT.write("done\n")

    def output(self):
        '''
        '''
        return luigi.LocalTarget(self.verificationFile)

class RunSmCounter(luigi.Task):
    ''' Run smCounter variant calling
    '''
    # Parameters
    samplesCfg = luigi.Parameter(description="Config file containing information on readSet to run")
    readSet = luigi.Parameter(description="The readSet to analyze from the config file above")    

    def __init__(self,*args,**kwargs):
        '''
        '''
        super(self.__class__.__name__,self).__init__(*args,**kwargs)
        currentDir = os.path.dirname(os.path.realpath(__file__))
        verificationDir = os.path.join(currentDir,"verification")
        self.verificationFile = os.path.join(verificationDir,
                                             self.__class__.__name__+'.verification.txt')

    def requires(self):
        '''
        '''
        return self.clone(SortBAM)
    
    def run(self):
        '''
        '''
        numVariants = varcall.sm_counter_wrapper.run(cfg, paramFile)        
        with open(self.verificationFile,'w') as OUT:
            OUT.write("done\n")

    def output(self):
        '''
        '''
        return luigi.LocalTarget(self.verificationFile)

class AnnotateVCF(luigi.Task):
    ''' Annotate a VCF file
    '''
    # Parameters
    samplesCfg = luigi.Parameter(description="Config file containing information on readSet to run")
    readSet = luigi.Parameter(description="The readSet to analyze from the config file above")    
    
    def __init__(self,*args,**kwargs):
        '''
        '''
        super(self.__class__.__name__,self).__init__(*args,**kwargs)
        # initialize logger
        core.run_log.init(readSet)        
        # Add readSet info to global config
        global cfg
        parser = ConfigParser.ConfigParser()
        parser.read(self.samplesCfg)        
        # Handle only 1 readSet at a time for now
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
        
        currentDir = os.path.dirname(os.path.realpath(__file__))
        verificationDir = os.path.join(currentDir,"verification")
        if not os.path.exists(verificationDir):
            os.makedirs(verificationDir)
        self.verificationFile = os.path.join(verificationDir,
                                             self.__class__.__name__+'.verification.txt')

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
        # close log file
        core.run_log.close()        
        # write the verification file
        with open(self.verificationFile,'w') as OUT:
            OUT.write("done\n")

    def output(self):
        '''
        '''
        return luigi.LocalTarget(self.verification_file)
