import os
import os.path
import subprocess

# our modules
import samtools

#---------------------------------------------------
def run(cfg,bamFileIn,bamFileOut):
   print("dedup: starting...")

   # get params
   readSet = cfg.readSet
   javaExe = cfg.javaExe
   picardJar = cfg.picardJar
   deleteLocalFiles = cfg.deleteLocalFiles
   tagNameUmiSeq = cfg.tagNameUmiSeq

   # sort the input bam
   samtools.sort(cfg,bamFileIn,readSet + ".dedup.temp0.bam")
   print("dedup: done with sorting BAM file")
   
   # delete input file if requested
   if deleteLocalFiles and len(os.path.dirname(bamFileIn)) == 0:
      os.remove(bamFileIn)
   
   # run Picard MarkDuplicates
   cmd = javaExe + " -jar " \
   + picardJar + " MarkDuplicates" \
   + " I=" + readSet + ".dedup.temp0.bam" \
   + " O=" + bamFileOut \
   + " M=" + readSet + ".dedup.markdup_metrics.txt" \
   + " BARCODE_TAG=" + tagNameUmiSeq \
   + " REMOVE_DUPLICATES=True" \
   + " >> " + readSet + ".dedup.shell.log 2>&1"
   subprocess.check_call(cmd,shell=True)
   print("dedup: done with Picard MarkDuplicates")
   os.remove(readSet + ".dedup.temp0.bam")
