import subprocess
import os
import os.path

def sort(cfg,bamFileIn,bamFileOut):
   # params
   deleteLocalFiles = cfg.deleteLocalFiles
   readSet          = cfg.readSet
   samtoolsDir      = cfg.samtoolsDir
   samtoolsMem      = cfg.samtoolsMem
   numCores         = cfg.numCores

   # sort
   cmd = samtoolsDir + "samtools sort"  \
   + " -m " + samtoolsMem \
   + " -@ " + numCores \
   + " -T " + readSet \
   + " -o " + bamFileOut \
   + "    " + bamFileIn \
   + " >> " + readSet + ".samtools.shell.log 2>&1"
   subprocess.check_call(cmd, shell=True)
   
   # index
   cmd = samtoolsDir + "samtools index " + bamFileOut
   subprocess.check_call(cmd, shell=True)

   # delete input file if requested
   if deleteLocalFiles and len(os.path.dirname(bamFileIn)) == 0:
      os.remove(bamFileIn)
   