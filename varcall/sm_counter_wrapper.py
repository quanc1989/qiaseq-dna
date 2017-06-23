import ConfigParser

# our modules
import sm_counter

def run(cfg, paramFile):
   # get read set name
   readSet = cfg.readSet

   # get mean base-level MT depth from disk (needed for smCounter)
   mtDepthMean = None
   for line in open(readSet + ".sumUniformityBase.summary.txt","r"):
      (metricVal,metricName) = line.strip().split("\t")
      if metricName == "mean MT depth":
         mtDepthMean = int(round(float(metricVal)))

   # get reads per MT from disk (needed for smCounter)
   numReadsPerMt = None
   for line in open(readSet + ".sumPrimerMts.summary.txt","r"):
      (metricVal,metricName) = line.strip().split("\t")
      if metricName == "read fragments per MT, mean":
         numReadsPerMt = int(round(float(metricVal)))
         
   # debug check
   if mtDepthMean == None or numReadsPerMt == None:
      raise Exception("could not find mean MT depth or RPMT on disk!")
      
   # get standard smCounter parameters from the main run-params.txt file
   parser = ConfigParser.SafeConfigParser()
   parser.optionxform = str
   parser.read(paramFile)
   cfgSmCounter = {}
   for (paramName, paramVal) in parser.items("smCounter"): 
      cfgSmCounter[paramName] = paramVal
      
   # set up config dictionary to pass to smCounter
   cfgSmCounter["outPrefix"] = readSet
   cfgSmCounter["bamFile"  ] = readSet + ".bam"
   cfgSmCounter["bedTarget"] = cfg.roiBedFile
   cfgSmCounter["mtDepth"  ] = mtDepthMean
   cfgSmCounter["rpb"      ] = numReadsPerMt
   cfgSmCounter["nCPU"     ] = cfg.numCores
   cfgSmCounter["refGenome"] = cfg.genomeFile
   
   # run smCounter variant caller
   smCounterThreshold = sm_counter.main(cfgSmCounter)
   
   # write smCounter threshold to disk file, for main summary table
   fileout = open(readSet + ".smCounter.summary.txt", "w")
   fileout.write("{}\tsmCounter variant calling threshold\n".format(smCounterThreshold))
   fileout.close()

   # make a file that contains variants with PI below-threshold but above 12
   fileout = open(readSet + ".smCounter.GT12PI.txt","w")
   firstLine = True
   idxPI = None
   for line in open(readSet + ".smCounter.all.txt", "r"):
      vals = line.strip().split("\t")
      if firstLine:
         firstLine = False
         idxPI = vals.index("PI")
         fileout.write("read set\t")
         fileout.write(line)
         continue
      predictionIndex = vals[idxPI]
      predictionIndex = int(float(predictionIndex)) if len(predictionIndex) > 0 else 0
      if 12 <= predictionIndex < smCounterThreshold:
         fileout.write(readSet)
         fileout.write("\t")
         fileout.write(line)
   fileout.close()
   
   # return number of primitive variants called
   numVariants = -1
   for line in open(readSet + ".smCounter.cut.txt","r"):
      numVariants += 1
   return numVariants
