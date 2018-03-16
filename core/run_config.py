import ConfigParser
from multiprocessing.dummy import cpu_count as cpu_count

#--------------------------------------------------------------------------------------
def run(readSet,paramFile):

   # read parameter file
   parser = ConfigParser.SafeConfigParser()
   parser.optionxform = str
   parser.read(paramFile)

   # copy all options to a config object - both general params, and params for this readSet
   cfg = lambda:0
   cfg.__dict__["readSet"] = readSet
   for section in ("general", readSet):
      for (paramName, paramVal) in parser.items(section):
         if paramName in cfg.__dict__:
            raise Exception("Config file contains duplicate specification of parameter: " + paramName)
         cfg.__dict__[paramName] = paramVal
         print(paramName, paramVal)

   # use all cores if numCores = 0
   if int(cfg.numCores) == 0:
      if cfg.sampleType.lower() in ['normal','tumor']: # use half the cores for tumor normal readSets
         cfg.numCores = str(cpu_count()/2)
      cfg.numCores = str(cpu_count())

   # convert some params to boolean
   cfg.deleteLocalFiles = cfg.deleteLocalFiles.lower() == "true"

   # return config object
   return cfg
