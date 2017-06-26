import os
import subprocess

def run(cfg, vcfFileIn, vcfFileOut):
   #-------------------------------
   # part 1 - annotate VCF
   #-------------------------------

   # get params
   javaExe     = cfg.javaExe
   snpEffPath  = cfg.snpEffPath
   snpEffConfig= cfg.snpEffConfig
   dbSnpFile   = cfg.dbSnpFile
   cosmicFile  = cfg.cosmicFile
   clinVarFile = cfg.clinVarFile
   readSet     = cfg.readSet
   logFile = readSet + ".vcf_annotate.snpEffSift.log"
      
   # snpEff eff (NOTE: removed the "-t" option for multi-thread operation, due to multiple customer cases of "java.util.ConcurrentModificationException")
   print("vcf_annotate: starting snpEff eff...")
   cmd = javaExe + " -Xmx6G" \
   + " -jar " + snpEffPath + "snpEff.jar " \
   + " -c " + snpEffConfig \
   + " -noLog -v -noStats -noMotif -noNextProt GRCh37.75 " \
   + vcfFileIn \
   + " > " + readSet + ".temp0.vcf" \
   + " 2> " + logFile
   print(cmd)
   subprocess.check_call(cmd, shell=True)

   # add dbSNP 
   print("vcf_annotate: starting snpSift annoate dbSNP...")
   cmd = javaExe + " -Xmx6G" \
   + " -jar " + snpEffPath + "SnpSift.jar annotate" \
   + " -c " + snpEffConfig \
   + " -id " + dbSnpFile \
   + " "   + readSet + ".temp0.vcf" \
   + " > " + readSet + ".temp1.vcf" \
   + " 2>> " + logFile
   subprocess.check_call(cmd, shell=True)

   # add cosmic
   print("vcf_annotate: starting snpSift annoate Cosmic...")
   cmd = javaExe + " -Xmx6G" \
   + " -jar " + snpEffPath + "SnpSift.jar annotate" \
   + " -c " + snpEffConfig \
   + " -id " + cosmicFile \
   + " "   + readSet + ".temp1.vcf" \
   + " > " + readSet + ".temp0.vcf" \
   + " 2>> " + logFile
   subprocess.check_call(cmd, shell=True)
   
   # add clinvar
   print("vcf_annotate: starting snpSift annoate ClinVar...")
   cmd = javaExe + " -Xmx6G" \
   + " -jar " + snpEffPath + "SnpSift.jar annotate" \
   + " -c " + snpEffConfig \
   + " -id " + clinVarFile \
   + " "   + readSet + ".temp0.vcf" \
   + " > " + readSet + ".temp1.vcf" \
   + " 2>> " + logFile
   subprocess.check_call(cmd, shell=True)

   # rename final VCF file
   os.rename(readSet + ".temp1.vcf", vcfFileOut)
   os.remove(readSet + ".temp0.vcf")

   #-------------------------------
   # part 2 - convert VCF to table
   #-------------------------------
   
   # open output file, prepare column headers
   fileout = open(readSet + ".smCounter.anno.txt", "w")   
   colNames = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]
   tagsNeeded = ("TYPE", "DP", "MT", "UMT", "PI", "THR", "VMT", "VMF", "VSM")
   colNames.extend(tagsNeeded)
   numAnnCols = None
      
   # dump the final VCF file in tab-delimited text format, for easy mulit-sample concatenation 
   for line in open(vcfFileOut, "r"):

      # parse snpEff ANN column headers
      if line.startswith("##INFO=<ID=ANN,"):
         idx = line.find("Functional annotations: '")
         line = line[idx:]
         idx = line.find("'")
         line = line[idx+1:]
         idx = line.find("'")
         line = line[:idx]
         vals = line.strip().split(" | ")
         numAnnCols = len(vals) + 1  # including "ANN other"
         colNames.extend(vals)
         colNames.extend(("ANN other","INFO other"))
         fileout.write("\t".join(colNames))
         fileout.write("\n")
         continue      

      # skip other headers
      if line.startswith("#"):
         continue

      # parse line (drop FORMAT and GENOTYPE - useless fakes)
      (chrom, pos, id, ref, alt, qual, filter, info, fmt, geno) = line.strip().split("\t")
      
      # remove ".0" in QUAL, added by snpEff
      idx = qual.find(".")
      if idx >= 0:
         qual = qual[:idx]
      
      # save main fields
      outvec = []
      outvec.extend((chrom, pos, id, ref, alt, qual, filter))

      # get INFO field, split to output table columns
      vals = info.split(";")
      for idx in range(len(tagsNeeded)):
         val = vals[idx]
         if val.find("=") == -1:
            print(idx,val,line)
            raise Exception("bad INFO tag parsing")
         (tagName, tagVal) = val.split("=")
         if tagName != tagsNeeded[idx]:
            print(idx,val,line)
            raise Exception("VCF header INFO column tags not in expected sort order!")
         outvec.append(tagVal)

      # split ANN field to output table columns
      annTagFound = False
      if len(vals) > len(tagsNeeded):
         val = vals[len(tagsNeeded)]
         if val.startswith("ANN="):
            annTagFound = True
            (tagName, tagVal) = val.split("=")
            anns = tagVal.split(",")
            outvec.extend(anns[0].split("|"))  # most-damaging transcript in own columns
            if len(anns) > 1:                  # all other transcripts in one column
               outvec.append(",".join(anns[1:]))
            else:
               outvec.append("")

      # handle missing ANN field (snpEff cannot annotate variants with "N" on chrMT!)
      if not annTagFound:
         for i in range(numAnnCols):
            outvec.append("")
      
      # INFO tags beyond the ANN tag, all in one cell
      otherInfo = vals[len(tagsNeeded)+1:]
      if len(otherInfo) > 0:
         outvec.append(";".join(otherInfo))
      else:
         outvec.append("_NONE_")  # avoid empty last column, more difficult to parse
            
      # output row
      if len(outvec) != len(colNames):
         print(len(outvec), len(colNames), line)
         raise Exception("Unexpected INFO tag count")
      fileout.write("\t".join(outvec))
      fileout.write("\n")

#----------------------------------------------------------------------------------------------
# pythonism to run from the command line
#----------------------------------------------------------------------------------------------
if __name__ == "__main__":
   import sys
   cfg = lambda:0
   cfg.javaExe     = sys.argv[1]
   cfg.snpEffPath  = sys.argv[2] 
   cfg.snpEffConfig= sys.argv[3] 
   cfg.dbSnpFile   = sys.argv[4] 
   cfg.cosmicFile  = sys.argv[5] 
   cfg.clinVarFile = sys.argv[6] 
   cfg.readSet     = sys.argv[7]
   vcfFileIn  = cfg.readSet + ".smCounter.cplx.vcf"
   vcfFileOut = cfg.readSet + ".smCounter.anno.vcf"
   run(cfg, vcfFileIn, vcfFileOut)
