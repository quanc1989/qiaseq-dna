from tests import test_sm_counter
from tests import test_consensus_calling
from tests import duplicate_removal 
from tests import utils



#################################################################################################
##                                                                                             ##
##         Simulates a test suite for the qiaseq-dna repository code                           ##
##                                                                                             ##
##         Will run the 3 pipelines in this repository , i.e.                                  ##
##           1.) Conensensus Calling                                                           ##
##           2.) smCounter (UMI aware variant calling)                                         ##
##           3.) Duplicate removal                                                             ##
##                                                                                             ##
##         This test suite should be run inside a docker container with the example readset    ##
##         Key metrics from the pipeline pertaining to variant calling, umi statistics         ##
##         will be used to validate whether the test was successful                            ## 
#################################################################################################
## To do :
##   1. Have a cleaner  way to initialize the test directories
##   2. Complete consensus and duplicate tests 

########## smCounter ##########
## Setup directory structure
os.system("mkdir -p /mnt/qiaseq-run/test_sm_counter/")
os.system("cd /mnt/qiaseq-run/test_sm_counter/")
os.system("cp /srv/qgen/code/qiaseq-dna/run_sm_counter.params.txt ./")
## Run tests
test_sm_counter.run_pipeline()
test_sm_counter.validate_umi_metrics()
########## Consensus Calling ##########
## switch to consensus dir 


########## Duplicate Removal ##########
