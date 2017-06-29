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


## Run this test as follows :
## sudo docker pull rpadmanabhan9/qiaseq-dna
## sudo docker run -v /home/your_fav_dir:/mnt/qiaseq-run/ rpadmanabhan9/qiaseq-dna /bin/bash -c "cd /srv/qgen/code/; git clone https://qiauser:anz2teu@github.com/qiaseq/qiaseq-dna.git; cd /mnt/qiaseq-run/; python /srv/qgen/code/qiaseq-dna/simulate_tests.py;" 


########## smCounter ##########
## Setup directory structure
os.system("mkdir -p /mnt/qiaseq-run/test_sm_counter/")
os.system("cd /mnt/qiaseq-run/test_sm_counter/")
os.system("cp /srv/qgen/code/qiaseq-dna/run_sm_counter.params.txt ./")
## Run tests
test_sm_counter.run_pipeline()
test_sm_counter.validate_umi_depth_metrics("NEB_S2.umi_depths.summary.txt")
test_sm_counter.validate_umi_read_frags_metrics("NEB_S2.umi_frags.summary.txt")
test_sm_counter.validate_umi_filter_metrics("NEB_S2.umi_filter.summary.txt")
test_sm_counter.validate_variant_calling_metrics("NEB_S2.vcf_complex.summary.txt")
########## Consensus Calling ##########
## switch to consensus dir
## run tests 

########## Duplicate Removal ##########
