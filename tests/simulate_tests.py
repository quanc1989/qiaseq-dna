import test_sm_counter
import utils
import os
import sys


#################################################################################################
##                                                                                             ##
##         Simulates a test suite for the qiaseq-dna repository code                           ##
##                                                                                             ##
##                                                                                             ##
##         This test suite should be run inside a docker container with the example readset    ##
##         Key metrics from the pipeline pertaining to variant calling, umi statistics         ##
##         will be used to validate whether the test was successful                            ##
#################################################################################################
## To do :
##   1. Have a cleaner  way to initialize the test directories

## Run this test as follows :
## sudo docker pull rpadmanabhan9/qiaseq-dna
## sudo docker run -v /home/your_fav_dir:/mnt/qiaseq-run/ rpadmanabhan9/qiaseq-dna /bin/bash -c "cd /srv/qgen/code/; git clone https://qiauser:anz2teu@github.com/qiaseq/qiaseq-dna.git; cd qiaseq-dna; python tests/simulate_tests.py True;"


########## smCounter ##########
print "Starting to run tests.......\n------------------"
## Setup directory structure
if sys.argv[1].lower() == "true":
    os.system("mkdir -p /mnt/qiaseq-run/test_sm_counter/")
    os.chdir("/mnt/qiaseq-run/test_sm_counter/")
    os.system("cp /srv/qgen/code/qiaseq-dna/run_sm_counter_v1.params.txt ./")
    print "Running sm_counter pipeline\n"
    test_sm_counter.run_pipeline()
# Run tests
test_sm_counter.validate_umi_depth_metrics("NEB_S2.umi_depths.summary.txt")
test_sm_counter.validate_umi_read_frags_metrics("NEB_S2.umi_frags.summary.txt")
test_sm_counter.validate_umi_filter_metrics("NEB_S2.umi_filter.summary.txt")
test_sm_counter.validate_variant_calling_metrics("NEB_S2.vcf_complex.summary.txt")
print "Passed all sm_counter tests\n------------------"
