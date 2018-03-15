import utils 

## To do :
## 1. Move base metric values to a config file


## Ground truth values 
## UMI related metrics 
UMIs = 434994
READS_PER_MT = 7.8
READS_FROM_INTERNAL_DOWNSTREAM_PRIMING = 21.9
## UMI depth related metrics 
MEAN_UMI_DEPTH = 1905.71
MEAN_READ_DEPTH = 24393.80 
## UMI filters 
READS_DROPPED_R1_R2_UNMAPPED = 119428  
READS_DROPPED_MAPPING_FILTERS = 143845
READS_DROPPED_OFF_TARGET = 246741 
## smCounter
VARIANTS_CALLED = 18 


def run_pipeline():
    """ Run a full smCounter pipeline on the test readset 
    """

    cmd = (
        """ python /srv/qgen/code/qiaseq-dna/run_qiaseq_dna.py run_sm_counter_v1.params.txt v1 """
        """ single NEB_S2 > run.log 2>&1 """
    )
    return_code =  utils.run_shell_cmd(cmd)
    assert return_code == 0, """ smCounter pipeline failed !\n
    The command : \n %s \n exited with non zero return code.
    Please check the logs ! """

def validate_umi_depth_metrics(metrics_file):
    """ Validate depth related umi metrics 

    metrics_file : str ; path to the metrics file 
    """
    
    ## Parse the values from the summary file first
    metrics = utils.parse_summary_file(metrics_file)
    ## Check values
    assert utils.compare_with_tolerance(metrics["mean UMI depth"],MEAN_UMI_DEPTH,0.02) is True
    assert utils.compare_with_tolerance(metrics["mean read depth"],MEAN_READ_DEPTH,0.02) is True
    
def validate_umi_read_frags_metrics(metrics_file):
    """ Validate metrics related to read frags and umi

    metrics_file : str ; path to the metrics file 
    """

    ## Parse the values from the summary file first
    metrics = utils.parse_summary_file(metrics_file)
    assert utils.compare_with_tolerance(metrics["UMIs"],UMIs,0.02) is True 
    assert utils.compare_with_tolerance(metrics["read fragments per MT, mean"],READS_PER_MT,0.02) is True
    assert utils.compare_with_tolerance(metrics["% of reads from internal downstream priming"],READS_FROM_INTERNAL_DOWNSTREAM_PRIMING,0.02) is True
    
def validate_umi_filter_metrics(metrics_file):
    """ Validate metrics related to how the reads were filtered

    metrics_file : str ; path to the metrics file 
    """

    ## Parse the values from the summary file first
    metrics = utils.parse_summary_file(metrics_file)
    assert utils.compare_with_tolerance(metrics["read fragments dropped, R1 or R2 not mapped"],READS_DROPPED_R1_R2_UNMAPPED,0.02) is True
    assert utils.compare_with_tolerance(metrics["read fragments dropped, not passing mapping filters (low mapq, split alignments, discordant pairs, etc.)"],READS_DROPPED_MAPPING_FILTERS,0.02) is True
    assert utils.compare_with_tolerance(metrics["read fragments dropped, off target"],READS_DROPPED_OFF_TARGET,0.02) is True
    
def validate_variant_calling_metrics(metrics_file):
    """ Validate metrics related 
    
    metrics_file : str ; path to the metrics file 
    """
    
    ## Parse the values from the summary file first
    metrics = utils.parse_summary_file(metrics_file)
    assert utils.compare_with_tolerance(metrics["variants called by smCounter"],VARIANTS_CALLED,0) is True
