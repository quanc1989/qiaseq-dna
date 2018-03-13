import os

def run(cfg):
    # get read set name
    readSet = cfg.readSet

    if cfg.outputDetail.lower() == "false":
        umiFilterFile = "umi_filter"
    else:
        umiFilterFile = "umi_filter.detail"

    # concatenate summary files
    with open(readSet + ".sumAll.summary.txt", "w") as OUT:
        for fileType in ("prep", umiFilterFile, "sum.primer.umis", "umi_frags.len-distrib.txt", "sum.uniformity.primer", "umi_depths", "smCounter", "vcf_complex"):
            '''
            prep --> trimming metrics ; core/prep.py
            umiFilterFile --> primer finding metrics : core/umi_filter.py
            sum.primer.umis --> umi information for each primer : metrics/sum_primer_umis.py
            umi_frags.len-distrib.txt --> fragment length and reads per umi info: metrics/umi_frags.py
            umi_depths.summary.txt -->  umi depth and lod info : metrics/umi_depths.py
            smCounter --> metrics from smCounter : sm_counter_wrapper.py
            vcf_complex --> metrics from reconstructing complex variants from primitives : annotate/vcf_complex.py
            '''
            fileName = readSet + "." + fileType + ".summary.txt"
            if os.path.isfile(fileName):
                for line in open(fileName):
                    OUT.write(line)
        
    
