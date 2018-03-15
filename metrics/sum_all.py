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
        for fileType in ("prep", "align", umiFilterFile, "umi_frags", "sum.uniformity.primer", "umi_depths", "smCounter", "vcf_complex"):
            '''
            prep, align --> trimming metrics ; core/prep.py, misc/tvc.py
            umiFilterFile --> primer finding metrics : core/umi_filter.py
            umi_frags --> reads per umi info : metrics/umi_frags.py
            umi_depths.summary.txt -->  umi depth and lod info : metrics/umi_depths.py
            smCounter --> metrics from smCounter : sm_counter_wrapper.py
            vcf_complex --> metrics from reconstructing complex variants from primitives : annotate/vcf_complex.py
            '''
            fileName = readSet + "." + fileType + ".summary.txt"
            if os.path.isfile(fileName):
                for line in open(fileName):
                    OUT.write(line)
        
    
