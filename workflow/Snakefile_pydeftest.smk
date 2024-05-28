# optional: provide config file here:
#configfile: "config.yml"
import os
from snakemake.utils import min_version
min_version("5.3.0")





#def get_raw_fqs(wildcards):
samples = []
with open(config["sample_list"], "r") as fhin:
    for line in fhin:
        # strip the newline character from the end of the line
        # get R1 files basename
        # strip the ending, including R1/R2 = sample id is left
        query="_R1_001_50k_fastq.gz"
        fq_path = line.strip()
        if query in fq_path:
            sample = os.path.basename(fq_path)
            sample = sample.strip(query)
            samples.append(sample)
print(samples)