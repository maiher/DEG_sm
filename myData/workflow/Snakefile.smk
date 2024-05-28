# optional: provide config file here:
#configfile: "config.yml"
from snakemake.utils import min_version
min_version("5.3.0")



rule fastp:
    """
    Run fastp on a FASTQ file.
    """
    output:
        'results/fastp/{sample}_R1.qc.fq.gz',
        'results/fastp/{sample}_R2.qc.fq.gz',
        'results/fastp/{sample}.fastp.html',
        'results/fastp/{sample}.fastp.json'
        #'{sample_id}_fastq.fastp.gz'
        # this would not run:
        #'results/0246_2_TP04_026_FC_MOD_1TR1_SL1_S26_L002_R1_001_50k_fastq.fastp.gz',
        #'results/0246_2_TP04_026_FC_MOD_1TR1_SL1_S26_L002_R2_001_50k_fastq.fastp.gz',
        #'results/fastp.html',
        #'results/fastp.json'
    input:
        'myData/{sample}_R1_001_50k_fastq.gz',
        'myData/{sample}_R2_001_50k_fastq.gz'
    log:
        # 2>&1 fwd both
        'logs/fastp/{sample}.log'
    shell:
        """
        time fastp -i {input[0]} -o {output[0]} -I {input[1]} -O {output[1]} 2>&1 {log}
        mv fastp.html {output[2]}
        mv fastp.json {output[3]}
        """

'''
# Run:
cd /home/AD/hermi/sm_tutorial/BYOC
snakemake -c 1 \
-s myData/workflow/Snakefile_fastp_wildcards_simple.smk \
results/fastp/0246_2_TP04_026_FC_MOD_1TR1_SL1_S26_L002_R1.qc.fq.gz \
results/fastp/0246_2_TP04_026_FC_MOD_1TR1_SL1_S26_L002_R1.qc.fq.gz
'''
