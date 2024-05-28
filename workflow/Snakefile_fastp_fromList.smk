# optional: provide config file here:
#configfile: "config.yml"
from snakemake.utils import min_version
min_version("5.3.0")


# get sample IDs based on fq paths:
sample_list = []
with open(config["sample_list"], "r") as fhin:
    for line in fhin:
    # strip the newline character from the end of the line
    # get R1 files basename
    # strip the ending, including R1/R2 = sample id is left
        query="_R1_001_50k_fastq.gz"    # config.yml!
        fq_path = line.strip()
        if query in fq_path:
            #print(fq_path)
            sample = os.path.basename(fq_path)
            #print(sample)
            sample = sample.replace(query, '')
            sample_list.append(sample)

# get blacklisted sample IDs based on fq paths:
black_list = []
with open(config["black_list"], "r") as fhin:
    for line in fhin:
    # strip the newline character from the end of the line
    # get R1 files basename
    # strip the ending, including R1/R2 = sample id is left
        fq_path = line.strip()
        for sample in sample_list:
            if sample in fq_path:
                print(sample)
                sample_list.remove(sample)
print(sample_list)




# if no target is given at the command line, Snakemake will define the first rule 
#   of the Snakefile as the target
rule all:
    input:
    # put output of targeted rule here!!!
    # here: output of 'rule fastp'
        expand('results/fastp/{sample}_R1.qc.fq.gz', sample = sample_list),
        expand('results/fastp/{sample}_R2.qc.fq.gz', sample = sample_list),
        expand('results/fastp/{sample}.fastp.html', sample = sample_list),
        expand('results/fastp/{sample}.fastp.json', sample = sample_list)
# to do: get {sample} wo path, expand in input!
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
        'results/logs/fastp/{sample}.log'
    params:
        config["fastp_params"]
    shell:
        """
        time fastp -i {input[0]} -o {output[0]} -I {input[1]} -O {output[1]} {params} &> {log}
        mv fastp.html {output[2]}
        mv fastp.json {output[3]}
        """

''''
# Run:
conda activate DEGpipe01
cd /home/AD/hermi/sm_tutorial/BYOC
snakemake -c 1 \
--configfile config.yml \
-s workflow/Snakefile_fastp_fromList.smk
'''
