# optional: provide config file here:
#configfile: "config.yml"
from snakemake.utils import min_version
min_version("5.3.0")

# set sample ids for whole workflow:
#SAMPLES = ["SRR935090", "SRR935091", "SRR935092"]


# define an empty 'samples' dictionary
fq_paths = []
sample_ids = []
mates = []

# read the sample list file and populate the dictionary
#with open("samples.tsv", "r") as fhin:
# populate from config-yml:
with open(config["sample_list"], "r") as fhin:
    for line in fhin:
        # strip the newline character from the end of the line
        # then split by tab character to get the sample id and url
        query="_R1_001_50k_fastq.gz"
        fq_path = line.strip()
        sample_id = 
        if query in fq_path:
            fq_paths.append(fq_path) 
            mates.append("R1")
        else:
            fq_paths.append(fq_path) 
            mates.append("R2")

# dict1: fq_path : sample_id
sample_id_dict = dict(zip(fq_paths, sample_ids))
# dict2: fq_path : mate
mate_dict = dict(zip(fq_paths, mates))

print(fq_paths)
print(sample_ids)
print(mates)
print(sample_id_dict)
print(mate_dict)


'''
rule fastqc:
    """
    Run FastQC on a FASTQ file.
    """
    output:
        '{sample_id_01}.qc.fq.gz',
        '{sample_id_02}.qc.fq.gz',
        '{sample_id_01}.fastp.html',
        '{sample_id_01}.fastp.json'
        #'{sample_id}_fastq.fastp.gz'
        # this would not run:
        #'results/0246_2_TP04_026_FC_MOD_1TR1_SL1_S26_L002_R1_001_50k_fastq.fastp.gz',
        #'results/0246_2_TP04_026_FC_MOD_1TR1_SL1_S26_L002_R2_001_50k_fastq.fastp.gz',
        #'results/fastp.html',
        #'results/fastp.json'
    input:
        '{sample_id_01}_fastq.gz',
        '{sample_id_02}_fastq.gz'
        #'{sample_id}.gz'
        #'0246_2_TP04_026_FC_MOD_1TR1_SL1_S26_L002_R1_001_50k_fastq.gz',
        #'0246_2_TP04_026_FC_MOD_1TR1_SL1_S26_L002_R2_001_50k_fastq.gz'
    #log:
        # 2>$1 fwd both
        #"logs/fastp/{sample_id}.log"
    shell:
        """
        time fastp -i {input[0]} -o {output[0]} -I {input[1]} -O {output[1]}
        mv fastp.html {output[2]}
        mv fastp.json {output[3]}
        """

#    input:
#        #"home/AD/hermi/sm_tutorial/BYOC/myData/0246_2_TP04_026_FC_MOD_1TR1_SL1_S26_L002_R1_001_50k_fastq.gz",
#        #"home/AD/hermi/sm_tutorial/BYOC/myData/0246_2_TP04_026_FC_MOD_1TR1_SL1_S26_L002_R2_001_50k_fastq.gz"
# .snakemake/log/2024-05-27T125611.315731.snakemake.log


'''