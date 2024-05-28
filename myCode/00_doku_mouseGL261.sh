# md5sums: DONE, good!
# -------------------------------------------------------------
/data2/240409_RKoenig_mouseGL261_round2/240409_getData.sh
#   output
/data2/240409_RKoenig_mouseGL261_round2/240409_md5_fastqs.txt
#   TRNO md5sums


# fastp QC, DONE
# -------------------------------------------------------------
/data2/240409_RKoenig_mouseGL261_round2/240409_qc.sh
#   input
/data2/240409_RKoenig_mouseGL261_round2/240409_rawFastq
#   output: only with -G option: polyG-trimming disabled
#       like it was in all scripts before: more reads. tested with human SAMHD1 data
/data2/240409_RKoenig_mouseGL261_round2/240409_qcFastq


# STAR mapping: done
# -------------------------------------------------------------
/data2/240409_RKoenig_mouseGL261_round2/240409_align.sh
# used STAR genome index by Csaba:
/data/public/RefSeq.GRCm39/STARindexes
#   input
/data2/240409_RKoenig_mouseGL261_round2/240409_qcFastq
#   output
/data2/240409_RKoenig_mouseGL261_round2/240409_align


# stats table - QC/mapping: done
# -------------------------------------------------------------
/data2/240409_RKoenig_mouseGL261_round2/240409_qc_map_stats.sh
#   input
#   all /data2/240122_RKoenig_humanHSPC/240125_align/ SAMPLE Log.final.out
#   all /data2/240122_RKoenig_humanHSPC/240125_qcFastq/ SAMPLE .log
#   output
/data2/240409_RKoenig_mouseGL261_round2/240409_QCstats.csv
# looks good



# featureCounts: get gene len and count data: DONE
# -------------------------------------------------------------
#Features and meta-features
#Each entry in the provided annotation file is taken as a feature (e.g. an exon). A meta-feature 
# is the aggregation of a set of features (e.g. a gene). The featureCounts program uses the gene_id 
# attribute available in the GTF format annotation (or the GeneID column in the SAF format annotation) 
# to group features into meta-features, ie. features belonging to the same meta-feature have the same 
# gene identifier.
# featureCounts can count reads at either feature level or at meta-feature level. When summarizing reads 
# at meta-feature level, read counts obtained for features included in the same meta-feature will be 
# added up to yield the read count for the corresponding meta-feature.

/data2/240409_RKoenig_mouseGL261_round2/240409_featureCounts.sh
#   output dir:
/data2/240409_RKoenig_mouseGL261_round2/240411_featureCounts
# real	173m4.672s



# normalize read counts: OK!
#  and edit colnames for output table!
# -----------------# --------------------
# conda base
cp /data2/240122_RKoenig_humanHSPC/240126_featureCount/240126_readCount_norm.R \
/data2/240409_RKoenig_mouseGL261_round2/240411_featureCount/240411_readCount_norm.R 
#   input
/data2/240409_RKoenig_mouseGL261_round2/240411_featureCount/mouseGL261_readCounts.csv
#   output:
/data2/240409_RKoenig_mouseGL261_round2/240411_featureCount/mouseGL261_readCounts_norm.csv





# DESeq2 + PCAs:
# ---------------
conda activate R-4.2_cpFix
outDir=/data2/240409_RKoenig_mouseGL261_round2/240411_DESeq2
mkdir $outDir

# R:
cp /data2/240122_RKoenig_humanHSPC/240126_DESeq2/240212_DESeq2_preprocessing.rmd $outDir
cp /data2/240122_RKoenig_humanHSPC/240126_DESeq2/240213_DESeq2_PCA_all.rmd $outDir

rstudio
# input
/data2/240409_RKoenig_mouseGL261_round2/240411_featureCount/mouseGL261_readCounts_norm.csv
/data2/240409_RKoenig_mouseGL261_round2/240409_RNAseq_mice_samples.csv

# prep dds:
/data2/240409_RKoenig_mouseGL261_round2/240411_DESeq2/240411_DESeq2_preprocessing.rmd
/data2/240409_RKoenig_mouseGL261_round2/240411_DESeq2/240411_DESeq2_preprocessing.html
/data2/240409_RKoenig_mouseGL261_round2/240411_DESeq2/240411_rawDds.rds

# PCA & DEG
# test1: all samples
# test1.1: AF_2_3 removed (outlier in PCA, off target control)
/data2/240409_RKoenig_mouseGL261_round2/240411_DESeq2/240411_test1_DESeq2_WTvsKOvsCtrl.rmd
/data2/240409_RKoenig_mouseGL261_round2/240411_DESeq2/240411_test1.1_DESeq2_WTvsKOvsCtrl_noAF2.3.rmd

/data2/240409_RKoenig_mouseGL261_round2/240411_DESeq2/240411_Test1_PCA_CtrlvsKOvsWT.pdf
/data2/240409_RKoenig_mouseGL261_round2/240411_DESeq2/240411_Test1.1_PCA_CtrlvsKOvsWT.pdf

/data2/240409_RKoenig_mouseGL261_round2/240411_DESeq2/240411_Test1_DESeq2_Off_target_control_vs_WT.csv
/data2/240409_RKoenig_mouseGL261_round2/240411_DESeq2/240411_Test1_DESeq2_SAMHD1_KO_vs_Off_target_control.csv
/data2/240409_RKoenig_mouseGL261_round2/240411_DESeq2/240411_Test1_DESeq2_SAMHD1_KO_vs_WT.csv

/data2/240409_RKoenig_mouseGL261_round2/240411_DESeq2/240411_Test1.1_DESeq2_Off_target_control_vs_WT.csv
/data2/240409_RKoenig_mouseGL261_round2/240411_DESeq2/240411_Test1.1_DESeq2_SAMHD1_KO_vs_Off_target_control.csv
/data2/240409_RKoenig_mouseGL261_round2/240411_DESeq2/240411_Test1.1_DESeq2_SAMHD1_KO_vs_WT.csv
# also xlsx

# volcanos
/data2/240409_RKoenig_mouseGL261_round2/240411_DESeq2/240411_test1_volcanos.rmd
/data2/240409_RKoenig_mouseGL261_round2/240411_DESeq2/240411_test1.1_volcanos_noAF2.3.rmd

/data2/240409_RKoenig_mouseGL261_round2/240411_DESeq2/240411_Test1_volcano_l2FCTH1_Off_target_control_vs_WT.pdf
/data2/240409_RKoenig_mouseGL261_round2/240411_DESeq2/240411_Test1_volcano_l2FCTH1_SAMHD1_KO_vs_Off_target_control.pdf
/data2/240409_RKoenig_mouseGL261_round2/240411_DESeq2/240411_Test1_volcano_l2FCTH1_SAMHD1_KO_vs_WT.pdf

/data2/240409_RKoenig_mouseGL261_round2/240411_DESeq2/240411_Test1.1_volcano_l2FCTH1_Off_target_control_vs_WT.pdf
/data2/240409_RKoenig_mouseGL261_round2/240411_DESeq2/240411_Test1.1_volcano_l2FCTH1_SAMHD1_KO_vs_Off_target_control.pdf
/data2/240409_RKoenig_mouseGL261_round2/240411_DESeq2/240411_Test1.1_volcano_l2FCTH1_SAMHD1_KO_vs_WT.pdf




# GO & KEGG NOT DONE
# ---------------
outDir=/data2/240409_RKoenig_mouseGL261_round2/240412_GO_KEGG
mkdir $outDir
# R: make tables

cp /data2/240122_RKoenig_humanHSPC/240319_GO_test3/240227_test3.4_GO_KEGG_table.Rmd $outDir
cp /data2/240122_RKoenig_humanHSPC/240319_GO_test3/240327_test3.4_GO_top25_l2fc_plots.Rmd $outDir
cp /data2/240122_RKoenig_humanHSPC/240319_GO_test3/240403_test3.4_GO_biomarkers_plots.Rmd $outDir
cp /data2/240122_RKoenig_humanHSPC/240319_GO_test3/240403_test3.4_GO_top25_pGR_plots.Rmd $outDir




#   'org.Mm.eg.db' not found
# mouse GO DB is necessary, but can't be installed with R 4.2.0!

# solution: transfered to windows, installation was fine with R 4.2.3
# I copied all into the project folder here to keep track but 


# !!! THIS won't run with R 4.2.0 !!!

# script:
/data2/240122_RKoenig_mouseGL261/240222_GO_KEGG/240222_GO_KEGG_table.rmd
# output:
/data2/240122_RKoenig_mouseGL261/240222_GO_KEGG/240222_GO_KEGG_table.html
/data2/240122_RKoenig_mouseGL261/240222_GO_KEGG/240222_mouseGL261_enrichGO_SAMHD1_KO_vs_WT.csv
/data2/240122_RKoenig_mouseGL261/240222_GO_KEGG/240222_mouseGL261_gseKEGG_SAMHD1_KO_vs_WT.csv
/data2/240122_RKoenig_mouseGL261/240222_GO_KEGG/240222_mouseGL261_enrichGO_SAMHD1_KO_vs_control.csv
/data2/240122_RKoenig_mouseGL261/240222_GO_KEGG/240222_mouseGL261_gseKEGG_SAMHD1_KO_vs_control.csv
/data2/240122_RKoenig_mouseGL261/240222_GO_KEGG/240222_mouseGL261_enrichGO_control_vs_WT.csv
/data2/240122_RKoenig_mouseGL261/240222_GO_KEGG/240222_mouseGL261_gseKEGG_control_vs_WT.csv


# biomarkers: extract from readcounts, DEG & GO results
# ------------------------
/data2/240122_RKoenig_mouseGL261/240226_biomarkers.sh




# qiime: check unmapped reads for contamination
# ------------------------

# 1. extract unmapped reads with samtools:

# Remove reads that mapped to reference
#   then sort resulting BAM by readname:
# -b		output=BAM
# -f		only include reads with all  of the FLAGs in INT present (12 = read & mate unpaired)
# -F		only include reads with none of the FLAGS in INT present (256=not primary alignment)
# -n		Sort by read name


conda activate bedtools

inDir=/data2/240122_RKoenig_mouseGL261/240215_align
cd $inDir

for i in *.bam; do
    echo $i;
    #echo "${i%%.*}"'.unmapped.bam';
    bam="${i%%.*}"'.nameSort.unmapped.bam';
    fq1="${i%%.*}"'.nameSort.unmapped_R1.fq';
    fq2="${i%%.*}"'.nameSort.unmapped_R2.fq'; 
    samtools view -hb -f 12 -F 256 $i |\
    samtools sort -n - -o $bam;
    bamToFastq -i $bam -fq $fq1 -fq2 $fq2
done  
# still empty! check star!!!


/opt/samtools/bin/samtools view -hb -f 12 -F 256 \
$PREFIX$OUT2 2>> $PREFIX$OUT1 | /opt/samtools/bin/samtools sort \
-n - -@ 32 -o $PREFIX$OUT3 \
2>> $PREFIX$OUT1


samtools sort -n -o $outFile $i;

inDir=/data2/240122_RKoenig_mouseGL261/240215_align
cd $inDir
for i in *.bam; do
    echo $i;
    #echo "${i%%.*}"'.unmapped.bam';
    samtools view -b -f 4 $i > "${i%%.*}"'.unmapped.bam'
done    

# 2. get fqs:
/data2/240122_RKoenig_mouseGL261/240215_align/0246_2_TP04_001_CE_NOR_2TR1_SL1_S2Aligned.unmapped.bam
/data2/240122_RKoenig_mouseGL261/240215_align/0246_2_TP04_001_CE_NOR_1TR1_SL1_S1Aligned.unmapped.bam
/data2/240122_RKoenig_mouseGL261/240215_align/0246_2_TP04_002_CE_MOD_1TR1_SL1_S3Aligned.unmapped.bam
/data2/240122_RKoenig_mouseGL261/240215_align/0246_2_TP04_002_CE_MOD_2TR1_SL1_S4Aligned.unmapped.bam
/data2/240122_RKoenig_mouseGL261/240215_align/0246_2_TP04_003_CE_MOD_1TR1_SL1_S5Aligned.unmapped.bam
/data2/240122_RKoenig_mouseGL261/240215_align/0246_2_TP04_003_CE_MOD_2TR1_SL1_S6Aligned.unmapped.bam
/data2/240122_RKoenig_mouseGL261/240215_align/0246_2_TP04_004_CE_MOD_1TR1_SL1_S7Aligned.unmapped.bam
/data2/240122_RKoenig_mouseGL261/240215_align/0246_2_TP04_004_CE_MOD_2TR1_SL1_S8Aligned.unmapped.bam

conda activate bedtools
for i in *unmapped.bam; do
    echo $i;
    #echo "${i%%.*}"'.nameSort.unmapped.bam';
    outFile="${i%%.*}"'.nameSort.unmapped.bam';
    fq1="${i%%.*}"'.nameSort.unmapped_R1.fq';
    fq2="${i%%.*}"'.nameSort.unmapped_R2.fq'; 
    samtools sort -n -o $outFile $i;
    bamToFastq -i $outFile -fq $fq1 -fq2 $fq2
done    


test=/data2/240122_RKoenig_mouseGL261/240215_align/0246_2_TP04_002_CE_MOD_2TR1_SL1_S4Aligned.sortedByCoord.out.bam
samtools view $test | head


inFile=/data2/240122_RKoenig_mouseGL261/240215_align/0246_2_TP04_004_CE_MOD_1TR1_SL1_S7Aligned.unmapped.bam
fq1=/data2/240122_RKoenig_mouseGL261/240215_align/0246_2_TP04_004_CE_MOD_1TR1_SL1_S7Aligned.unmapped_R1.fq
fq2=/data2/240122_RKoenig_mouseGL261/240215_align/0246_2_TP04_004_CE_MOD_1TR1_SL1_S7Aligned.unmapped_R2.fq 
bamToFastq [OPTIONS] -i <BAM> -fq <FASTQ>  -fq2 <FASTQ> 

# first sort BAM by query name:
samtools sort -n -o aln.qsort.bam aln.bam
# same version in aligne and bedtools!
#       samtools 1.9
#       Using htslib 1.9



# 3. put them into qiime









# --------------- # ---------------# ---------------# ---------------# ---------------# ---------------
# Try & fail to install mouse GO db:
# ---------------
# INSTALL mouse databank: not possible in R: version/dependencies. Try mamba:y
# in env R-4.2_cpFix:
mamba install bioconductor-org.mm.eg.db

# export all packages:
conda env export > /home/hermi/Documents/conda_envs/R-4.2_cpFix.yml
# manually delete most enties, incuding versions and 
mamba env create -f /home/hermi/Documents/conda_envs/R-4.2_min.yml
#ERROR conda.core.link:_execute(740): An error occurred while installing package 
#'bioconda::bioconductor-genomeinfodbdata-1.2.9-r42hdfd78af_0'.
#Rolling back transaction: done
#class: LinkError

# try new mamba env:
#mamba create -n R-4.2_mouseDB r-base=4.2 rstudio bioconductor-org.mm.eg.db -c conda-forge -c bioconda
#mamba install r-knitr r-markdown r-dplyr r-data.table r-stringr bioconductor-clusterprofiler -c conda-forge -c bioconda

# try BiocManager in rmd:
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GenomeInfoDbData")
# nope
# --------------- # ---------------# ---------------# ---------------# ---------------# ---------------








# ----------# ----------# ----------
# ----------# ----------# ----------
# options: for all cmds here
# ----------# ----------# ----------
# ----------# ----------# ----------


# QC loop:
# --------
for i in $(ls *.fastq.gz | grep -Ev '*.*_pA_*.*fastq.gz'); do
    COND=$(basename $i .fastq.gz);
    echo $COND;
done 
# &> $log: necessary! stdout has no information!
# grep -E, --extended-regexp
#              Interpret  PATTERN  as  an extended regular expression (ERE, see
#              below).  (-E is specified by POSIX.)
#       -v, --invert-match
#              Invert the sense of matching, to select non-matching lines.  (-v
#              is specified by POSIX.)
#           > exclude matches!





# from STAR manual:
# -----------------
#   3.2.2    ENCODE options An example of ENCODE standard options 
#            for long RNA-seq pipeline is given below:
# --outFilterType BySJout
#       reduces the number of ”spurious” junctions
# --outFilterMultimapNmax 20
#       max  number  of  multiple  alignments  allowed  for  a  read:
#         if  exceeded,  the  read  is  considered unmapped
# --alignSJoverhangMin 8
#       minimum overhang for unannotated junctions
# --outFilterMismatchNmax 999
#       maximum number of mismatches per pair, large number switches off this filter
# --outFilterMismatchNoverReadLmax 0.04     # comment Csaba, not relevant for me: ######### Here is different for pA protocol !!!!!!!!!!!!!!!
#       max number of mismatches per pair relative to read length:
#         for 2x100b, max number of mis-matches is 0.04*200=8 #for the paired read
# --alignIntronMin 20
#       minimum intron length
# --alignIntronMax 1000000
#       maximum intron length
# --alignMatesGapMax1000000
#       maximum genomic distance between mates
#       # set by me, not by Csaba. There are no mates anyway..? But It shouldn't hurt?
# --limitBAMsortRAM    default: 0
#      int>=0: maximum available RAM (bytes) for sorting BAM. If =0, it will be set
#      to the genome index size. 0 value can only be used with –genomeLoad
#      NoSharedMemory option.
# --quantMode GeneCounts 
#     not used! program "featureCounts" used instead!
#       count reads per gene
# --outSAMmapqUnique 60
#   sets unique mapping value. This ensures compatibility with downstream tools such as GATK
# --twopassMode Basic 
#       Annotated junctions will be included in both the 1st and 2nd passes. To run STAR 2-pass mapping for
#       each sample separately, use --twopassMode Basic option. STAR will perform the 1st pass mapping,
#       then it will automatically extract junctions, insert them into the genome index, and, finally, re-map
#       all reads in the 2nd mapping pass. This option can be used with annotations, which can be included
#       either at the run-time (see #1), or at the genome generation step.
#       DE analysis: https://github.com/alexdobin/STAR/issues/1616: 
#       "the 2-pass approaches mainly affect the alignment of the reads to unannotated 
#       junctions, which should not have a strong effect on DGE. However, it also reduces 
#       false alignments to pseudogenes and other sequence similar regions, so it may have 
#       some effect on some genes with similar sequences. I would recommend running it
#       with and without the --twopassMode Basic to see how much difference it makes. The 
#       more complicated multi-sample 2-pass option is overkill for DGE in most cases."
# zcat = gunzip -c
# --limitBAMsortRAM 80000000000 \
#       changed from 8 GB to 80 GB! I guess Csaba wanted to put in one 0 more, too.  
#       this is pretty small. default would be 4x this!
#       genome index size = 317000000000 bytes (31,7 GB)
#       RAM on system = 125,8 GiB
# --sjdbOverhang 100 
#       only needed when annotation gtf was not included in the genome Create step, unnecessary!
# --genomeLoad LoadAndKeep \
#     removed for 2passmode
# --outSAMmapqUnique 60 \
# --twopassMode Basic 
# --alignMatesGapMax 1000000 \
#       used, it's in manual example!








