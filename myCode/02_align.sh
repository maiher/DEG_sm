# alignment: STAR
# ----------------
# see doku on used options below.
conda activate DEGpipe01

baseDir=/home/AD/hermi/sm_tutorial/BYOC/testrun
# store bams at:
outDir=$baseDir"/240526_align"
mkdir $outDir
cd $outDir

# GRCm39 mouse genome:
# belongs to Csaba!
# parameters: /data/public/RefSeq.GRCm39/STARindexes/genomeParameters.txt
#  same as mine in /data/public/RefSeq.GRCh38.p14/index/STAR.2.7.1a/genomeParameters.txt!
#genDir=/data/public/RefSeq.GRCm39/STARindexes
genDir=/data/public/RefSeq.GRCm39/STARindexes

fqDir=/home/AD/hermi/sm_tutorial/BYOC/testrun/240526_qcFastq

# all samples:
query='*_R1_qc.fq.gz'
remove_right='_R1*'

# loop: test & adjust
for i in $(cd $fqDir && ls $query); do
    sampleID=${i//$remove_right/};
    fq1=$fqDir'/'$i;
    fq2=$fqDir'/'$sampleID'_R2_qc.fq.gz';
    echo $fq1;
    echo $fq2;
    echo $sampleID;
    COND=$sampleID
done
#/data2/240122_RKoenig_mouseGL261/240125_qcFastq/0246_2_TP04_004_CE_MOD_2TR1_SL1_S8_R1_qc.fq.gz
#/data2/240122_RKoenig_mouseGL261/240125_qcFastq/0246_2_TP04_004_CE_MOD_2TR1_SL1_S8_R2_qc.fq.gz
#0246_2_TP04_004_CE_MOD_2TR1_SL1_S8

# loop: final
for i in $(cd $fqDir && ls $query); do
    sampleID=${i//$remove_right/};
    fq1=$fqDir'/'$i;
    fq2=$fqDir'/'$sampleID'_R2_qc.fq.gz';
    echo $fq1;
    #echo $fq2;
    #echo $sampleID;
    log=$sampleID'.log';
    time STAR \
    --runThreadN 8 \
    --readFilesIn $fq1 $fq2 \
    --outFileNamePrefix $sampleID \
    --genomeDir $genDir \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outFilterType BySJout \
    --outFilterMultimapNmax 20 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.6 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --limitBAMsortRAM 800000000000 \
    --outSAMmapqUnique 60 \
    --twopassMode Basic \
    &> $outDir/$sampleID'.log';
done


#OLD:



#/data2/240122_RKoenig_mouseGL261/240125_qcFastq/0246_2_TP04_001_CE_NOR_1TR1_SL1_S1_R1_qc.fq.gz
#real	128m45.630s

#/data2/240122_RKoenig_mouseGL261/240125_qcFastq/0246_2_TP04_001_CE_NOR_2TR1_SL1_S2_R1_qc.fq.gz
#real	105m2.858s

#/data2/240122_RKoenig_mouseGL261/240125_qcFastq/0246_2_TP04_002_CE_MOD_1TR1_SL1_S3_R1_qc.fq.gz
#real	133m23.496s

#/data2/240122_RKoenig_mouseGL261/240125_qcFastq/0246_2_TP04_002_CE_MOD_2TR1_SL1_S4_R1_qc.fq.gz
#real	74m26.267s

#/data2/240122_RKoenig_mouseGL261/240125_qcFastq/0246_2_TP04_003_CE_MOD_1TR1_SL1_S5_R1_qc.fq.gz
#real	97m29.737s

#/data2/240122_RKoenig_mouseGL261/240125_qcFastq/0246_2_TP04_003_CE_MOD_2TR1_SL1_S6_R1_qc.fq.gz
#real	95m6.641s

#/data2/240122_RKoenig_mouseGL261/240125_qcFastq/0246_2_TP04_004_CE_MOD_1TR1_SL1_S7_R1_qc.fq.gz
#real	133m39.011s

#/data2/240122_RKoenig_mouseGL261/240125_qcFastq/0246_2_TP04_004_CE_MOD_2TR1_SL1_S8_R1_qc.fq.gz
#real	108m6.235s






# ----------# ----------# ----------
# ----------# ----------# ----------
# options:
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
#       - not used! program "featureCounts" used instead!
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








