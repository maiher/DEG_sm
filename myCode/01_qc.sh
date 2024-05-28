# test fastqs: troubleshooting
# see for ref:
#/data2/221114_RBrown_PBMC



#-------------------------#
#   fastp QC              #
#-------------------------#
#
#
#    -5, --cut_front                      move a sliding window from front (5') to tail, drop the bases in the window if its mean quality < threshold, stop otherwise.
#    -3, --cut_tail                       move a sliding window from tail (3') to front, drop the bases in the window if its mean quality < threshold, stop otherwise.
#    -r, --cut_right                      move a sliding window from front to tail, if meet one window with mean quality < threshold, drop the bases in the window and the right part, and then stop.
#   If --cut_right is enabled, then there is no need to enable --cut_tail, since the former is more aggressive. 
#   If --cut_right is enabled together with --cut_front, --cut_front will be performed first before --cut_right to avoid dropping whole reads due to the low quality starting bases.
#    -G:            disable polyG tail trimming, by default trimming is automatically enabled for
#                   Illumina NextS
#                   legacy prompt from Csaba, not sure why used?
#    -V:            be vervose (output log)   
#    -w, --thread                         worker thread number, default is 2 (int [=2])
#
#
#-------------------------#

# WD for all cmds:
baseDir=/home/AD/hermi/sm_tutorial/BYOC/testrun
# store qc fqs at:
outDir=$baseDir"/240526_qcFastq"
mkdir -p $outDir
# raw fqs:
fqDir=/home/AD/hermi/sm_tutorial/BYOC/myData
cd $fqDir

conda activate DEGpipe01
#_R1_001_50k_fastq.gz
remove_right='_R1*'
query="*_R1_001_50k_fastq.gz"

for i in $query; do
    sampleID=${i//$remove_right/};
    r1=$fqDir'/'$i;
    r2=$fqDir'/'$sampleID'_R2_001_50k_fastq.gz';
    #echo $r1;
    #echo $r2;
    #/data2/240409_RKoenig_mouseGL261_round2/240409_rawFastq/0246_2_TP04_026_FC_MOD_1TR1_SL1_S26_L002_R1_001.fastq.gz
    #/data2/240409_RKoenig_mouseGL261_round2/240409_rawFastq/0246_2_TP04_026_FC_MOD_1TR1_SL1_S26_L002_R2_001.fastq.gz
    #/0246_2_TP04_005_FC_MOD_1TR1_SL1_S9_L001_R1_001-pooled.fastq.gz
    #/0246_2_TP04_005_FC_MOD_1TR1_SL1_S9_L001_R2_001-pooled.fastq.gz
    COND=$(basename $i _L002_R1_001_50k_fastq.gz);
    echo $COND;
    #0246_2_TP04_005_FC_MOD_1TR1_SL1_S9_
    fq1=$outDir'/'$COND'_R1_qc.fq.gz'
    fq2=$outDir'/'$COND'_R2_qc.fq.gz'
    #echo $fq1;
    #echo $fq2;
    #printf "# QC on separate R1 and R2 fastqs!\n";
    time fastp \
    -i $r1 \
    -I $r2 \
    -o $fq1 \
    -O $fq2 \
    --adapter_sequence=TCAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
    -3 \
    -5 \
    -G \
    -V \
    -w 16 \
    &> $outDir/$COND'.log';
    mv fastp.html $outDir/$COND'.html';
done
rm fastp.json

