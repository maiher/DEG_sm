# featureCounts: get gene len and count data
# -------------------------------------------------------------
# man:
#/home/hermi/Documents/prog_manuals/221101_featureCounts.txt
# Usage: featureCounts [options] -a <annotation_file> -o <output_file> input_file1 [input_file2] ... 

# custom list of bams: space- or newline+' \'-separated list of paths
#bamList="/path/1/file.bam \
#/path/2/file.bam \
#/path/3/file.bam"


baseDir=/home/AD/hermi/sm_tutorial/BYOC/testrun


# list of bams
bamDir=/home/AD/hermi/sm_tutorial/BYOC/testrun/240526_align
for file in $(cd $bamDir && ls *.bam); do
    echo "$bamDir/$file"' \';
done

# store as variable:
bamList=$(for file in $(cd $bamDir && ls *.bam); do
            echo "$bamDir/$file ";
          done)

# test:
for f in $bamList; do
    echo $f;
done    

outDir=$baseDir"/240526_featureCount"
mkdir $outDir
cd $outDir
# gtf from Csabas STAR genome index:
gtf=/data/public/RefSeq.GRCm39/GCF_000001635.27_GRCm39_genomic.gtf.gz
prefix=mouseGL261_readCounts

# featureCounts cmd: 
#   1. mit -O: overlapping meta features: intron etc overlapping reads are counted too!
# not usefull right now!
#conda activate DGE
#time featureCounts \
#-p \
#-O \
#-C \
#-a $gtf -o $prefix'_O.csv' \
#$bamList &> $prefix'_O.log'
#real	136m40.658s

#   2. without "-O: overlapping meta features": only exon-mapping reads are counted
conda activate DEGpipe01
time featureCounts \
-p \
-C \
-a $gtf -o $prefix'.csv' \
$bamList &> $prefix'.log'

# real	88m2.469s







# see
https://subread.sourceforge.net/SubreadUsersGuide.pdf
#Generally, we recommend that featureCounts be run with default settings, which does not include --primary. Just follow the documentation. The specific use depends on your intended application.
# featureCounts automatically detects the format of input read files (SAM/BAM). (but not paired/single end)
#   > didn't with my files... wo -p option.

# Parameters specific to paired end reads
# -------------------------------------------------------------
#
#  -p                  If specified, fragments (or templates) will be counted
#                      instead of reads. This option is only applicable for
#                      paired-end reads; single-end reads are always counted as
#                      reads.
#
#  -B                  Only count read pairs that have both ends aligned.
#
#  -P                  Check validity of paired-end distance when counting read 
#                      pairs. Use -d and -D to set thresholds.
#
#  -d <int>            Minimum fragment/template length, 50 by default.
#
#  -D <int>            Maximum fragment/template length, 600 by default.
#
#  -C                  Do not count read pairs that have their two ends mapping 
#                      to different chromosomes or mapping to same chromosome 
#                      but on different strands.
#
#  --donotsort         Do not sort reads in BAM/SAM input. Note that reads from 
#                      the same pair are required to be located next to each 
#                      other in the input.
# Overlap between reads and features
#
#  -O                  Assign reads to all their overlapping meta-features (or 
#                      features if -f is specified).
# A read is said to overlap a feature if at least one read base is found to overlap the feature. 
#   For paired-end data, a fragment (or template) is said to overlap a feature if any of the two 
#   reads from that fragment is found to overlap the feature.
# By default, featureCounts does not count reads overlapping with more than one feature (or more than one meta-feature when summarizing at meta-feature level). Users can use the -O option to instruct featureCounts to count such reads (they will be assigned to all their overlapping features or meta-features).
# Note that, when counting at the meta-feature level, reads that overlap multiple features of the same meta-feature are always counted exactly once for that meta-feature, provided there is no overlap with any other meta-feature. For example, an exon-spanning read will be counted only once for the corresponding gene even if it overlaps with more than one exon.
I



