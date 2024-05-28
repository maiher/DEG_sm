# raw fastqs - DONE


sinkDir=/data2/240409_RKoenig_mouseGL261_round2/240409_rawFastq
mkdir $sinkDir
# manually copied files   


# get md5sums:
outFile=/data2/240409_RKoenig_mouseGL261_round2/240409_md5_fastqs.txt
for i in $(find $sinkDir -type f); do
    echo $i;
    md5sum $i >> $outFile
done    


