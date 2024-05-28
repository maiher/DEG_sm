# outfile:
# ---------
# total reads + original ref:
statsTab=/home/AD/hermi/sm_tutorial/BYOC/testrun/240526_QCstats.tsv

# only important for virus mapping:
# leave like that wo virus
statsPref=''

fields="sId\tno_raw_reads\tno_qc_reads\tno_mapped_reads\tperc_qc_of_raw\tperc_mapped_of_qc\tav_len_mapped\tno_unmapped_reads\tno_virus_mapped_reads\tperc_virus_mapped_reads\n"
printf "$fields" > $statsTab

fqDir=/home/AD/hermi/sm_tutorial/BYOC/testrun/240526_qcFastq
bamDir=/home/AD/hermi/sm_tutorial/BYOC/testrun/240526_align



# eg Calu3p28TRIzol_1_S1_R1_001_inf1Aligned.sortedByCoord.out.bam
# loop test: get files
cd $bamDir
for log in $(cd $bamDir && ls *Log.final.out); do
    echo $log;
    prefix=$(basename $log Log.final.out);
    echo $prefix;
    # get fastp log:
    fp=$fqDir'/'$prefix'.log'
    echo $fp;
done


# loop for real:
for log in $(cd $bamDir && ls *Log.final.out); do
    echo $log;
    prefix=$(basename $log Log.final.out);
    echo $prefix;
    # get fastp log:
    fp=$fqDir'/'$prefix'.log'
    echo $fp;
    noIn=$(grep 'Number of input reads' "${log}" | cut -f 2);
    noMap=$(grep 'Uniquely mapped reads number' "${log}" | cut -f 2);
    percMap=$(grep 'Uniquely mapped reads %' "${log}" | cut -f 2 | tr -d '%');
    lenMap=$(grep 'Average mapped length' "${log}" | cut -f 2);
    noRaw=$(grep -m 1 'total reads' "${fp}" | cut -d " " -f 3);

    # no unmapped reads
    noMulitMap=$(grep 'Number of reads mapped to multiple loci' "${log}" | cut -f 2);
    noUnMap=$(echo ${noIn} | awk -v a=${noMap} -v b=${noMultimap} '{print ($0-(a+b))}');
    
    # % qc filt of raw reads:
    pFilt=$(echo ${noIn} | awk -v a=${noRaw} '{print ($0/a*100)}');
    # output:
    vals="$prefix\t$noRaw\t$noIn\t$noMap\t$pFilt\t$percMap\t$lenMap\t$noUnMap";
    printf "$vals" >> $statsTab;
    # virus mapping
    statsFile=$statsPref'/'$prefix'_HBV_AD38.stats';   
    if [ -f $statsFile ]; then 
        pMapVir=$(echo ${noMapVir} | awk -v a=${noUnMap} '{print ($0/a*100)}');
        vals="\t$noMapVir\t$pMapVir\n";
        printf "$vals" >> $statsTab;
    else
        printf "\tNA\tNA\n" >> $statsTab;
    fi;
done




    

