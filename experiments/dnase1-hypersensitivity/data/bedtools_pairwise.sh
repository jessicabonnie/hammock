#  file_labels=`ls *.bed | sed -e 's/.hotspot.twopass.fdr0.05.merge.bed//g' \
#                                        -e 's/.hg19//g'`
#     echo name" "$file_labels >> bedtools_pairwise_jaccard.txt

echo "file1 file2 intersection union jaccard n_intersections" > bedtools_pairwise_jaccard.txt
    for file1 in `ls *.bed`
    do
        # make reasonable file labels
        # file1_short=`echo $file1 \
                    # | sed -e 's/.hotspot.twopass.fdr0.05.merge.bed//g' \
                    # -e 's/.hg19//g'`
        # echo -n $file1_short >> bedtools_pairwise_jaccard.txt

        for file2 in `ls *.bed`;
        do
            # compute the jaccard stat for these two files.
            bvalues=`bedtools jaccard \
                       -a $file1 \
                       -b $file2 \
                       | awk 'NR==2 {print $0}'`


            echo $file1 $file2 $bvalues >> bedtools_pairwise_jaccard.txt
            
            # # report the jaccard stat for these two files
            # echo -n " "$jaccard >> bedtools_pairwise_jaccard.txt
        done
        echo >> pairwise_jaccard.txt
    done