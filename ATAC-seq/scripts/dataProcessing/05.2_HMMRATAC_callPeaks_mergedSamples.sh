## Call peaks with HMMCRATAC

dir=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2019/ATAC-seq

java -jar /user01/group_folders/Personal/Ximena/software/HMMRATAC_V1.2.5_exe.jar -b $dir/HMMRATAC/mergedSamples.bam -i $dir/HMMRATAC/mergedSamples.bam.bai -g $dir/HMMRATAC/genome.info -o $dir/HMMRATAC/allGQsamples_HMMRATAC -e $dir/data/mm10-blacklist.v2.bed -p True --bedgraph True --window 2500000 --modelonly True
