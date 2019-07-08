## use MACS2 to call peaks in each sample
## this serves as a proxy for the signal-to-noise ratio

dir=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2019/ATAC-seq

bams=$(ls -d -1 $dir/data/BWA/* | grep bam)

for b in $bams
do
  sample=${b/$dir\/data\/BWA\//}
  sample=${sample/.noDUPs.GQ.bam/}
  echo Processing sample $sample
  macs2 callpeak -f BAMPE -g mm --keep-dup all --broad -t $b --outdir $dir/peaks/individualReps/ -n $sample
done


