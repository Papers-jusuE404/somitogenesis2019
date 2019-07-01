## Merge ATAC-seq aligned reads from all samples for peak calling with HMMRATAC

dir=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2019/ATAC-seq

echo "merge..."
samtools merge -@ 14 $dir/HMMRATAC/mergedSamples.bam $dir/data/BWA/e1_SI-2.noDUPs.GQ.bam $dir/data/BWA/e1_SII-2.noDUPs.GQ.bam $dir/data/BWA/e2_SII-2.noDUPs.GQ.bam $dir/data/BWA/e3_SII-2.noDUPs.GQ.bam $dir/data/BWA/e5_SI-2.noDUPs.GQ.bam $dir/data/BWA/e9_SI-2.noDUPs.GQ.bam $dir/data/BWA/e9_SIII-2.noDUPs.GQ.bam $dir/data/BWA/e13_SI-2.noDUPs.GQ.bam $dir/data/BWA/e13_SII-2.noDUPs.GQ.bam $dir/data/BWA/e13_SIII-2.noDUPs.GQ.bam $dir/data/BWA/e14_SII-2.noDUPs.GQ.bam $dir/data/BWA/e15_SI-2.noDUPs.GQ.bam $dir/data/BWA/e15_SII-2.noDUPs.GQ.bam $dir/data/BWA/e16_SII-2.noDUPs.GQ.bam $dir/data/BWA/e16_SIII-2.noDUPs.GQ.bam $dir/data/BWA/e17_SI-2.noDUPs.GQ.bam $dir/data/BWA/e17_SII-2.noDUPs.GQ.bam $dir/data/BWA/e17_SIII-2.noDUPs.GQ.bam $dir/data/BWA/e19_SI-2.noDUPs.GQ.bam $dir/data/BWA/e19_SII-2.noDUPs.GQ.bam $dir/data/BWA/e20_SI-2.noDUPs.GQ.bam $dir/data/BWA/e20_SII-2.noDUPs.GQ.bam $dir/data/BWA/e20_SIII-2.noDUPs.GQ.bam $dir/data/BWA/e21_SI-2.noDUPs.GQ.bam $dir/data/BWA/e21_SII-2.noDUPs.GQ.bam $dir/data/BWA/e21_SIII-2.noDUPs.GQ.bam $dir/data/BWA/e22_SI-2.noDUPs.GQ.bam $dir/data/BWA/e22_SII-2.noDUPs.GQ.bam $dir/data/BWA/e22_SIII-2.noDUPs.GQ.bam $dir/data/BWA/e23_SI-2.noDUPs.GQ.bam $dir/data/BWA/e24_SI-2.noDUPs.GQ.bam $dir/data/BWA/e24_SII-2.noDUPs.GQ.bam $dir/data/BWA/e24_SIII-2.noDUPs.GQ.bam $dir/data/BWA/e25_SI-2.noDUPs.GQ.bam $dir/data/BWA/e26_SI-2.noDUPs.GQ.bam $dir/data/BWA/e26_SII-2.noDUPs.GQ.bam $dir/data/BWA/e26_SIII-2.noDUPs.GQ.bam $dir/data/BWA/e27_SI-2.noDUPs.GQ.bam $dir/data/BWA/e27_SII-2.noDUPs.GQ.bam $dir/data/BWA/e27_SIII-2.noDUPs.GQ.bam $dir/data/BWA/e28_SII-2.noDUPs.GQ.bam $dir/data/BWA/e28_SIII-2.noDUPs.GQ.bam $dir/data/BWA/e29_SII-2.noDUPs.GQ.bam $dir/data/BWA/e30_SI-2.noDUPs.GQ.bam $dir/data/BWA/e30_SII-2.noDUPs.GQ.bam $dir/data/BWA/e31_SI-2.noDUPs.GQ.bam $dir/data/BWA/e32_SI-2.noDUPs.GQ.bam $dir/data/BWA/e32_SII-2.noDUPs.GQ.bam $dir/data/BWA/e32_SIII-2.noDUPs.GQ.bam $dir/data/BWA/e33_SI-2.noDUPs.GQ.bam

echo "index..."
samtools index $dir/HMMRATAC/mergedSamples.bam
