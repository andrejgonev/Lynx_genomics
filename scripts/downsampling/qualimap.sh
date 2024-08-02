# bam folder
bam_folder=/mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_bams/blx_new

# sbatch QualiMap of all samples
bams=$(ls ${bam_folder}/*er.bam)

for bam in $bams ; do

   sampleid=$(basename "$bam" | sed 's/_mLynLyn_ref.sorted.rg.merged_sorted.rmdup.indelrealigner.bam$//')
   echo "sbatch qualimap of ${sampleid}"
   sbatch -o logs/alignment/qualimap/qualimap-${sampleid}.out -e logs/alignment/qualimap/qualimap-${sampleid}.err \
   scripts/sbatch_qualimap_bam_out.sh ${bam} data/qualimap/blx_new

done