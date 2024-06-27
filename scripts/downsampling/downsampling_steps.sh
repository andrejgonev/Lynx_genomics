####################### Part 1

for input_bam in $(ls /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_bams/downsampling/*.bam); do 
  job_id=$(sbatch -c 12 --mem=10GB -t 00:15:00 /mnt/netapp2/Store_csebdjgl/agonev/scripts/mean_cov_mosdepth.sh ${input_bam} /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_bams/downsampling/mosdepth | awk '{print $4}')
    echo "${job_id} ${input_bam}" >> /mnt/netapp2/Store_csebdjgl/agonev/logs/mosdepth/downsampling/job_ids_mean_cov_mosdepth.txt
done

### for the New balkan lynx samples
for input_bam in $(ls /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_bams/blx_new/*.bam | grep -E 'c_ll_ba_0263|c_ll_ba_0264|c_ll_ba_0266'); do 
  job_id=$(sbatch -c 12 --mem=10GB -t 00:15:00 /mnt/netapp2/Store_csebdjgl/agonev/scripts/mean_cov_mosdepth.sh ${input_bam} /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_bams/blx_new/mosdepth | awk '{print $4}')
    echo "${job_id} ${input_bam}" >> logs/mosdepth/job_ids_mean_cov_mosdepth.txt
done



### Part 2. get a table with the mean coverage per sample and the downsampling factor for each of the target coverages.

# downsample mosdepth coverage to 6, 8, and 10x

for input_mosdepth in $(ls /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_bams/mosdepth/*.summary.txt); do
  cov=$(grep -w "total" ${input_mosdepth} | cut -f4);
  sample=$(basename ${input_mosdepth} _mLynLyn_ref.sorted.rg.merged_sorted.rmdup.indelrealigner_mosdepth.summary.txt)
  target_6=$(LC_ALL=C printf "%.4f" $(echo "scale=4; 6/${cov}" | bc -l))
  target_8=$(LC_ALL=C printf "%.4f" $(echo "scale=4; 8/${cov}" | bc -l))
  target_10=$(LC_ALL=C printf "%.4f" $(echo "scale=4; 10/${cov}" | bc -l))
  echo -e "${sample}\t${cov}\t${target_6}\t${target_8}\t${target_10}" >> /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_bams/mosdepth/downsampling_factor_table.txt
done

### Balkan lynx samples

for input_mosdepth in $(ls /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_bams/blx_new/mosdepth/*.summary.txt); do
  cov=$(grep -w "total" ${input_mosdepth} | cut -f4);
  sample=$(basename ${input_mosdepth} _mLynLyn_ref.sorted.rg.merged_sorted.rmdup.indelrealigner_mosdepth.summary.txt)
  target_6=$(LC_ALL=C printf "%.4f" $(echo "scale=4; 6/${cov}" | bc -l))
  target_8=$(LC_ALL=C printf "%.4f" $(echo "scale=4; 8/${cov}" | bc -l))
  target_10=$(LC_ALL=C printf "%.4f" $(echo "scale=4; 10/${cov}" | bc -l))
  echo -e "${sample}\t${cov}\t${target_6}\t${target_8}\t${target_10}" >> /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_bams/blx_new/mosdepth/downsampling_factor_table.txt
done


####################### Downsampling 

#  list all bam files that need to be downsampled. They all contain the above listed strings

ls /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_bams/*ner.bam | grep -E 'c_ll_ba_0224|c_ll_ba_0226|c_ll_ba_0227|c_ll_ba_0228|c_ll_ba_0229|c_ll_ba_0230|c_ll_ba_0233|c_ll_ca_0240|c_ll_ca_0241|c_ll_ca_0242|c_ll_ca_0243|c_ll_ca_0244|c_ll_ca_0245|c_ll_ca_0247|c_ll_ca_0248|c_ll_ca_0252|c_ll_ca_0254|c_ll_ca_0259|c_ll_ca_0260|c_ll_cr_0205|c_ll_cr_0206|c_ll_cr_0207|c_ll_cr_0208|c_ll_cr_0209|c_ll_cr_0211|c_ll_cr_0212|c_ll_ki_0090|c_ll_no_0065|c_ll_po_0150' > bams_downsample.txt

#  run the downsampling script for each bam file
for input_bam in $(cat /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_bams/downsampling/bams_downsample.txt); do 
    job_id=$(sbatch -c 5 --mem=10GB -t 00:45:00 /mnt/netapp2/Store_csebdjgl/agonev/scripts/downsampling_samtools.sh ${input_bam} /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_bams/mosdepth/downsampling_factor_table.txt | awk '{print $4}')
      echo "${job_id} ${input_bam}" >> /mnt/netapp2/Store_csebdjgl/agonev/logs/downsampling/job_ids_downsampling_samtools.txt
  done
  
## The same for the three new Balkan lynx samples

mkdir -p /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_bams/blx_new/downsampling

ls /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_bams/blx_new/*ner.bam | grep -E 'c_ll_ba_0263|c_ll_ba_0264|c_ll_ba_0266' > /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_bams/blx_new/downsampling/blx_new_bams_downsample.txt

#  run the downsampling script for each bam file
for input_bam in $(cat /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_bams/blx_new/downsampling/blx_new_bams_downsample.txt); do 
    job_id=$(sbatch -c 5 --mem=10GB -t 00:45:00 /mnt/netapp2/Store_csebdjgl/agonev/scripts/downsampling_samtools.sh ${input_bam} /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_bams/blx_new/downsampling/downsampling_factor_table.txt | awk '{print $4}')
      echo "${job_id} ${input_bam}" >> /mnt/netapp2/Store_csebdjgl/agonev/logs/downsampling/job_ids_downsampling_samtools.txt
  done




####################### Part 3 Check coverage again after downsampling

for input_bam in $(ls /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_bams/blx_new/downsampling/*.bam); do 
  job_id=$(sbatch -c 5 --mem8GB -t 00:10:00 /home/csic/eye/lmf/scripts/Phasing_and_imputation/ref_panel_validation/mean_cov_mosdepth.sh ${input_bam} /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/novogene_lp_sept23/downsampling/mosdepth | awk '{print $4}')
    echo "${job_id} ${input_bam}" >> /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/mosdepth/job_ids_mean_cov_mosdepth.txt
done

mkdir -p /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_bams/blx_new/downsampling/mosdepth

for input_bam in $(ls /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_bams/blx_new/downsampling/*.bam); do 
  job_id=$(sbatch -c 12 --mem=10GB -t 00:15:00 /mnt/netapp2/Store_csebdjgl/agonev/scripts/mean_cov_mosdepth.sh ${input_bam} /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_bams/blx_new/downsampling/mosdepth | awk '{print $4}')
    echo "${job_id} ${input_bam}" >> logs/mosdepth/downsampling/job_ids_mean_cov_mosdepth.txt
done



## get just the coverage (after downsampling)
for input_mosdepth in $(ls /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_bams/blx_new/downsampling/mosdepth/*.summary.txt); do
  cov=$(grep -w "total" ${input_mosdepth} | cut -f4);
  sample=$(basename ${input_mosdepth} mLynLyn_ref.sorted.rg.merged_sorted.rmdup.indelrealigner*)
  echo -e "${sample}\t${cov}" >> /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_bams/downsampling/mosdepth/downsampling_coverage_table.txt
done



#### Part 4: Qualimap of all samples

# bam folder
bam_folder=/mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_bams/downsampling

# sbatch QualiMap of all samples
bams=$(ls ${bam_folder}/*.bam | grep -E 'c_ll_ba_0263|c_ll_ba_0264|c_ll_ba_0266')

for bam in $bams ; do

    sampleid=$(basename "$bam" | cut -d'_' -f1,2,3,4)
    echo "sbatch qualimap of ${sampleid}"
    sbatch -o logs/downsampling/qualimap/qualimap-${sampleid}.out -e logs/downsampling/qualimap/qualimap-${sampleid}.err \
    scripts/sbatch_qualimap_bam_out.sh ${bam} data/downsampling/qualimap

done

