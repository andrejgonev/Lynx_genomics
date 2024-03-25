for yaml in $(ls config/alignment_2/); do  

    sampleid=$(echo $yaml | sed 's/.mLynLyn.alignment.yaml$//')    
    echo $sampleid
    sbatch -o logs/alignment/${sampleid}.out -e logs/alignment/${sampleid}.err scripts/sbatch_alignment_of_sample_from_yaml.sh config/alignment_2/$yaml

done