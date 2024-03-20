# I included a batch size of 50 because that's the maximum I can submit to cesga. The counter starts at 0 and increases until it reaches 50, then it resets to 0 and waits for the 50 jobs to finish. This is a good way to avoid overloading the cluster with too many jobs at once. 
# The script is called with the following command:
# I didn't run this script, because I was not sure if it would work. I ran the sbatch command instead.


batch_size=50
counter=0

for yaml in $(ls config/alignment/); do

    sampleid=$(echo $yaml | sed 's/.mLynLyn.alignment.yaml$//')
    echo $sampleid
    sbatch -o logs/alignment/${sampleid}.out -e logs/alignment/${sampleid}.err scripts/sbatch_alignment_of_sample_from_yaml.sh config/alignment/$yaml

    counter=$((counter + 1))

    if [ $counter -eq $batch_size ]; then
        counter=0
        wait
    fi

done