#!/bin/bash
 # Be sure to have placed the .fastq.gz files in the current folder 


 source "$(conda info --base)/etc/profile.d/conda.sh" 
 conda activate metagenomics_env

# Merge the files from different lanes 

samples=$(ls reads/*.fastq.gz | xargs -n1 basename | awk -F'_S' '{print $1}' | sort -u)

for sample in $samples; do
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Starting concatenation for ${sample}"
    
    for read in R1 R2; do
        echo "$(date '+%Y-%m-%d %H:%M:%S') - Concatenating ${read} files for ${sample}..."

        cat reads/${sample}_S*_L*_${read}_001.fastq.gz > reads/${sample}_${read}.fastq.gz
    done

    echo "$(date '+%Y-%m-%d %H:%M:%S') - Removing original chunk files for ${sample}"
    rm reads/${sample}_S*_L*_*_001.fastq.gz

    echo "$(date '+%Y-%m-%d %H:%M:%S') - Finished merging files for ${sample}"
done

# Decompress
for sample in reads/*.gz; do
gunzip $sample
echo "$(date '+%Y-%m-%d %H:%M:%S') - ${sample} has been decompressed"
done

# Use KneadData

TRF_DIR="$CONDA_PREFIX/bin"

samples=$(for file in reads/*.fastq; do basename "$file" | cut -d'_' -f1; done | sort -u)

for sample in $samples; do
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Starting Kneaddata cleaning for ${sample}"
    kneaddata \
        --input1 reads/"${sample}"_R1.fastq \
        --input2 reads/"${sample}"_R2.fastq \
        --reference-db ./human_genome \
        --output "${sample}" \
        --trimmomatic ./Trimmomatic-0.39 \
        --trf "$TRF_DIR" \
        --threads 10

    echo "$(date '+%Y-%m-%d %H:%M:%S') - Finished Kneaddata for ${sample} ; formatting the output in a single file"

    cat "${sample}"/*paired_1.fastq "${sample}"/*paired_2.fastq > "${sample}.fastq"

    echo "$(date '+%Y-%m-%d %H:%M:%S') - Finished merging the files, removing the intermediary folder"

rm -r "${sample}"
done

mkdir -p ./Read_counts 

for file in *.fastq; do
    count=$(grep -c "^@" "$file")
    filename=$(basename "$file" .fastq)
    echo "$filename: $count" >> Read_counts/read_counts.txt
done

# MetaPhlan call 

for i in *fastq; do 
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Starting MetaPhlAn for ${i}"
    metaphlan $i --bowtie2db ./metaphlan_database --index mpa_vJan21_CHOCOPhlAnSGB_202103 --input_type fastq --nproc 10 --add_viruses -t rel_ab_w_read_stats > ${i%.fastq}.txt; 
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Finished MetaPhlAn for ${i}"
    done

# MetaPhlAn merging
    echo "$(date '+%Y-%m-%d %H:%M:%S') - All MetaPhlAn done, merging the outputs"

rm *bowtie2*
rm *fastq*
rm reads/*fastq*

if [ ! -s data/data_GM_Response.txt ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Generating data/data_GM_Response.txt..."
    python ./scripts/merge_metaphlan_tables_abs.py *.txt > data/data_GM_Response.txt
else
    echo "$(date '+%Y-%m-%d %H:%M:%S') - data/data_GM_Response.txt already exists and is not empty."
fi

rm *txt*
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Merged table created !"
