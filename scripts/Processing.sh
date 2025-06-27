#!/bin/bash
 # Be sure to have placed the .fastq.gz files in the current folder 
 mkdir -p ./Read_counts
 
 source "$(conda info --base)/etc/profile.d/conda.sh"
 
 conda init
 conda activate metagenomics_env
# Merge the files from different lanes 
samples=$(ls reads/*.fastq.gz | awk -F'_S' '{print $1}' | sort -u)
for sample in $samples; do
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Starting concatenation for ${sample}"
    for read in R1 R2; do
        echo "$(date '+%Y-%m-%d %H:%M:%S') - Concatenating ${read} files for ${sample}..."
        cat ${sample}_*_${read}_*.fastq.gz > ${sample}_${read}.fastq.gz
    done
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Removing original chunk files for ${sample}"
	rm ${sample}_*_001.fastq.gz

    echo "$(date '+%Y-%m-%d %H:%M:%S') - Finished merging files for ${sample}"
done

# Decompress
for sample in *.gz; do
gunzip $sample
echo "$(date '+%Y-%m-%d %H:%M:%S') - ${sample} has been decompressed"
done

# Use KneadData

TRF_DIR="$(conda info --base)/envs/metagenomics_env.bin"

samples=$(for file in *.fastq; do echo "$file" | cut -d'_' -f1; done | sort -u) #
for sample in $samples; do
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Starting Kneaddata cleaning for ${sample}"
           kneaddata     --input1 "$sample"_R1.fastq \
--input2 "$sample"_R2.fastq   \
--reference-db ./human_genome  \
--output "$sample"     \
--trimmomatic ./Trimmomatic-0.39 \
--trf "$TRF_Dir"  \
--threads 10
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Finished Kneaddata for ${sample} ; formatting the output in a single file"
cat ${sample}/${sample}*paired_1.fastq ${sample}/${sample}*paired_2.fastq > ${sample}.fastq
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Finished merging the files, removing the intermediary folder"
rm -r "${sample}" 
rm -f "${sample}_R1.fastq" "${sample}_R2.fastq"
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Finished !"
done

for file in *.fastq; do
    count=$(cat "$file" | grep "^@" -c)
    filename=$(basename "$file" .fastq)
    echo "$filename: $count" >> Read_counts/read_counts.txt
done


# Application de la fonction Metaphlan pour obtenir des tables d'abondance

for i in *fastq; do 
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Starting MetaPhlAn for ${i}"
    metaphlan $i --bowtie2db ./metaphlan_database --index mpa_vJan21_CHOCOPhlAnSGB_202103 --input_type fastq --nproc 10 --add_viruses -t rel_ab_w_read_stats > ${i%.fastq}.txt; 
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Finished MetaPhlAn for ${i}"
    done

# Fusion des différents fichiers pour obtenir une table unique.
    echo "$(date '+%Y-%m-%d %H:%M:%S') - All MetaPhlAn done, merging the outputs"

rm *bowtie2*
rm *fastq*

if [ ! -s data/data_GM_Response.txt ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Generating data/data_GM_Response.txt..."
    python ./Scripts/merge_metaphlan_tables_abs.py *.txt > data/data_GM_Response.txt
else
    echo "$(date '+%Y-%m-%d %H:%M:%S') - data/data_GM_Response.txt already exists and is not empty."
fi

rm *txt*
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Merged table created !"
