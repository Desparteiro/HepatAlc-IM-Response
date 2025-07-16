#!/bin/bash
set -e
#The whole script takes several hours, mainly due to MetaPhlAn database installation.

# Ensure conda environment is available
source "$(conda info --base)/etc/profile.d/conda.sh"

# Create necessary folders
echo "$(date '+%Y-%m-%d %H:%M:%S') - Creating folder and environment"
mkdir -p human_genome metaphlan_database

# Create a new conda environment
conda create --name metagenomics_env -y -c bioconda python=3.7
conda activate metagenomics_env

# Install core packages
conda install -y -c bioconda metaphlan=4.0.1 kneaddata=0.12.2 trf=4.10.0

# Install Trimmomatic manually (due to issues with Bioconda package)
echo "$(date '+%Y-%m-%d %H:%M:%S') - Installing trimmomatic"
wget https://github.com/usadellab/Trimmomatic/files/5854859/Trimmomatic-0.39.zip -O Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip 

# Download the human genome database for kneaddata
echo "$(date '+%Y-%m-%d %H:%M:%S') - Downloading the human genome"
wget http://huttenhower.sph.harvard.edu/kneadData_databases/Homo_sapiens_hg37_and_human_contamination_Bowtie2_v0.1.tar.gz -O hg37.tar.gz
tar -xvzf hg37.tar.gz -C human_genome

# Download the MetaPhlAn database
echo "$(date '+%Y-%m-%d %H:%M:%S') - Downloading MetaPhlAn database"

metaphlan --install --index mpa_vJan21_CHOCOPhlAnSGB_202103 --bowtie2db ./metaphlan_database

rm hg37.tar.gz
rm Trimmomatic-0.39.zip

# Download and prepare the MetaPhlAn merging script

curl -L -o ./scripts/merge_metaphlan_tables_abs.txt https://forum.biobakery.org/uploads/short-url/3PBFnEdw1MuvQDuSh7kx9uSViGp.txt

mv scripts/merge_metaphlan_tables_abs.txt scripts/merge_metaphlan_tables_abs.py

chmod +x scripts/merge_metaphlan_tables_abs.py

