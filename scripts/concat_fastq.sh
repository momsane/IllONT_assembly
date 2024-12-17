#!/bin/bash

## This script takes an input a tab-delimited table reporting for each strain, its sequencing barcode (first column) and its name (2nd column)
## The table does not have a header. See config/strains_barcodes.txt for an example.
## The inputs files are organized into one folder per barcode
## Then the script concatenates all reads files for each barcode/strain then names the output with the strain instead of barcode.

# Define inputs

# Path to dir where raw fastq files are stored
dir=/Volumes/RECHERCHE/FAC/FBM/DMF/pengel/general_data/D2c/datasets/NGS_data/20240730_syncom_Meline_ONT/20240815_X3_ES171_IFIK/fastq_pass

# Path to desired output dir
outdir=/Volumes/RECHERCHE/FAC/FBM/DMF/pengel/general_data/D2c/datasets/NGS_data/20240730_syncom_Meline_ONT/20240815_X3_ES171_IFIK/concatenated_reads

# Path to table with strain number and barcode
table=/Volumes/RECHERCHE/FAC/FBM/DMF/pengel/general_data/D2c/datasets/NGS_data/20240730_syncom_Meline_ONT/strains_barcodes.txt

cd ${dir}
mkdir -p ${outdir}

for bc in $(cat $table | cut -f 1); do
    barcode=$(echo "barcode$bc") 
    echo "processing $barcode"
    file=${outdir}/$(awk -v bc=$bc '$1 == bc {print $2}' "$table")_concatenated.fastq.gz # make output file name with strain number
    cat ${barcode}/* > "$file" # concatenate all fastq files present in the barcode folder
    mv ${file} "$(echo $file | sed s'/\r//g')" # rename file to remove ? character, don't know why it appears
done



