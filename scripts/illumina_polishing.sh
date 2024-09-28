#!/bin/bash

# module purge

# module load gcc
# module load bowtie2/2.5.1
# module load samtools/1.17
# module load pilon
# # This script doesn't work in a conda env for some reason, so I run it using the cluster softwares

threads=$1
reference_fasta=$2
R1=$3
R2=$4
outdir=$5
sample_name=$6

###===========================
##Priming the variables used and starting the for-loop
echo -e "-----0. Priming the variables used and starting the for-loop"
###===========================

overall_polishing_counter=3

for Illumina_polishing_counter in $(echo "first second third"); do

  echo -e "\n-------"${Illumina_polishing_counter}" round of pilon polishing"

  ###===========================
  ##bowtie2 mapping
  echo -e "---------1. First we map with bowtie2"
  ###===========================
  mkdir -p ${outdir}/0${overall_polishing_counter}_${Illumina_polishing_counter}_pilon_bowtie2/{bowtie2,assembly}
  bam=${outdir}/0${overall_polishing_counter}_${Illumina_polishing_counter}_pilon_bowtie2/bowtie2/${Illumina_polishing_counter}_pilon_ONT2FinalAssembly_sorted.bam

  bowtie2-build --quiet ${reference_fasta} ${reference_fasta}

  bowtie2 -1 ${R1} \
          -2 ${R2} \
          -x  ${reference_fasta} \
          --threads ${threads} \
          --local --very-sensitive-local | samtools sort -O BAM -o ${bam}

  ###--------------------
  ##index bam file
  ###--------------------
  samtools index ${bam}

  ###===========================
  ##pilon polishing
  echo -e "---------2. Then polish with pilon"
  ###===========================
  echo -e "${outdir}/0${overall_polishing_counter}_${Illumina_polishing_counter}_pilon_bowtie2/assembly/"

  pilon --genome ${reference_fasta} \
  --frags ${bam} \
  --output ${sample_name}_${Illumina_polishing_counter}_Pilon_Assembly \
  --outdir ${outdir}/0${overall_polishing_counter}_${Illumina_polishing_counter}_pilon_bowtie2/assembly/ \
  --changes

  tree ${outdir}/0${overall_polishing_counter}_${Illumina_polishing_counter}_pilon_bowtie2/

  ###===========================
  ##reset the parameters
  echo -e "---------3. Resetting the parameters"
  ###===========================

  reference_fasta=$(echo ${outdir}/0${overall_polishing_counter}_${Illumina_polishing_counter}_pilon_bowtie2/assembly/${sample_name}_${Illumina_polishing_counter}_Pilon_Assembly.fasta)
  echo ${reference_fasta}
  overall_polishing_counter=$((overall_polishing_counter+1))

done
#polishing round
