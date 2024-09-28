#!/bin/bash

###===========================
##Priming the variables used
echo -e "---0. Priming the variables used"
###===========================
sample_name=$1
ont_reads=$2
reference_fasta=$3
outdir=$4

overall_polishing_counter=1

###===========================
##Starting the for-loop
echo -e "\n-----0. Starting the for-loop"
###===========================

for ONT_polishing_counter in $(echo "first second"); do
  echo -e "\n-------"${ONT_polishing_counter}" round of Racon polishing"
  ###===========================
  ##graphmap mapping
  echo -e "--------1. First we map with graphmap"
  ###===========================

  mkdir -p $outdir/0${overall_polishing_counter}_${ONT_polishing_counter}_Racon_graphmap/{Graphmap_ONT,assembly}
  
  # Store path to output sam into a variable
  overlap=${outdir}/0${overall_polishing_counter}_${ONT_polishing_counter}_Racon_graphmap/Graphmap_ONT/${ONT_polishing_counter}_Racon_ONT2FinalAssembly_sorted.sam
  
  graphmap align --rebuild-index --circular -r ${reference_fasta} -d ${ont_reads} -o ${overlap}
  
  ###===========================
  ##racon polishing
  echo -e "--------2. Then we polish with Racon"
  ###===========================

  echo -e "Reference assembly:  " ${reference_fasta} "  of size:" $(stat --printf="%s" ${reference_fasta}) "bytes"
  echo -e "Reads for correction:  " ${ont_reads} "  of size:" $(stat --printf="%s" ${ont_reads}) "bytes"
  echo -e "Overlap set:  " ${overlap} "  of size:" $(stat --printf="%s" ${overlap}) "bytes"
  
  racon ${ont_reads} ${overlap} ${reference_fasta} > \
  ${outdir}/0${overall_polishing_counter}_${ONT_polishing_counter}_Racon_graphmap/assembly/${sample_name}_${ONT_polishing_counter}_Racon_Assembly.fasta

  ###===========================
  ##Reset the variables and reference
  echo -e "--------3. Reset the variables and reference assembly"
  ###===========================

  reference_fasta=${outdir}/0${overall_polishing_counter}_${ONT_polishing_counter}_Racon_graphmap/assembly/${sample_name}_${ONT_polishing_counter}_Racon_Assembly.fasta
  overall_polishing_counter=$((overall_polishing_counter+1))

done
