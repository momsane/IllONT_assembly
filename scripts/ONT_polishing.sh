###===========================
##Priming the variables used
echo -e "-------0. Priming the variables used"
###===========================
ont_reads=$2
reference_fasta=$3
sample_name=$1

Overall_output_directory=$4

overall_polishing_counter=1
###===========================
##Starting the for-loop
echo -e "-------0. Starting the for-loop"
###===========================

for ONT_polishing_counter in $(echo "first second"); do
  echo -e "----"${ONT_polishing_counter}" round of Racon polishing"

  ###===========================
  ##graphmap mapping
  echo -e "-------1. First we map with graphmap"
  ###===========================

  mkdir -p $Overall_output_directory/0${overall_polishing_counter}_${ONT_polishing_counter}_Racon_graphmap/{Graphmap_ONT,assembly}

  graphmap align --rebuild-index --circular  \
      -r $reference_fasta \
      -d $ont_reads -o ${Overall_output_directory}/0${overall_polishing_counter}_${ONT_polishing_counter}_Racon_graphmap/Graphmap_ONT/${ONT_polishing_counter}_Racon_ONT2FinalAssembly_sorted.sam

  ###===========================
  ##racon polishing
  echo -e "-------2. Second polish with Racon"
  ###===========================


  touch ${Overall_output_directory}/0${overall_polishing_counter}_${ONT_polishing_counter}_Racon_graphmap/assembly/${sample_name}_${ONT_polishing_counter}_Racon_Assembly.fasta
  echo -e "####### REF ${reference_fasta} ###################"
  echo -e "####### Reads ${ont_reads} ###################"
  echo -e "####### OV ${Overall_output_directory}/0${overall_polishing_counter}_${ONT_polishing_counter}_Racon_graphmap/Graphmap_ONT/${ONT_polishing_counter}_Racon_ONT2FinalAssembly_sorted.sam ###################"
  racon ${ont_reads} \
  ${Overall_output_directory}/0${overall_polishing_counter}_${ONT_polishing_counter}_Racon_graphmap/Graphmap_ONT/${ONT_polishing_counter}_Racon_ONT2FinalAssembly_sorted.sam \
  ${reference_fasta} > \
  ${Overall_output_directory}/0${overall_polishing_counter}_${ONT_polishing_counter}_Racon_graphmap/assembly/${sample_name}_${ONT_polishing_counter}_Racon_Assembly.fasta


  ###===========================
  ##Reset the variables and reference
  echo -e "-------3. Reset the variables and reference"
  ###===========================

  reference_fasta=$(echo ${Overall_output_directory}/0${overall_polishing_counter}_${ONT_polishing_counter}_Racon_graphmap/assembly/${sample_name}_${ONT_polishing_counter}_Racon_Assembly.fasta)
  overall_polishing_counter=$((overall_polishing_counter+1))

done
