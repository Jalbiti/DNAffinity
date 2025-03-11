#!/bin/bash

script="selex.R"

FASTQ_FOLDER="fastqs"
OUTPUT_FOLDER="selex_results"

###############################################################################

# Create Config Files for all proteins:
for d in $FASTQ_FOLDER/*/ ; do
	
	protein=$(echo $d | cut -d "/" -f 2)
	filename=config.xml
	touch $filename

	#Copy R0 equally for all proteins
	sed "s/protein/$protein/" initial.txt >> $filename
	
	#Append Rounds depending on length
	j=0
	for file in `ls --color=no $FASTQ_FOLDER/$protein/*`; do
		if [ $j != 0 ] 
		then
			sed -e "s/protein/$protein/" -e "s/NUM/$j/" -e "s/NUM2/$j/" sample.txt >> $filename
	        fi
		
		((j=j+1))
	
	done
	echo "</SELEXSequencingConfig>" >> $filename
	mv $filename $FASTQ_FOLDER/$protein
done

###############################################################################

for d in $FASTQ_FOLDER/*/ ; do

   protein=$(echo $d | cut -d "/" -f 2)

   input=$FASTQ_FOLDER/$protein
   output=$OUTPUT_FOLDER/$protein

   mkdir $output
   out_file=$output/R_stats.txt

Rscript $script $protein $input $output > $out_file

done;



