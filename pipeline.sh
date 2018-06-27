#!/bin/bash

###
### DESCRIBES ROOT PIPELINE FOR GENERATING ARGS AND PAIRWISE DISTANCE MATRICES
###

#### CHECK FOR VALID INPUT AND PREPARE INPUT FILES FOR ARGWEAVER

# check if correct parameters are supplied
if [ $# -lt 2 ]
then
	echo "USAGE: $0 <aligned fasta file> <label appended to fasta name for smc files> [any additional argweaver parameters]"
	exit 1
fi


# formatting within the fastq alignment file (remove gaps and *'s)
name=`basename $1` 
if [[ $name == *.fa ]]
then
	output_label=${name%.fa}_$2
elif [[ $name == *.fasta ]]
then
	output_label=${name%.fasta}_$2
else
	echo -e "\nError: Unknown fasta extenstion: $1.\nAccepted fasta file extensions are .fa and .fasta"
	exit 1
fi

sed -e 's/-/N/g' $1 > $output_label.fasta # replace - with N
sed -i -e 's/*/_/g' $output_label.fasta # replace * with _

# make output directory and cd to that location
mkdir $output_label
cd $output_label

# convert aligned fasta to a sites file (used as argweaver input)
python ../fasta2sites.py ../$output_label.fasta $output_label.sites
echo "finished fasta2sites" > pipeline.log
####

#### RUN ARGWEAVER

arg-sample --sites $output_label.sites \
	   -o $output_label \
	   --overwrite \
	   ${@:3} \
	   | tee $output_label.log
echo -e "arg sample command run:\n
      \targ-sample --sites $output_label.sites \\n
	   -o $output_label \\n
	   --overwrite \\n
	   ${@:3}" >> pipeline.log
echo "finished arg-sample" >> pipeline.log

####


#### COMPUTE PAIRWISE DISTANCES AND ADD CORRECT LABELS

SMC_file=$output_label.1000.smc.gz
#Rscript ../pairwise_parser.R $SMC_file $output_label.min_dist_matrix.tsv
Rscript ../smc_parser.R $SMC_file $output_label.sites
echo "finished smc_parser" >> pipeline.log

####


#### CLEANUP

mkdir smc_files_$output_label
mv *smc.gz smc_files_$output_label

####

echo "FINISHED! Data file is found at treeList.RData" >> pipeline.log
