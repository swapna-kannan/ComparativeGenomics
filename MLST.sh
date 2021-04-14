#!/bin/bash

print_help() { echo "
USAGE
	MLST.sh [OPTIONS...] <INPUT_READS_DIRECTORY> 
DESCRIPTION
A script to generate an allelic profile matrix for sequence typing samples from paired-end read fastq files using stringMLST.
The script is optimized to run for sequencing data from Campylobacter jejuni isolates.
Script can be run to install the necessary tools and dependencies if running for the first time.

Required Arguments:
	-r 										path to directiory of paired end reads

Optional Arguments:
	-h 										display help
	-s 	<SPECIES_NAME>						name of species (default is "Campylobacter jejuni")
	-i 										install pipeline
	-o 	<OUTPUT_NAME> 						output file name
"
}

reads=""
species="Campylobacter jejuni"
output="7gMLST.tsv"
install=false

while getopts "is:r:o:h" option
do
	case $option in
		i) install=true;;
		s) species=$OPTARG;; #species arg may be needed here if we chose to allow Typing of other species... for now I'll just leave the option in here and it can be removed later
		r) reads=$OPTARG;;
		o) output=$OPTARG;;
		h) print_help
			exit;;
		*) echo "UNKNOWN OPTION $OPTARG PROVIDED"
			exit;;
	esac
done

#installation instructions if install flag is called:
if $install
then
	echo "Installing stringMLST and its dependencies"
	conda install -c bioconda stringmlst
fi
echo"downloading database for $species from pubMLST..."
stringMLST.py --getMLST -P datasets/ --species "Campylobacter jejuni"

echo "configuring database for $species..."
stringMLST.py --buildDB --config datasets/Campylobacter_jejuni_config.txt -k 35 -P CJ

echo "running sequence typing for paired end reads..."
stringMLST.py --predict -d $reads -p --prefix CJ -k 35 -o $output

echo "Done!"


