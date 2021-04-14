#!/bin/bash

reads=""
species="Campylobacter jejuni"
output="7gMLST.tsv"
install=false

while getopts "is:r:o:" option
do
	case $option in
	i) install=true;;
	s) species=$OPTARG;;
	r) reads=$OPTARG;;
	o) output=$OPTARG;;
	esac
done

#installation instructions if install flag is called:
if $install
then
	conda install -c bioconda stringmlst
fi

stringMLST.py --getMLST -P datasets/ --species "Campylobacter jejuni"
stringMLST.py --buildDB --config datasets/Campylobacter_jejuni_config.txt -k 35 -P CJ
stringMLST.py --predict -d $reads -p --prefix CJ -k 35 -o $output


