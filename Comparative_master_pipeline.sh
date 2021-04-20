#!/bin/bash
#Usage that the user is going to see
print_help() { echo "
USAGE
        Comparative_master_pipeline [OPTIONS...] < -o output name > [-t] [-i ASSEMBLED_INPUT_READS_DIRECTORY] [-I RAW_INPUT_READS_DIRECTORY] [-g GFF_FILES_DIRECTORY] [-r PARSNP_REFERENCE_FILE] [-b] [-m] [-M] [-p] [-P] [-V] [-R]


DESCRIPTION
This is a script to install and run a pipeline for Comparative Genomics Analysis.
There are multiple tools that you can choose from and run.
The script takes in FASTA files for ANI, SNP analysis, Virulence and PlasmidFinder. It takes in raw reads for MLST and it takes in gff annonated files to find resistance genes.

PREREQUISTITES:
        git
        conda
        tools folder with SRST2 and VF database

TOOLS INSTALLED/INVOKED:
        SNP Level: parSNP
        Whole Genome Level:pyANI, stringMLST
        Virulence Level: VFDB,SRST2,BLAST
        Accessory DNA: PlasmidFinder
OPTIONS
        -t	installs tools 
        -o      Output name for the files
        -i      PATH for input of  assembled .fasta files
        -I      PATH for input of raw reads. This is used for stringMLST
        -g	PATH for annotated gff files. This is used for resistance
        -r      fasta file for the reference genome
        -b      run pyANI with ANIb
        -m      run pyANI with ANIm
        -M      run stringMLST
        -p      run parSNP NOTE: when running parSNP, pass the -r command and select a sample for your reference, then pass the -s command and select a DIFFERENT sample for rooting the tree
        -P      run PlasmidFinder
        -V      run SRST2 and BLAST
        -R      Find resistance genes
"
}

#our GETOPTS BLOCK
assembled_input=$false
raw_input=$false
ref_genome=$false
gff_files=false
ANIb=false
ANIm=false
stringMLST=false
PlasmidFinder=false
parSNP=false
virulence=false
resistance=false
tools=false
sample=$false

while getopts "hi:I:r:g:o:s:bmMpPVRt" option
do
	case $option in
		h) print_help
			exit;;
		i) 	assembled_input=$OPTARG;;
		I) 	raw_input=$OPTARG;;
		r) 	ref_genome=$OPTARG;;
		g) 	gff_files=$OPTARG;;
		o) 	output=$OPTARG;;
		s) 	sample=$OPTARG;;
		b) 	ANIb=true;;
		m) 	ANIm=true;;
		M) 	stringMLST=true;;
		p) 	parSNP=true;;
		P) 	PlasmidFinder=true;;
		V) 	virulence=true;;
		R) 	resistance=true;;
		t) 	tools=true;;
		*) 	echo "UNKNOWN OPTION PROVIDED"
			exit;; #isnt echoing the option provided
	esac
done

# make general CompGen directory with tools and outputs 
mkdir -p CompGen/tools CompGen/output

# running aniB
if $ANIb; then
	
	#check if input files exist
	if [ $ANIb -a -z $assembled_input ]
	then 
		echo "Assemblies do not exist. Please call the -i flag and provide a path to the input directory of assembled reads"
		exit
	fi

	#make tools and output directory for ANIm
	mkdir -p CompGen/tools/ANIb CompGen/tools/ANIb/extra CompGen/output/ANIb
	
	 conda create --name pyani_env python=3.8 -y
	 eval "$(conda shell.bash hook)"
	 source ./anaconda3/bin/activate pyani_env
	 conda activate pyani_env
		
	#download tools
	if $tools; then
		echo "Installing pyani..."
		conda install -y biopython
		conda install -y -c bioconda pyani
		conda install -y -c bioconda blast-legacy
		conda install -y -c bioconda/label/cf201901 blast-legacy
		conda install -y -c biocore blast-legacy
	fi
	
	#run the ANIm command
	echo "Calculating average nucleotide identity using ANIm..."
	average_nucleotide_identity.py -o CompGen/output/ANIb/$output -i $assembled_input -m ANIb -g -f -v
	
	 deactivate pyani_env
	 conda deactivate
	 conda env remove -n pyani_env

fi

#running aniM
if $ANIm; then
	echo "starting ANIm processes"
	#check if input files exist
	if [ $ANIm -a -z $assembled_input ];
	then
		echo "Assemblies do not exist. Please call the -i flag and provide a path to the input directory of assembled reads"
		exit
	fi

	#make tools and output directory for ANIm
	mkdir -p CompGen/tools/ANIm CompGen/tools/ANIm/extra CompGen/output/ANIm
	
	 conda create --name pyani_env python=3.8 -y
	 eval "$(conda shell.bash hook)"
	 source ./anaconda3/bin/activate pyani_env
	 conda activate pyani_env

	#download tools
	if $tools; then
		echo "Installing pyani..."
		conda install -y biopython
		conda install -y -c bioconda pyani
		conda install -y -c bioconda blast-legacy
		conda install -y -c bioconda mummer blast legacy-blast
		conda install -y -c bioconda/label/cf201901 blast-legacy
		conda install -y -c biocore blast-legacy
	fi 
	

	#run the ANIm command
	echo "Calculating average nucleotide identity using ANIm..."
	average_nucleotide_identity.py -o CompGen/output/ANIm/$output -i $assembled_input -m ANIm -g -f -v --maxmatch
	
	 conda deactivate pyani_env
	 conda deactivate
	 conda env remove -n pyani_env

fi

#running MLST
if $stringMLST; then

	# check if input files exist
	if [ $stringMLST -a -z $raw_input ];
	then 
		echo "Raw reads input for provided. Please call the -I flag and provide path to a directory of raw reads."
		exit
	fi

	# check if Sample name provided alongside stringMLST flag
	if [ $stringMLST -a -z $sample ]; then
		echo "No sample name provided. Please call the -s flag and enter a sample name"
		exit
	fi

	mkdir -p CompGen/tools/stringMLST CompGen/tools/stringMLST/extra CompGen/output/stringMLST

	# Install stringMLST
	if $tools; then
		echo "Installing stringMLST and its dependencies"
		conda install -c bioconda stringmlst -y

		# Install GrapeTree
		echo "Installing GrapeTree and its dependencies"
		pip install grapetree 
		conda install pandas -y
		conda install requests -y

		# Install Toytree
		echo "Installing Toytree and its dependencies"
		conda install toytree -c conda-forge -y
	fi
	
	# run stringMLST
	echo "downloading database for $species from pubMLST..."
	stringMLST.py --getMLST -P $PWD/CompGen/tools/stringMLST/datasets/ --species "Campylobacter jejuni"
	
	echo "configuring database for $species..."
	stringMLST.py --buildDB --config $PWD/CompGen/tools/stringMLST/datasets/Campylobacter_jejuni_config.txt -k 35 -P CJ
	
	echo "running sequence typing for paired end reads..."
	stringMLST.py --predict -d $raw_input -p --prefix CJ -k 35 -o $PWD/CompGen/output/stringMLST/7gMLST_${output}.tsv
	
	mv $PWD/CJ* $PWD/CompGen/tools/stringMLST/extra/
	
	echo "Done with stringMLST!"

	# run GrapeTree to generate newick file for cluster visualization
	echo "generating Newick file from allele profile"
	grapetree -p $PWD/CompGen/output/stringMLST/7gMLST_${output}.tsv -m "MSTreeV2" > $PWD/CompGen/output/stringMLST/7gMLST_${output}.newick
	
	echo "Newick file generated"
	echo "Generating tree PDF"

python3 - << EOF
import toytree
import toyplot
import toyplot.pdf
import numpy as np
import subprocess as sp

with open('$PWD/CompGen/output/stringMLST/7gMLST_${output}.newick', 'r') as fh:
	newick = fh.read()
print(newick)
tre = toytree.tree(newick)

rtre = tre.root(wildcard="${sample}")
#rtre.draw(tip_labels_align=True);

canvas, axes, mark = rtre.draw()

toyplot.pdf.render(canvas, "$PWD/CompGen/output/stringMLST/MLSTtree_${output}.pdf")
EOF

	echo "MLST profile tree generated"

fi

#running parsnp
if $parSNP; then

	#check if reference exists
	if [ $parSNP -a -z $ref_genome ];
	then 
		echo "Reference genome not provided. Please call the -r flag and specify an assembled genome to serve as reference."
		exit
	fi

	if [ $parSNP -a -z $sample ]; then
		echo "No sample name provided. Please call the -s flag and enter a sample name"
		exit
	fi

	if [ $sample -eq $ref_gemome ]; then
		echo "For the -s flag, please provide a sample that is different from the reference."
		exit
	fi
	# check if input is correct 
	if [ $parSNP -a -z $assembled_input ]
	then 
		echo "Assembled genomes not provided. Please call the -i flag and specify the file path to assembled genomes directory"
		exit
	fi


	#make tools and output directory for parsnp
	mkdir -p CompGen/tools/parsnp CompGen/tools/parsnp/extra CompGen/output/parsnp

	#download tools
	if $tools; then
		echo "Installing parSNP..."
		conda install -y -c bioconda parsnp

	# Install Toytree
		echo "Installing Toytree..."
		conda install toytree -c conda-forge -y
	fi 

	#run parsnp
	parsnp -r $assembled_input/$ref_genome -d $assembled_input -o CompGen/output/parsnp/$output

	# generate tree with Toytree
	python3 - << EOF
import toytree
import toyplot
import toyplot.pdf
import numpy as np
import subprocess as sp

with open('$PWD/CompGen/output/parsnp/$output/parsnp.tree', 'r') as fh:
	newick = fh.read()
tre = toytree.tree(newick)

rtre = tre.root(wildcard="${sample}")
#rtre.draw(tip_labels_align=True);

canvas, axes, mark = rtre.draw()

toyplot.pdf.render(canvas, "$PWD/CompGen/output/parsnp/${output}/${output}_tree.pdf")
EOF

	echo "SNP tree generated"

fi

# running virulence
if $virulence; then
	
	# check if input is correct
	if  [$virulence -a -z $assembled_input ]
	then
		echo "Assemblies do not exist. Please call the -i flag and provide a path to the input directory of assembled reads"
		exit
	fi

	#make tools and output directory for virulence
	mkdir -p CompGen/tools/virulence CompGen/tools/virulence/extra CompGen/output/virulence

	# create a virtual env with python 2 to run srst2
	conda create -y --name temp_py2 python=2.7
	eval "$(conda shell.bash hook)"
	source ./anaconda3/bin/activate temp_py2 #need to check path to conda
	conda activate temp_py2
	
	#virulence needs the tools installed everytime bc it's making a new environment
	echo "Installing biopython"
	conda install -y biopython
	echo "Installing srst2"
	conda install -y -c bioconda srst2
	echo "Installing cd-hit"
	conda install -y -c bioconda cd-hit
	echo "Installing blast"
	conda install -y -c bioconda blast

	# generate the reference files for genus
	python ./tools/srst2/database_clustering/VFDBgenus.py --infile ./tools/VFs.ffn --genus Campylobacter #need to have tool folder and reference
	cd-hit -i Campylobacter.fsa -o Campylobacter_cdhit90 -c 0.90 > Campylobacter_cdhit90.stdout 
	python ./tools/srst2/database_clustering/VFDB_cdhit_to_csv.py --cluster_file Campylobacter_cdhit90.clstr --infile Campylobacter.fsa --outfile Campylobacter_cdhit90.csv
	python ./tools/srst2/database_clustering/csv_to_gene_db.py -t Campylobacter_cdhit90.csv -o Campylobacter_VF_clustered.fasta -s 5

	# deactivate and delete the virtual env
	conda deactivate
	#conda env remove -n temp_py2

	# make the blast database
	makeblastdb -in Campylobacter_VF_clustered.fasta -dbtype 'nucl' -out CompGen/tools/virulence/extra/Campylobacter_database

	# run blast on all the files and put results in output folder
	for file in $(ls $assembled_input); do
		blastn -db CompGen/tools/virulence/extra/Campylobacter_database -query $assembled_input/$file -perc_identity .98 -out CompGen/tools/virulence/$file -outfmt "6 stitle"
	done


	# get all the gene names
	for file in $(ls CompGen/tools/virulence/*fasta); do
		awk '{print $3}' $file >> CompGen/tools/virulence/extra/VF_all_$output.txt
	done

	# get unique genes and add them as columns in a new file
	sort CompGen/tools/virulence/extra/VF_all_$output.txt | uniq >> CompGen/tools/virulence/extra/VF_unique_$output.txt
	long_line=""
	for line in $(cat CompGen/tools/virulence/extra/VF_unique_$output.txt); do
        	long_line+="	$line"
	done

	echo "$long_line" > CompGen/tools/virulence/VF_table_$output.txt

	# check if gene is in each file and add info to VF_table.txt
	for file in $(ls CompGen/tools/virulence/*fasta); do
        	data="$file     "
        	for line in $(cat CompGen/tools/virulence/extra/VF_unique_$output.txt); do
                	if grep -q $line $file; then
                        	data+="X	"
               		else
                        	data+="	"
                	fi
        	done
        	echo "$data" >> CompGen/output/virulence/VF_table_$output.txt
	done

	mv Campylobacter* CompGen/tools/virulence/extra
fi


#running resistance
if $resistance; then
	
	#make the output directory 
	mkdir -p CompGen/output/Deeparg CompGen/tools/Deeparg

	# get all the gene names
	for file in $(ls $gff_files); do
		echo $file
        	grep DeepARG $gff_files/$file  | awk '{print $9}' | sed 's/|/,/g' >> CompGen/tools/Deeparg/res_temp_$output.txt
	done

	for line in $(cat CompGen/tools/Deeparg/res_temp_$output.txt); do
        	readarray -d , -t strarr <<< "$line"
        	echo ${strarr[2]} >> CompGen/tools/Deeparg/res_all_$output.txt
	done 

	# get unique genes and add them as columns in a new file
	sort CompGen/tools/Deeparg/res_all_$output.txt | uniq >> CompGen/tools/Deeparg/res_unique_$output.txt
	long_line=""
	for line in $(cat CompGen/tools/Deeparg/res_unique_$output.txt); do
        	long_line+="	$line"
	done

	echo "$long_line" > CompGen/output/Deeparg/res_table_$output.txt

	# check if gene is in each file and add info to VF_table.txt
	for file in $(ls $gff_files); do
        	data="$file	"
        	for line in $(cat CompGen/tools/Deeparg/res_unique_$output.txt); do
                	if grep -q $line $gff_files/$file; then
                        	data+="X	"
                	else
                        	data+="	"
                	fi
        	done
        echo "$data" >> CompGen/output/Deeparg/res_table_$output.txt
	done

fi 

#running PlasmidFinder
if $PlasmidFinder; then
	
	# check if input is correct 
	if [ $PlasmidFinder -a -z $assembled_input ]
	then
		echo "Assemblies do not exist. Please call the -i flag and provide a path to the input directory of assembled reads"
		exit
	fi
	
	#making the directories
	mkdir -p CompGen/tools CompGen/tools/PlasmidFinder/extra CompGen/output/PlasmidFinder
	cd CompGen/tools/PlasmidFinder
	
	#installing PlasmidFinder
	if $tools; then
		echo "Installing Plasmidfinder"
		conda install -c bioconda plasmidfinder -y
		download-db.sh
	fi
	
	echo "Running plasmidfinder"  
    	for file in $(ls $assembled_input); do
    		echo "$file"
    		v=$(echo $file | cut -d "." -f 1)
    		v1="${output}${v}"
    		plasmidfinder.py -i $assembled_input/$file > CompGen/output/PlasmidFinder/$v1 
    	done

	mv data.json CompGen/output/PlasmidFinder/${output}_data.json
	mv tmp $PWD/CompGen/tools/PlasmidFinder
fi