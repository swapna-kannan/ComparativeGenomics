#!/bin/bash
#Usage that the user is going to see
print_help() { echo "
USAGE
        #need to add the final usage of the pipeline here

DESCRIPTION
This is a script to install and run a pipeline for Comparative Genomics Analysis.
There are multiple tools that you can choose from and run. The script takes in either FASTA files or raw reads,
for input and you can run it through tools such as ANI, MLST, SNP Analysis, VFDB, PlasmidFinder.

PREREQUISTITES:
        git
        conda
        tools folder with SRST2
        #ADD MORE

TOOLS INSTALLED/INVOKED:
        SNP Level: parSNP
        Whole Genome Level:pyANI, stringMLST
        Virulence Level: VFDB,SRST2,BLAST
        Accessory DNA: PlasmidFinder
OPTIONS
        -i      PATH for input of  assembled .fasta files.
        -o      Output name for the files.
        -I      PATH for input of raw reads. This is used for stringMLST.
	-g	PATH for annotated gff files. This is used for resistance. 
	-N 	To install any of the tools 
        -r      fasta file for the reference genome
        -b      run pyANI with ANIb
        -m      run pyANI with ANIm
        -M      run stringMLST
        -p      run parSNP
        -P      run PlasmidFinder
        -V      run SRST2 and BLAST
        -R      Find resistance genes
"
}

#our GETOPTS BLOCK
assembled_input=false
raw_input=false 
ref_genome=false
gff_files=false
ANIb=false
ANIm=false
stringMLST=false
PlasmidFinder=false
parSNP=false
virulence=false
resistance=false
tools=false

while getopts "hi:I:r:g:o:bmMpPVRt" option; do
	case "$option" in
		h)print_help exit;;
		i)assembled_input=$OPTARG;;
		I)raw_input=$OPTARG;;
		r)ref_genome=$OPTARG;;
		g)gff_files=$OPTARG;; 
		o)output=$OPTARG;;
		b)ANIb=true;;
		m)ANIm=true;;
		M)stringMLST=true;;
		p)parSNP=true;;
		P)PlasmidFinder=true;;
		V)virulence=true;;
		R)resistance=true;;
		t)tools=true;;
		*)echo "UNKNOWN OPTION $OPTARG PROVIDED" exit;; #isn't echoing the option provided 
	esac
done

# make general CompGen directory with tools and outputs 
mkdir -p CompGen/tools CompGen/output

# running aniB
if $ANIb; then
	
	#check if input files exist
	if [[ -f $assembled_input ]]
	then 
		echo "Assemblies do not exist"
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
		#conda install -c bioconda mummer blast legacy-blast -y
		conda install -y -c bioconda/label/cf201901 blast-legacy
		conda install -y -c biocore blast-legacy
	fi 
	

	#run the ANIm command
	echo "Calculating average nucleotide identity using ANIm..."
	average_nucleotide_identity.py -o CompGen/output/ANIb/$output -i $assembled_input -m ANIb -g -f -v
	
	# deactivate pyani_env 
	conda deactivate
	conda env remove -n pyani_env

fi

#running aniM
if $ANIm; then

	#check if input files exist
	if [[ -f $assembled_input ]]
	then 
		echo "Assemblies do not exist"
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
		#conda install -c bioconda mummer blast legacy-blast -y
		conda install -y -c bioconda/label/cf201901 blast-legacy
		conda install -y -c biocore blast-legacy
	fi 
	

	#run the ANIm command
	echo "Calculating average nucleotide identity using ANIm..."
	average_nucleotide_identity.py -o CompGen/output/ANIm/$output -i $assembled_input -m ANIm -g -f -v
	
	# deactivate pyani_env 
	conda deactivate
	conda env remove -n pyani_env

fi

#running MLST
if $stringMLST; then

	#check if input files exist
	if [[ -f $raw_reads ]]
	then 
		echo "Raw reads do not exist"
		exit
	fi

	mkdir -p CompGen/tools/stringMLST CompGen/tools/stringMLST/extra CompGen/output/stringMLST
	cd CompGen/tools/stringMLST

	# Install stringMLST
	if $install; then 
		echo "Installing stringMLST and its dependencies"
		conda install -c bioconda stringmlst

		#Install GrapeTree
		echo "Installing GrapeTree and its dependencies"
		pip install grapetree
	fi 

	# run stringMLST
	echo"downloading database for $species from pubMLST..."
	stringMLST.py --getMLST -P datasets/ --species "Campylobacter jejuni"

	echo "configuring database for $species..."
	stringMLST.py --buildDB --config datasets/Campylobacter_jejuni_config.txt -k 35 -P CJ

	echo "running sequence typing for paired end reads..."
	stringMLST.py --predict -d $raw_input -p --prefix CJ -k 35 -o $outputDir #<-- needs to be updated
	mv ./CompGen/tools/stringMLST/CJ* ./CompGen/tools/stringMLST/extra/
	echo "Done!"

	# run GrapeTree to generate newick file for cluster visualization
	echo "generating Newick file from allele profile"
	#grapetree -p <mlst_output_name> -m "MSTreeV2" > <mlst_output_name.newick> #<-- this needs to be updated... I'll talk about this more with yall tomorrow regarding output names and pathing
fi

#running parsnp
if $parSNP; then

	#check if reference exists
	if [[ -f $ref_genome ]]
	then 
		echo "reference genome do not exist"
		exit
	fi
     
     	# check if input is correct 
     	if [[ -f $assembled_input ]]
	then 
		echo "Assemblies do not exist"
		exit
	fi


	#make tools and output directory for parsnp
	mkdir -p CompGen/tools/parsnp CompGen/tools/parsnp/extra CompGen/output/parsnp

	#download tools
	if $install; then 
		echo "Installing parSNP..."
		conda install -y -c bioconda parsnp
	fi 

	#run parsnp
	parsnp -r $assembled_input/$ref_genome -d $assembled_input -o CompGen/output/parsnp/$output

	#FIND OUTPUT FILE NAME AND MOVE TO OUTPUT FOLDER 
fi



# running virulence
if $virulence; then
	
	# check if input is correct 
	if [[ -f $assembled_input ]]
	then
		echo "Assemblies do not exist"
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
	if [[ -f $assembled_input ]]
	then
		echo "Assemblies do not exist"
		exit
	fi
	
	#making the directories 
	mkdir -p CompGen/tools CompGen/tools/PlasmidFinder/extra CompGen/output/PlasmidFinder
	cd CompGen/PlasmidFinder/tools
	
	#installing PlasmidFinder 
	if $install; then 
		echo "Installing Plasmidfinder"
		conda install -c -y bioconda plasmidfinder 
		download-db.sh
	fi  
	
	echo "Running plasmidfinder"  
    	for file in $(ls $assembled_input);  
    	do 
    		echo "$file"
    		v=$(echo $file | cut -d "." -f 1)
    		v1="${output}${v}"
    		plasmidfinder.py -i $assembled_input/$file > CompGen/output/PlasmidFinder/$v1 
    	done

	
	#run PlasmidFinder
	#echo "Running Plasmidfinder"
	#plasmidfinder.py -i $assembled_input > CompGen/output/PlasmidFinder/log.txt
	# -o CompGen/output/PlasmidFinder/$output > 
fi 
