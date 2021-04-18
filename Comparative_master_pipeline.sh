#!/bin/bash
#Usage that the user is going to see
print_help() { echo "
USAGE
        #need to add the final usage of the pipeline here

DESCRIPTION
This is a script to install and run a pipeline for Comparative Genomics Analysis.
There are multiple tools that you can choose from and run. The script takes ineither FASTA files or raw reads,from Illumina samples
(#should we include this here or not?
I think we should remove this. I think reads from most platforms will work -Nikesh ),
for input and you can run it through tools such as ANI, MLST, SNP Analysis, VFDB, PlasmidFinder.

PREREQUISTITES:
        git
        conda
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
outputDir=""
toolsDir=""
ref_genome=false
ANIb=false
ANIm=false
stringMLST=false
kSNP3=false
PlasmidFinder=false
parSNP=false
virulence=false
resistance=false
install=false
while getopts "hiIro:bmMpPVRN" option
do
        case $option in
                h)print_help exit;;
                i)assembled_input=$OPTARG;;
                I)raw_input=$OPTARG;;
                r)ref_genome=$OPTARG;;
                o)output=$($OPTARG);;
                b)ANIb=true;;
                m)ANIm=true;;
                M)stringMLST=true;;
                p)parSNP=true;;
                P)PlasmidFinder=true;;
                V)virulence=true;;
                R)resistance=true;;
		N)install=true;;
                *) echo "UNKNOWN OPTION $OPTARG PROVIDED" exit;;
        esac
done

mkdir -p CompGen/tools Compgen/output

#NEED ADD CONDA ENVIRONMENT


# running aniB
if $ANIb; then

	#ADD CHECK FOR ASSEMBLED READ INPUT 

	#make tools and output directory for ANIb
	mkdir -p CompGen/tools/ANIb CompGen/tools/ANIb/extra CompGen/output/ANIb
	cd CompGen/tools/ANIb

	if $install; then
		echo "Installing  pyANI..."
		#downloading the tools 
		conda install pyani 
		conda install mummer blast legacy-blast -y
	fi   

	# run the ANIb command 
	echo "Calculating average nucleotide identity using ANIb..."
	average_nucleotide_identity.py -i $assembled_input -o CompGen/output/aniB/${output} -m ANIb -g -f -v 
	#move files into specific folder 
	mv ${output}ANIb_alignment_coverage.* CompGen/output/ANIb
	mv ${output}ANIb_alignment_length.* CompGen/output/ANIb
	mv ${output}ANIb_hadamard.* CompGen/output/ANIb
	mv ${output}ANIb_percentage_identity.* CompGen/output/ANIb
	mv ${output}ANIb_similarity_errors.* CompGen/output/ANIb
fi

#running aniM
if $ANIm; then

	#ADD CHECK FOR ASSEMBLED READ INPUT 
	
	#make tools and output directory for ANIm
        mkdir -p CompGen/tools/ANIm CompGen/tools/ANIm/extra Compgen/output/ANIm
	cd CompGen/tools/ANIm

        #move to ANIm tools folder so any junk files get outputted there 
        cd CompGen/tools/ANIm
	
	#download tools
	if $install; then
		echo "Installing pyani..." 
        	conda install pyani
        	conda install mummer blast legacy-blast -y

	#run the ANIm command
	echo "Calculating average nucleotide identity using ANIm..."
	average_nucleotide_identity.py -i $assembled_input -o ../../output/ANIm/$output -m ANIm -g -f -v
	
	#moving the files into the appropriate directories
	#I JUST REALIZED WE DONT NEED TO DO THIS...FAIL 
	mv ${output}ANIm_alignment_coverage.* CompGen/output/ANIm
	mv ${output}ANIm_alignmnet_length.* CompGen/output/ANIm
	mv ${output}ANIm_hadamard.* CompGen/output/ANIm
	mv ${output}ANIm_percentage_identity.* CompGen/output/ANIm
	mv ${output}ANIm_similarity_errors.* CompGen/output/ANIm

fi

#running MLST
if $stringMLST; then

	#ADD CHECK FOR RAW READ INPUT 

	mkdir -p CompGen/tools/stringMLST Compgen/tools/stringMLST/extra CompGen/output/stringMLST
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

	#make tools and output directory for parsnp
        mkdir -p CompGen/tools/parsnp CompGen/tools/parsnp/extra Compgen/output/parsnp

        #move to parsnp tools folder so any junk files get outputted there 
        cd CompGen/tools/parsnp

	#download tools
	if $install; then 
		echo "Installing parSNP..."
		conda install -y -c bioconda parsnp
	fi 

	#check if reference exists
	if [ ! -f $ref_genome ]; then
        	echo "Reference does not exist"
        	exit
	fi

	# check if input files exist
	if [ ! -f $assembled_input ]; then
        	echo "Assemblies do not exist"
        exit
	fi

	#run parsnp WHY IS THIS IN $?
	$(parsnp -r $ref_genome -d $assembled_input -o CompGen/parsnp/$output)

	#FIND OUTPUT FILE NAME AND MOVE TO OUTPUT FOLDER 
fi



# running virulence
if $virulence; then
	#make tools and output directory for virulence 
	mkdir -p CompGen/tools/virulence CompGen/tools/virulenece/extra CompGen/output/virulence
	cd CompGen/tools/virulence 
	# create a virtual env with python 2 to run srst2
	conda create -y --name temp_py2 python=2.7
	eval "$(conda shell.bash hook)"
	source ./anaconda3/bin/activate temp_py2 #need to check path to conda 
	conda activate temp_py2
	if $install; then 
		echo "Installing virulence #WHAT IS THE TOOL SPECIFICALLY CALLED"
		#download conda tools
		conda install -y biopython
		conda install -y -c bioconda srst2
		conda install -y -c bioconda cd-hit
	fi
	# generate the reference files for genus
	python ./tools/srst2/database_clustering/VFDBgenus.py --infile ./tools/VFs.ffn --genus Campylobacter #need to have tool folder and reference
	cd-hit -i Campylobacter.fsa -o Campylobacter_cdhit90 -c 0.90 > Campylobacter_cdhit90.stdout 
	python ./tools/srst2/database_clustering/VFDB_cdhit_to_csv.py --cluster_file Campylobacter_cdhit90.clstr --infile Campylobacter.fsa --outfile Campylobacter_cdhit90.csv
	python ./tools/srst2/database_clustering/csv_to_gene_db.py -t Campylobacter_cdhit90.csv -o Campylobacter_VF_clustered.fasta -s 5

	# deactivate and delete the virtual env
	conda deactivate
	conda env remove -n temp_py2

	# make the blast database
	makeblastdb -in Campylobacter_VF_clustered.fasta -dbtype 'nucl' -out Campylobacter_database

	# run blast on all the files and put results in output folder
	for file in $(ls $assembled_input); do
        	blastn -db Campylobacter_database -query $assembled_input/$file -perc_identity .98 -out $outputDir/$file -outfmt "6 stitle"
	done

	#remove any old files if they exist
	rm -f VF_all.txt
	rm -f VF_unique.txt
	rm -f VF_table.txt

	# get all the gene names
	for file in $(ls $outputDir/*.fasta); do
        	awk '{print $3}' $file >> VF_all.txt
	done

	# get unique genes and add them as columns in a new file
	sort VF_all.txt | uniq >> VF_unique.txt
	long_line=""
	for line in $(cat VF_unique.txt); do
        	long_line+="    $line"
	done

	echo "$long_line" > VF_table.txt

	# check if gene is in each file and add info to VF_table.txt
	for file in $(ls $outputDir); do
        	data="$file     "
        	for line in $(cat VF_unique.txt); do
                	if grep -q $line $assembled_input/$file; then
                        	data+="X        "
               		else
                        	data+=" "
                	fi
        	done
        	echo "$data" >> VF_table.txt
	done

	mv VF_table.txt $outputDir

	rm Campylobacter*
	rm VF_all.txt
	rm VF_unique.txt
fi


#running resistance
if $resistance; then
	#making the directories
	#DO WE NEED TO DO THIS?
	mkdir -p CompGen/output/Deeparg
	# remove any old files if they exist
	rm -f res_all.txt
	rm -f res_unique.txt
	rm -f res_table.txt

	# get all the gene names
	for file in $(ls ./merged_gff/*gff); do=
        	grep DeepARG $file  | awk '{print $9}' | sed 's/|/,/g' >> res_temp.txt
	done

	for line in $(cat res_temp.txt); do
        	readarray -d , -t strarr <<< "$line"
        	echo ${strarr[2]} >> res_all.txt
	done 

	# get unique genes and add them as columns in a new file
	sort res_all.txt | uniq >> res_unique.txt
	long_line=""
	for line in $(cat res_unique.txt); do
        	long_line+="    $line"
	done

	echo "$long_line" > res_table.txt

	# check if gene is in each file and add info to VF_table.txt
	for file in $(ls ./merged_gff/*gff); do
        	data="$file     "
        	for line in $(cat res_unique.txt); do
                	if grep -q $line $file; then
                        	data+="X        "
                	else
                        	data+=" "
                	fi
        	done
        echo "$data" >> res_table.txt
	done

	rm res_all.txt
	rm res_unique.txt

	mv res_table.txt $outputDir
fi 

#running PlasmidFinder

if $PlasmidFinder
then
	#making the directories 
	mkdir -p CompGen/tools CompGen/tools/PlasmidFinder/extra CompGen/output/PlasmidFinder
	cd CompGen/PlasmidFinder/tools
	#installing PlasmidFinder 
	if $install
	then 
		echo "Installing Plasmidfinder"
		conda install -c bioconda plasmidfinder 
	fi  
	echo "Running Plasmidfinder"
	plasmidfinder.py -i $assembled_input > CompGen/output/PlasmidFinder
fi
