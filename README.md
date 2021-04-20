# Team 1 Comparative Genomics

#### MEMBERS: 
Anneke Augenbroe, Davina Campbell, Harini Adepu, Jintian Lyu, Marcus Valancius, Nikesh Kumar, Swapna Kannan

This repository includes a script to run a Comparative Genomics pipeline. <Information about what comparative genomics is?>

**```Comparative_master_pipeline.sh```**
This script is to install and run a pipeline which takes assembled or raw fasta sequences and runs them through various tool to gain more insight into the genetic relatedness of the samples. The script can be run with and without tool installations. The script can also run as few as one tool and as many as all of the tools. Further information about the specific outputs are listed below. 
Further details available at our [class wiki page].

https://github.gatech.edu/comgenomics2021/Team1-ComparativeGenomics
  
##### PREREQUISITES:
-    git
-    conda 
-    tools folder with SRST2 and VF database 

##### OPTIONS
        -t	    installs tools 
        -o      Output name for the files
        -i      PATH for input of  assembled .fasta files
        -I      PATH for input of raw reads. This is used for stringMLST
        -g	    PATH for annotated gff files. This is used for resistance
        -r      fasta file for the reference genome
        -s      sample name to root phylogenetic tree; MUST BE DIFFERENT THAN REFERENCE 
        -b      run pyANI with ANIb
        -m      run pyANI with ANIm
        -M      run stringMLST
        -p      run parSNP
        -P      run PlasmidFinder
        -V      run SRST2 and BLAST
        -R      Find resistance genes

##### TOOLS INSTALLED/INVOKED:
        SNP Level: parSNP
        Whole Genome Level:pyANI, stringMLST
        Virulence Level: VFDB,SRST2,BLAST
        Accessory DNA: PlasmidFinder

## ANIb
ANIb is run on Pyani. This is a Python3 module for calculating ANI and  other relatedness measures (alignment coverage, alignment lengths) for whole genome comparisons and creating summary figures. When running ANIb calculations you are using BLAST for robust analysis.

The command used was:
```
Comparative_master_pipeline.sh -t -i <assembled_fasta_input> -m -o <output name>
``` 
## ANIm
ANIm is run on Pyani. This is a Python3 module for calculating ANI and  other relatedness measures (alignment coverage, alignment lengths) for whole genome comparisons and creating summary figures. When running ANIm calcualtions you are using Mummer to calculate for robust analysis.

The command used was:
  
```
Comparative_master_pipeline.sh -t -i <assembled_fasta_input> -m -o <output name>
``` 
  
## stringMLST
stringMLST is a rapid k-mer based 7 gene MLST tool for our analyses. This tool can be easily installed and it has an input of raw sequene reads. This tool uses access to the PubMLST databases as a reference. 
  
The command used was:
```
Comparative_master_pipeline.sh -t -I <raw input fasta file> -M -o <output name> -s <sample to root tree>
```
## parSNP
PARSNP is a core genome SNP analysis tool. It requires reference genome. It identifies maximum unique matches (MUMs) and utilizes compressed suffix graph (CSG) to index the reference genome to identify multi-MUMs. It aligns genomes of any size and it is optimized for microbial genomes. Once the reference is built then other genomes are ran through the CSG. It takes Fasta file as input and output core-genome alignments, variant calls and SNP trees. 

The command used was:
  
```
Comparative_mater_pipeline.sh -t -i <assembled_fasta_input> -p -o <output name> -r <reference file> -s <sample to root tree>
``` 

## Virulence 
<ADD DESCRIPTION>
The files for the Virulence database is in the github folder. If you clone this repository, the tools folder has the required information.
  
The command used was:

```
Comparative_master_pipeline.sh -t -i <assembled_fasta_input> -V -o <output name>
```

## Antimicrobial Resistance Genes 
Taking an annotated GFF file for Antimicrobial Resistance genes, this portion of the script shows you if there is presence of a certain AMR resistance category in a particular isolate. 

```
Comparative_master_pipeline.sh -g <GFF PATH> -R -o <output name>
``` 
## Plasmidfinder 
This tool aims to identify and detect plasmid replicons and assign to Incompatibility groups. It takes sequence data (raw assembled) as input. The default parameters: Detection of replicons with 80%+ nucleotide identity, 60%+ coverage. Its outputs are matched sequences and location on sequence. Its advantage over query to blastn: immediate classification of plasmid to existing plasmid lineages.

```
Comparative_master_pipeline.sh -t -P -i <assembled_fasta_input> -o <output name> 
``` 
