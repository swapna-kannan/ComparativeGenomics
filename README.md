# Team 1 Comparative Genomics

#### MEMBERS: 
Anneke Augenbroe, Davina Campbell, Harini Adepu, Jintian Lyu, Marcus Valancius, Nikesh Kumar, Swapna Kannan

This repository includes a script to run a Comparative Genomics pipeline. <Information about what comparative genomics is?>

**```comparative_master_pipeline.sh```**
This script is to install and run a pipeline which takes assembled or raw fasta sequences and runs them through various tool to gain more insight into the genetic relatedness of the samples. The script can be run with and without tool installations. The script can also run as few as one tool and as many as all of the tools. Further information about the specific outputs are listed below. 
Further details available at our [class wiki page].
<INSERT CLASS WIKI PAGE>
  
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
<ADD DESCRIPTION>


The command used was:
```
comparative_master_pipeline.sh -t -i <assembled_fasta_input> -m -o <output name>
``` 
## ANIm
<ADD DESCRIPTION>
The command used was:
  
```
comparative_master_pipeline.sh -t -i <assembled_fasta_input> -m -o <output name>
``` 
  
## stringMLST
<ADDDESCRIPTION>
  
The command used was:
```
comparative_master_pipeline.sh -t -I <raw input fasta file> -M -o <output name> -s <sample to root tree>
```
## parSNP
<ADD DESCRIPTION>
The command used was:
  
```
comparative_mater_pipeline.sh -t -i <assembled_fasta_input> -p -o <output name> -r <reference file> -s <sample to root tree>
``` 

## Virulence 
<ADD DESCRIPTION>
The command used was:

```
comparative_master_pipeline.sh -t -i <assembled_fasta_input> -V -o <output name>
```

## Antimicrobial Resistance Genes 
<ADD DESCRIPTION>

```
comparative_master_pipeline.sh -g <GFF PATH> -R -o <output name>
``` 
## Plasmidfinder 
<ADD DESCRIPTION> 

```
comparative_master_pipeline.sh -t -P -i <assembled_fasta_input> -o <output name> 
``` 
