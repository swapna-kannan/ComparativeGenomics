sage that the user is going to see
print_help() {echo "
USAGE
        #need to add the final usage of the pipeline here

DESCRIPTION
This is a script to install and run a pipeline for Comparative Genomics Analysis. There are multiple tools that you can choose from and run. The script takes ineither FASTA files or raw reads,from Illumina samples(#should we include this here or not?), for input and you can run it through tools such as ANI, MLST, SNP Analysis, VFDB, PlasmidFinder.

PREREQUISTITES:
        git
        conda
        #ADD MORE

TOOLS INSTALLED/INVOKED:
        SNP Level: kSNP3, parSNP
        Whole Genome Level:pyANI, stringMLST
        Virulence Level: VFDB,SRST2,BLAST
        Accessory DNA: PlasmidFinder
OPTIONS
        -i      PATH for input of  assembled .fasta files.
        -o      Output directory for the files.
        -t      Folder where you want the tools to be downloaded.
        -I      PATH for input of raw reads. This is used for stringMLST.
        -A      run pyANI
        -b      run pyANI with ANIb
        -m      run pyANI with ANIm
        -M      run stringMLST
        -K      run kSNP3
        -p      run parSNP
        -P      run PlasmidFinder
        -R      .fasta file for the refernce genome
"
}
#our GETOPTS BLOCK
outputDir = ""
toolsDir = ""
pyANI = false
ANIb = false
ANIm = false
stringMLST = false
kSNP3 = false
PlasmidFinder = false
parSNP = false
ref_genome = false
while getopts "hiIo:t:AbmMKpPR" option
do
        case $option in
                h)print_help exit;;
                i)assembled_input=$OPTARG;;
                I)raw_input=$OPTARG;;
                o)outputDir=$($OPTARG);;
                t)toolsDir=$($OPTARG);;
                A)pyANI=true;;
                b)ANIb=true;;
                m)ANIm=true;;
                M)stringMLST=true;;
                K)kSNP=true;;
                p)parSNP=true;;
                P)PlasmidFinder=true;;
                R)ref_genome=true;;
 		*) echo "UNKNOWN OPTION $OPTARG PROVIDED" exit;;
        esac
done

	
