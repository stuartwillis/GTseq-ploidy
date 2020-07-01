#!/bin/bash

#this script will read all "genos" files in the working directory, collect read numbers from each locus,
#estimate ploidy from default or specified ranges, genotype those at estimated ploidy,
#and create CSV files of genotypes and other information, with one optionally filtered
#to import into Progeny "Lxxxx_ProgGenos90.csv"

set -o pipefail


#Define arguments
while getopts ":a:b:c:d:g:r:l:p:" opt; do
  case $opt in
    a) ploidy_min="$OPTARG"
    ;;
    b) ploidy_max="$OPTARG"
    ;;
    c) llr="$OPTARG"
    ;;
    d) depth="$OPTARG"
    ;;
    g) geno="$OPTARG"
    ;; 
    r) cov="$OPTARG"
    ;;  
    l) library2="$OPTARG"
    ;;  
    p) parallel="$OPTARG"
    ;;  
    \?) 	printf "\nGTseq_GenoCompile_ploidy_v1.sh : Estimates ploidy and genotypes using the read counts in genos files \n \n"
	printf " Requires that the script is run with the directory containing the genos files as the working directory ('pwd')\n "
	printf "Usage:   GTseq_GenoCompile_ploidy_v1.sh [optional arguments] \n"
	printf "\n Argument                      Description   	                              Default\n \n"
	printf " - a   [ploidy min]              : Minimum ploidy to consider                     [4] \n"
	printf " - b   [ploidy max]              : Maximum ploidy to consider                     [6] \n" 
	printf " - c   [confidence]              : Minimum likelihood ratio to accept ploidy      [3] \n"
	printf " - d   [depth/locus]             : Minimum number of reads to genotype a locus    [20] \n" 
	printf " - g   [genotyping rate]         : Minimum genotyping rate to import to Progeny   [90] \n"
	printf " - l   [library name]            : Library name (from wd: working directory)      [wd] \n"
	printf " - p   [parallel]                : Use multiple cores (requires 'parallel') (y/n) [n] \n\n"
	printf " - r   [reads on target]         : Minimum total reads on target to genotype      [0] \n\n"
	exit    

    ;;
  esac
done

##set variable defaults if unspecified by command line arguments
if [ -z "$ploidy_min" ] ; then
	ploidy_min=4
fi

if [ -z "$ploidy_max" ] ; then
	ploidy_max=6
fi

if [ -z "$llr" ] ; then
	llr=3
fi

if [ -z "$depth" ] ; then
	depth=20
fi

if [ -z "$geno" ] ; then
	geno=90
fi

if [ -z "$cov" ] ; then
	cov=0
fi

if [ -z "$parallel" ] ; then
	parallel="n"
fi


##get workding directory
working_dir=$(pwd)
#echo $working_dir

##try to get library name from working directory if not already set
library=$(echo $working_dir | sed 's/.*\(L[0-9]*\)_Genos.*/\1/')

if [ -z "$library2" ] ; then

	library2=$(echo $library | grep "L[0-9]*")

	if [ -z "$library2" ] ; then
		echo -e "\n$directory does not contain properly formatted library name, e.g. \path\Lxxxx_Genos\n"
		exit
	fi
fi

#echo $library

echo -e "\nOperational parameters:"
echo -e "minimum ploidy: $ploidy_min"
echo -e "maximum ploidy: $ploidy_max"
echo -e "confidence threshold: $llr"
echo -e "minimum depth/locus: $depth"
echo -e "minumum genotyping rate: $geno"
echo -e "minimum reads on target: $cov"
echo -e "library name: $library2\n"


if [[ "$parallel" =~(y)$ ]]; then
	RSPATH=$(which ploidy_genotype_inference_shell_parallel.R)
	Rscript $RSPATH $ploidy_min $ploidy_max $llr $depth $geno $cov $working_dir $library2
else
	RSPATH=$(which ploidy_genotype_inference_shell.R)
	Rscript $RSPATH $ploidy_min $ploidy_max $llr $depth $geno $cov $working_dir $library2
fi


mv *.csv ../..
echo -e "\nLook for the following files:"
ls ../../*.csv
echo -e "\nGTseq_GenoCompile_ploidy_v1.sh script completed\n"