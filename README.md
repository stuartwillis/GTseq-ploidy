# GTseq-ploidy
Modified version of the GTseq pipeline (**Campbell et al. 2015**) for polysomic polyploids, as described in **Willis et al. (2020)** *Single nucleotide genotypes and ploidy estimates for ploidy variable species generated with massively parallel amplicon sequencing*

Scripts are provided to run the GT-seq pipeline from fastq files containing sequences from the GT-seq style amplicon-sequencing project.

An R script containing functions used in Willis et al. (submitted) is also provided, "ploidy_genotype_etc_from_genos.R". It is intended to be run piecemeal in e.g. Rstudio. Example genos input files containing 1) parents and offspring from a known cross containing ploidy variation and 2) samples of known but variable ploidy. "Genos" files are an intermediate product of the GT-seq pipeline containing probe (allele) counts for each locus, one file per individual. The script is designed to be run with the folder containing the "genos" files as the working directory.

Please list issues with the pipeline or R script in the Github tab. Please note that we do not have sufficient resources to provide support for basic Linux, perl, bash, R issues that are dataset-specific.

## AtrGTseq-PLOIDYv1 (325) ############################
###This code assumes that fastq files have been correctly demultiplexed prior to pipeline implementation. Genotype the samples using the following scripts. Make sure to verify the genotyping is complete before moving on. Use the command TOP and do not proceed until finished.


###Navicate to the folder containing the fastq files, eg.
`cd /data/GTseq_Libraries/NS_0000/AtrGTseqV2.0/Individuals/`

###run the shell script which creates a command to count reads in each fastq; 
Assumes that GTseq_Genotyper_v3_Tetra.pl and Atr325_ProbeSeqs.csv are in the current working directory or in the path

`sh Atr325_GTseq_GenotyperShellGenerator.sh > NS_0000_Atr325_Genotyper.sh`

###run the script just created which contains the commands; creates *.genos files
`sh NS_0000_Atr325_Genotyper.sh`

###Replace 'Lxxxx' with the appropriate library name
make directory for the genos files
`mkdir Lxxxx_Genos/`
###move the genos file to this directory
`mv i*genos Lxxxx_Genos/`
###change working directory to this genos directory
`cd Lxxxx_Genos/`
###see options for the ploidy estimation and genotyping script
`bash Atr_GTseq_GenoCompile_ploidy_v1.sh -h`
###run the script (here with default options, but using multiple cores)
`bash Atr_GTseq_GenoCompile_ploidy_v1.sh -p y`
###make directory for the negative control genos
`mkdir N_T_C/`
###move negative control genos to that directory
`mv *N_T_C*genos N_T_C/`
`cd N_T_C/`
###create genotypes of negative controls
`NTC_GTseq_GenoCompile_v3_Tetra.pl > Lxxxx_N_T_C.csv`
`cd ..`
###create summary files
`GTseq_SummaryFigures_v3_Tetra.py`
> '/data/GTseq_Libraries/NS_0000/AtrGTseqV2.0/Individuals/Lxxxx_Genos/'
> 'Lxxxx'
`convert *png Lxxxx.pdf`
`rm -f *png`

##References

**Campbell, N. R., Harmon, S. A., & Narum, S. R.** (2015). *Genotyping-in-Thousands by sequencing (GT-seq): A cost effective SNP genotyping method based on custom amplicon sequencing*. Molecular Ecology Resources, 15(4), 855–867. doi: 10.1111/1755-0998.12357

**Willis, S.C., T. Delomas, B. Parker, D. Miller, S. Narum** (2020) *Single nucleotide genotypes and ploidy estimates for ploidy variable species generated with massively parallel amplicon sequencing*. bioRxiv