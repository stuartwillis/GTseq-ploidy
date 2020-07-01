#!/bin/sh
# Atr325_GTseq_GenotyperShellGenerator.sh

ls *fastq | sed 's/i.*fastq/GTseq_Genotyper_v3_Tetra.pl Atr325_ProbeSeqs.csv & > &/' | sed 's/fastq$/genos \&/' 
