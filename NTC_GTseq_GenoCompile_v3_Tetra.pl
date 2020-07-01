#!/usr/bin/perl
#GenoCompile_v3_Tetra.pl by Nate Campbell
#Compile genotypes from individual genotype files from GTseq_Genotyper_v2 output ".genos" files...
#This version utilizes the expanded output from the GTseq_Genotyper_v2 script to gather summary data and does not
#require the individual fastq files.
#For optional output formats use arguments. N for numeric genotypes or C for allele counts or F for frequencies or G for GenAlEx(ish) or AB for AABB format.  
#Defaults to SNP genotypes.
#Optional filtered output requires 2 argument values.  Genotype output type and a genotyping threshold [S,N,C,F,G, or AB] [90]
# example: $ GTseq_GenoCompile_v2.pl S 90 (outputs SNP genotypes for individual .genos files with 90% or higher genotyping percentage)
# genotypes for individuals with less than the threshold genotyping percentage are converted to "0000".

use strict; use warnings;

my $flag = "S";
my $geno_thresh = 0;
if (@ARGV == 1) {$flag = $ARGV[0];}
if (@ARGV == 2) {$flag = $ARGV[0]; $geno_thresh = $ARGV[1];}

my %Ave_Freq = (); #declare hash for the average allele frequency at each locus...
my %Num_Geno = (); #declare hash for the number of true genotypes retained for each locus...

my @Files = `ls *genos`;
chomp ( @Files );
print "Sample,Raw Reads,On-Target Reads,%On-Target,%GT,IFI";

open (FILE1, "<$Files[0]") or die;
	while (<FILE1>) {
		if ($. > 1){
			my @info1 = split ",", $_;
			my $assay1 = $info1[0];
			$Ave_Freq{$assay1} = 0;  #initialize average frequency to zero for each locus.
			$Num_Geno{$assay1} = 0;  #initialize number of returned genotypes at each locus.
			if ($flag =~ m/G/) {print ",$assay1,,,";}
			else {print ",$assay1";}
			}
		}
close FILE1;

foreach my $genos (@Files) {
	open (FILE2, "<$genos") or die "Error opening $genos\n";
		while (<FILE2>) {
			my $freq = 0;
			my @junk = split ",", $_;
			if (($. > 1) && ($junk[5] !~ m/NA/)) {
				$Num_Geno{$junk[0]}++;
				if ($junk[5] =~ m/A1HOM/) {$freq = 0.00;}
				elsif ($junk[5] =~ m/HET1/) {$freq = 0.25;}
				elsif ($junk[5] =~ m/HET2/) {$freq = 0.50;}
				elsif ($junk[5] =~ m/HET3/) {$freq = 0.75;}
				elsif ($junk[5] =~ m/A2HOM/) {$freq = 1.00;}
				$Ave_Freq{$junk[0]} = $Ave_Freq{$junk[0]} + $freq;
				}
			}
			close FILE2;
		}

foreach my $assays (keys %Ave_Freq) {
	if ($Num_Geno{$assays} == 0) {$Num_Geno{$assays} = 0.01}
	$Ave_Freq{$assays} = $Ave_Freq{$assays}/$Num_Geno{$assays};
	}

print "\n";

foreach my $samples (@Files) {
	my $raw = 0;
	my $on_target = 0;
	my $GT_pct = 0;
	my $num_targets = 0;
	my $IFI = 0;
	my $sample_name = $samples;
	$sample_name =~ s/.genos//;
	print "$sample_name,";
	open (FILE, "<$samples") or die;
	while (<FILE>) {
	if ($. == 1) {chomp; my @summary = split ",", $_; $summary[1] =~ s/Raw-Reads\://; print "$summary[1],"; $raw = $summary[1];
		$IFI = $summary[4]; $IFI =~ s/IFI_score\://;}
	elsif ($. > 1) {
		$num_targets++;
		chomp;
		my @info = split ",", $_;
		my $assay = $info[0];
		my $geno = $info[4];
		if ($geno =~ m/NA|0000/) {$GT_pct++}
		my $count1 = $info[1];
		$count1 =~ s/.*=//;
		my $count2 = $info[2];
		$count2 =~ s/.*=//;
		$on_target = $on_target + $count1 + $count2;
				}
			}
		close FILE;
		$GT_pct = 100-($GT_pct/$num_targets*100);
		my $OT_pct = $on_target/$raw*100;
		my $Print_GT_pct = sprintf("%.2f", $GT_pct);
		$OT_pct = sprintf("%.2f", $OT_pct);
		print "$on_target,$OT_pct,$Print_GT_pct,$IFI,";
		
	open (FILE, "<$samples") or die;
	while (<FILE>) {
		if ($. > 1) {
		chomp;
		my @info1 = split ",", $_;
		my $geno = $info1[4];
		my $L_count = 0;
		$info1[1] =~ s/.*=//;
		$info1[2] =~ s/.*=//;
		$L_count = $info1[1] + $info1[2];
		my $NumGT = "0000";
		my $AB_GT = "0000";
		my $FGT = 9;
		if($info1[5] =~ m/A1HOM/) {$NumGT = "1111"; $FGT = sprintf("%.2f", 0.001); $AB_GT = "AAAA";}  #Frequency not exactly zero to avoid errors with log transformations in analyses such as PCA.
		elsif($info1[5] =~ m/HET1/) {$NumGT = "1112"; $FGT = sprintf("%.2f",0.250); $AB_GT = "AAAB";}
		elsif($info1[5] =~ m/HET2/) {$NumGT = "1122"; $FGT = sprintf("%.2f",0.500); $AB_GT = "AABB";}
		elsif($info1[5] =~ m/HET3/) {$NumGT = "1222"; $FGT = sprintf("%.2f",0.750); $AB_GT = "ABBB";}
		elsif($info1[5] =~ m/A2HOM/) {$NumGT = "2222"; $FGT = sprintf("%.2f",1.000); $AB_GT = "BBBB";}
		elsif($info1[5] =~ m/NA/) {$NumGT = "0000"; $FGT = sprintf("%.2f",$Ave_Freq{$info1[0]}); $AB_GT = "0000";}
		
		if(($flag =~ m/S/) && ($GT_pct >= $geno_thresh)) {print "$geno,";}
		elsif(($flag =~ m/S/) && ($GT_pct < $geno_thresh)) {print "0000,";}
		elsif($flag =~ m/C/) {print "$L_count,";}
		elsif(($flag =~ m/N/) && ($GT_pct >= $geno_thresh)) {print "$NumGT,";}
		elsif(($flag =~ m/N/) && ($GT_pct < $geno_thresh)) {print "0000,";}
		elsif($flag =~ m/F/) {print "$FGT,";}
		elsif(($flag =~ m/AB/) && ($GT_pct >= $geno_thresh)) {print "$AB_GT,";}
		elsif(($flag =~ m/G/) && ($GT_pct >= $geno_thresh)) {my @alleles = split "", $NumGT; print "$alleles[0],$alleles[1],$alleles[2],$alleles[3],";}
		elsif(($flag =~ m/G/) && ($GT_pct < $geno_thresh)) {print "0,0,0,0,";}
		else {print "0000,";}
				}
			}
		print "\n";
		close FILE;
	}
