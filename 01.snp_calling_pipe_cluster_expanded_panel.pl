#!/usr/bin/perl
#Build workflow to analyze data from F2 pools. The steps are outlined below:
#1. Align reads from each populations to the reference Col-0 genome
#2. Call variant positions using shoremap
#3. Subset variant bases for sites with high confidence (ie. previously identified in the 80 genomes resequencing project.
#Note - all paths are specific to DKS folder structure
#Note - Must run shore import before running this script - Use separate script to trim restriction enzyme site from zipped (.xz) shore format

use strict; use warnings;
die "usage: snp_calling_pipeline.pl <POPID> <P1 column> <P2 column> <index directory>" unless @ARGV==4;

#####
#build index for alignment
/ebio/abt6/dseymour/APPLICATIONS/shore-0.9.0/shore -T /ebio/abt6_projects7/ath_seg_distortion/tmp/ preprocess -f /ebio/abt6_projects7/ath_seg_distortion/data/analysis/refseq/TAIR10_All_Chr.fa -i /ebio/abt6_projects7/ath_seg_distortion/data/analysis/refseq

#/ebio/abt6/dseymour/APPLICATIONS/shore-0.9.0/shore -T /ebio/abt6_projects7/ath_seg_distortion/tmp/ preprocess --upgrade -g -U /ebio/abt6_projects7/ath_seg_distortion/data/analysis/refseq/TAIR10_All_Chr.fa.shore


#####
#choose population
my $i = $ARGV[0]."_trim";

#choose P1 and P2 columns from the genome matrix
my $P1 = $ARGV[1];
my $P2 = $ARGV[2];

my $indexdir = $ARGV[3];



#####
#MAPPING
#First put reads into its own directory in order to keep track of mapping outputs
`mkdir /ebio/abt6_projects7/ath_seg_distortion/data/analysis/shore_analysis_extended_samples/shore_trim_extended_samples/$i`;

`mv /ebio/abt6_projects7/ath_seg_distortion/data/analysis/shore_analysis_extended_samples/shore_trim_extended_samples/$i.shore /ebio/abt6_projects7/ath_seg_distortion/data/analysis/shore_analysis_extended_samples/shore_trim_extended_samples/$i/`;

#use shore to map files to TAIR10 reference
`/ebio/abt6/dseymour/APPLICATIONS/shore-0.9.0/shore -T /ebio/abt6_projects7/ath_seg_distortion/tmp mapflowcell -i $indexdir -f /ebio/abt6_projects7/ath_seg_distortion/data/analysis/shore_analysis_extended_samples/shore_trim_extended_samples/$i/$i.shore -c 4 -n 2.5% -g 2.5% -s 500`;

print "Alignment complete for $i\n";

#move resulting output files and rename
`mv /ebio/abt6_projects7/ath_seg_distortion/data/analysis/shore_analysis_extended_samples/shore_trim_extended_samples/$i/map.list.xz /ebio/abt6_projects7/ath_seg_distortion/data/analysis/shore_analysis_extended_samples/shore_align_extended_samples/$i.map.list.xz`;

print "Mapping files moved to the proper directories\n";




#####
#use shore to call variant bases

`/ebio/abt6/dseymour/APPLICATIONS/shore-0.9.0/shore -T /ebio/abt6_projects7/ath_seg_distortion/tmp qVar -n $i -f $indexdir -i /ebio/abt6_projects7/ath_seg_distortion/data/analysis/shore_analysis_extended_samples/shore_align_extended_samples/$i.map.list.xz -s /ebio/abt6/dseymour/APPLICATIONS/shore-0.9.0/scoring_matrices/scoring_matrix_het.txt -o /ebio/abt6_projects7/ath_seg_distortion/data/analysis/shore_analysis_extended_samples/shore_variant_dir_extended_samples/$i -e -a 0.2 -x`;

print "SNP calling complete for $i\n";




#####
#concatenate reference.shore and snp.shore and filter for a concorcance less than 1
`xz -d -k /ebio/abt6_projects7/ath_seg_distortion/data/analysis/shore_analysis_extended_samples/shore_variant_dir_extended_samples/$i/ConsensusAnalysis/reference.shore.xz`;

`cat /ebio/abt6_projects7/ath_seg_distortion/data/analysis/shore_analysis_extended_samples/shore_variant_dir_extended_samples/$i/ConsensusAnalysis/reference.shore /ebio/abt6_projects7/ath_seg_distortion/data/analysis/shore_analysis_extended_samples/shore_variant_dir_extended_samples/$i/ConsensusAnalysis/snp.shore > /ebio/abt6_projects7/ath_seg_distortion/data/analysis/shore_analysis_extended_samples/shore_variant_dir_extended_samples/$i/ConsensusAnalysis/all_sites.shore`;

`awk \'\$8<1\' /ebio/abt6_projects7/ath_seg_distortion/data/analysis/shore_analysis_extended_samples/shore_variant_dir_extended_samples/$i/ConsensusAnalysis/all_sites.shore > /ebio/abt6_projects7/ath_seg_distortion/data/analysis/shore_analysis_extended_samples/shore_variant_dir_extended_samples/$i/ConsensusAnalysis/variable_sites_temp.shore`;

`sort -k2,2d -k3,3n /ebio/abt6_projects7/ath_seg_distortion/data/analysis/shore_analysis_extended_samples/shore_variant_dir_extended_samples/$i/ConsensusAnalysis/variable_sites_temp.shore > /ebio/abt6_projects7/ath_seg_distortion/data/analysis/shore_analysis_extended_samples/shore_variant_dir_extended_samples/$i/ConsensusAnalysis/variable_sites.shore`;

`rm /ebio/abt6_projects7/ath_seg_distortion/data/analysis/shore_analysis_extended_samples/shore_variant_dir_extended_samples/$i/ConsensusAnalysis/reference.shore`;

`rm /ebio/abt6_projects7/ath_seg_distortion/data/analysis/shore_analysis_extended_samples/shore_variant_dir_extended_samples/$i/ConsensusAnalysis/all_sites.shore`;

`rm /ebio/abt6_projects7/ath_seg_distortion/data/analysis/shore_analysis_extended_samples/shore_variant_dir_extended_samples/$i/ConsensusAnalysis/variable_sites_temp.shore`;

print "Variant SNP file produced successfully\n";




#######
#Identify SNP positions called in the 80 resequencing genome matrix between the two parents
open (GM, "</ebio/abt6_projects7/ath_seg_distortion/data/sequence/80x80_TAIR10_SNPs_split.txt") or die "Cannot open Genome Matrix\n";

my %F2;

while (<GM>) {
	my @temp = split("\t", $_);
	#print $temp[42]."\t".$temp[36]."\n";
	if (($temp[$P1]=~ m/[ATGC]/) && ($temp[$P2]=~ m/[ATGC]/) && ($temp[$P1] ne $temp[$P2])) {
		$F2{"Chr".$temp[0]."_".$temp[1]} = $temp[$P1].$temp[$P2]
	}
}

#foreach my $a (keys %F2) {
#	print $a."\t".$F2{$a}."\n";
#}

my $snpnum = scalar keys %F2;

close(GM);

print "Parental SNPs identified: $snpnum\n";




#########
#Subset snp calls in RAD data by those SNPs idenified between two parents in the genome matrix. Pick the allele frequency for Parent 1
open (SHORE, "</ebio/abt6_projects7/ath_seg_distortion/data/analysis/shore_analysis_extended_samples/shore_variant_dir_extended_samples/$i/ConsensusAnalysis/variable_sites.shore") or die "Cannot open Variable sites file\n";

open (OUT, ">/ebio/abt6_projects7/ath_seg_distortion/data/analysis/shore_analysis_extended_samples/shore_allelefreq_extended_samples/$i.P1_allele_freq.txt") or die "Cannot open Output file\n";

my $count=0;
my $count1=0;
my $count2=0;
my $count3=0;
while (<SHORE>) {
	my @shorepos = split("\t", $_);
	if (exists $F2{$shorepos[1]."_".$shorepos[2]}) {
		$count++;
		if ($shorepos[4] eq substr($F2{$shorepos[1]."_".$shorepos[2]},0,1)) {
			$count1++;
			if($shorepos[7] ne 1) {
				$count2++;
				if($shorepos[8] eq 1) {
					$count3++;
					print OUT $_;
				}
			}
		}
	}
}

print "$count SNPs found\n";
print "$count1 SNPs for P1 allele\n";
print "$count2 SNPs less than 100% frequency\n";
print "$count3 good SNPs identified\n";

close (SHORE);
close (OUT);

print "Pipeline complete\n";
