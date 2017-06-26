#!/usr/bin/perl
#Generate bash scripts to run on cluster for extended distortion panel

use strict; use warnings;

open (POP, "</ebio/abt6_projects7/ath_seg_distortion/data/sequence/2014-10-21_POP_ID_GM_all.txt") or die "Cannot open Population description file\n";


while (<POP>) {
	chomp $_;
	my @pop = split("\t",$_);
	#choose population to start with
	my $i = $pop[0];
	print $i,"\n";

	#choose P1 and P2 columns from the genome matrix
	my $P1 = $pop[1];
	my $P2 = $pop[2];
	print $P1."\t".$P2."\n";
	
	open (OUT, ">/ebio/abt6_projects7/ath_seg_distortion/data/analysis/shore_analysis_extended_samples/extended_cluster_cmds/$i.command.sh") or die "Cannot open output file\n";
	#make file to run pipeline on the cluster
	print OUT "#\\!/bin/bash"."\n";
	print OUT "#\$ -pe parallel 4"."\n";
	print OUT "#\$ -l h_vmem=20G"."\n";
	print OUT "#\$ -cwd"."\n";
	print OUT "#\$ -V"."\n"."\n";
	print OUT "source ~/.bashrc"."\n";
	print OUT "src=\$(q-src --debug /ebio/abt6_projects7/ath_seg_distortion/data/analysis/refseq) || exit 1"."\n";
	print OUT "echo \$src"."\n"."\n";
	print OUT "/ebio/abt6_projects7/ath_seg_distortion/code/snp_calling_pipe_cluster_expanded_panel.pl $i $P1 $P2 \$src/TAIR10_All_Chr.fa.shore"."\n";
	close (OUT)
}
close (POP);













