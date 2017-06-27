#READ ME

##########
#The following script takes a shore formated raw read file and aligns it to the Col-0 reference genome and calls SNPs. The script then filters those SNPs for those previously discovered in the parental resequencing project (Cao, 2011). The allele frequency of the maternal parent is outputed

#01.snp_calling_pipe_cluster_expanded_panel.pl
#usage: snp_calling_pipeline.pl <POPID> <P1 column> <P2 column> <index directory>

#POPID and P1 and P2 column ids are found in 2014-10-21_POP_ID_GM_all.txt (provided here).
#Column ids refer to the position of the parental genotype in the genome matrix 80x80_TAIR10_SNPs_split.txt (also provided here)

##########
#The following script can be used to generate cluster submission scripts for each sample.
#02.generate_cluster_cmds_extended.pl

##########
#The R script in this folder will model the allele frequencies in specific window sizes.
#2015-11-20_BB_modeling_of_allele_freq.R

#########
#Output allele frequencies for each sample are located in the allele_freq folder. The header information is as follows:
#POPID	Chr	Pos	Ref	Alt(Maternal)	Mapping_Quality	Coverage(Maternal) Frequency(Maternal)	SNP_Quality

#########
#Output allele frequencies for bulked segregant analysis are in the allele_freq_BSA folder. Header information is the same as above. csv files indicate which populations were combined for each pool.