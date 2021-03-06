#!/usr/bin/perl -w

## 
## SOLiD Exome Seq. Annotation pipeline
## Software:
##          bgzip, tabix, vcf-merge, annovar
##
##
##


use Getopt::Std;

sub Usage{
	print STDERR <<EOF;
	SOLiD Exome Seq. Annotation Pipeline script

	annotation_pipeline_for_SOLiD.pl [arguments]

	-n  sample name
	-p  output path
	-a  annvoar input file (.avinput)
	-db annovar human database (default:/opt/annovar_2013feb21/humandb/)
	-h  help

	Example:
            annotation_pipeline_for_SOLid.pl -n 5132661NT -p /mnt/mibNGS_abi_bio_ExonCaptures/Analysis/Running/Sample_5132661NT/5132661NT \
            -a /mnt/mibNGS_abi_bio_ExonCaptures/Analysis/Running/Sample_5132661NT/5132661NT.avinput -db
	Version: 2013-04-01
EOF
exit;
}

my %opt;
getopt("n:p:a:db:h:", \%opt);
my $sample_name = $opt{n} or &Usage();
my $output_path = $opt{p} or &Usage();
my $avinput_file = $opt{a} or &Usage();
my $avdb = $opt{db} or &Usage();
$opt{h}  and &Usage();

## vcf files to compress use bgzip
if(-e '"$sample_name"_gatk.vcf.gz' and '"$sample_name"_samtools.vcf.gz'){
	next;
}else{
	`bgzip "$sample_name"_gatk.vcf`;
	`bgzip "$sample_name"_samtools.vcf`;
}

## creat index use tabix
if(-e '"$sample_name"_gatk.vcf.gz.tbi' and '"$sample_name"_samtools.vcf.gz.tbi'){
	next;
}else{
	`tabix -p vcf "$sample_name"_gatk.vcf.gz`;
	`tabix -p vcf "$sample_name"_samtools.vcf.gz`;
}

## VCFtools -> vcf-merge (gatk and samtools vcf merge)
if(-e "$sample_name.vcf"){
	next;
}else{
	`vcf-merge "$sample_name"_gatk.vcf.gz "$sample_name"_samtools.vcf.gz > $sample_name.vcf`;
}

## Annovar -> convert to annovar input format
if(-e "$sample_name.avinput"){
	next;
}else{
	`convert2annovar.pl -format vcf4 --includeinfo $sample_name.vcf > $sample_name.avinput`;
}

## Annovar -> summarize annovar
##            hg19
##            dbsnp137
##            ucsc knowngene
##            ESP 6500
##            1KG 1000g2012apr
##
## website: http://www.openbioinformatics.org/annovar/

if(-e "$sample_name.exome_summary.csv" and "$sample_name.genome_summary.csv"){
	next;
}else{
	`summarize_annovar.pl -buildver hg19 --verdbsnp 137 --ver1000g 1000g2012apr --veresp 6500 --genetype knowngene --outfile $output_path $avinput_file $avdb`;
}
