#!/usr/bin/perl -w

## 
## SOLiD Whole Genome Sequencing SNV calling Pipeline
## 
## Software update:
##                  GATK 2.7-2
##                  Picard 1.98
##                  samtools 0.1.98
##
## Reference: /mnt/mibNGSapps/Reference/humanDna.fasta
## Known site: /mnt/mibNGSapps/Reference/hg19.dbsnp.vcf
## 
## server: Cluster server
## 
##



use Getopt::Std;

sub Usage{
    print STDERR <<EOF;
    SOLiD Exome Seq. SNP calling Pipeline script
    
    -r reference (path: /mnt/mibNGSapps/Reference/humanDna.fasta)
    -k known site (path: /mnt/mibNGSapps/Reference/hg19.dbsnp.vcf)
    -b bam file prefix
    -h Help
    
EOF
    exit;
}
my %opt;
getopt("r:k:b:h:", \%opt);
my $ref = $opt{r} or &Usage();
my $knownsite = $opt{k} or &Usage();
my $prefix = $opt{b} or &Usage();
$opt{h} and &Usage();

##############################
print "Check SOLiD BAM file index";
`/mnt/mibNGSapps/opt/samtools-0.1.19/samtools index "$prefix"_Merge.bam`;

## GATK RealignerTargetCreator
## Emits intervals for the Local Indel Realigner to target for realignment.
`java -Xmx2g -jar /mnt/mibNGSapps/opt/GenomeAnalysisTK-2.7-2/GenomeAnalysisTK.jar -T RealignerTargetCreator -log $prefix.RealignerTargetCreator.log -nt 4 -R $ref -I "$prefix"_Merge.bam -known $knownsite -o "$prefix"_Merge.bam.intervals`;

## GATK IndelRealigner
## Performs local reaglignment of reads to correct misalignments due to the presence of indels.
`java -Xmx2g -jar /mnt/mibNGSapps/opt/GenomeAnalysisTK-2.7-2/GenomeAnalysisTK.jar -T IndelRealigner -log $prefix.IndelRealigner.log -R $ref -I "$prefix"_Merge.bam -targetIntervals "$prefix"_Merge.bam.intervals -known $knownsite -o $prefix.realigned.bam`;

## Picard Mark Duplicates
`java -Xmx2g -jar /mnt/mibNGSapps/opt/picard-tools-1.98/MarkDuplicates.jar I=$prefix.realigned.bam O=$prefix.realigned.rmdup.bam METRICS_FILE=$prefix.duplicate_report.txt VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=false ASSUME_SORTED=true`;

## SAMtools BAM file index
`/mnt/mibNGSapps/opt/samtools-0.1.19/samtools index $prefix.realigned.rmdup.bam`;

## SAMtools BAM file stat
`/mnt/mibNGSapps/opt/samtools-0.1.19/samtools flagstat $prefix.realigned.rmdup.bam > $prefix.realigned.rmdup.bam.stat`;

## GATK BaseRecalibrator
## First pass of the base quailty score recalibration -- Generates recalibration table based on various user-specified covariates.
## (such as read group, reported quality score, machine cycle, and nucleotide context.)
`java -Xmx2g -jar /mnt/mibNGSapps/opt/GenomeAnalysisTK-2.7-2/GenomeAnalysisTK.jar -T BaseRecalibrator -log $prefix.BaseRecalibrator.log -nt 4 -R $ref -I $prefix.realigned.rmdup.bam -knownSites $knownsite -o $prefix.recal_data.table --covariate QualityScoreCovariate --covariate ReadGroupCovariate --covariate ContextCovariate --covariate CycleCovariate --solid_nocall_strategy PURGE_READ --solid_recal_mode SET_Q_ZERO_BASE_N`;

## GATK PrintReads
## Renders, in SAM/BAM format, all reads from the input data set in the order in which they appear in the input file.
`java -Xmx2g -jar /mnt/mibNGSapps/opt/GenomeAnalysisTK-2.7-2/GenomeAnalysisTK.jar -T PrintReads -log $prefix.PrintReads.log -R $ref -I $prefix.realigned.rmdup.bam -BQSR $prefix.recal_data.table -o $prefix.realigned.rmdup.recali.bam`;

## GATK UinfiedGenotyper
## A variant caller which unifies the approaches of several disparate callers -- Works for single-sample and multi-sample data.
`java -Xmx2g -jar /mnt/mibNGSapps/opt/GenomeAnalysisTK-2.7-2/GenomeAnalysisTK.jar -T UnifiedGenotyper -log $prefix.UnifiedGenotyper.log -nt 4 -R $ref -I $prefix.realigned.rmdup.recali.bam -o "$prefix"_gatk.vcf -stand_call_conf 30.0 -stand_emit_conf 10.0 -glm both -D $knownsite`;

##############################
## SAMtools mpileup
`/mnt/mibNGSapps/opt/samtools-0.1.19/samtools mpileup -ugf $ref $prefix.realigned.rmdup.recali.bam | /mnt/mibNGSapps/opt/samtools-0.1.19/bcftools/bcftools view -bcvg - > $prefix.samtools.var.raw.bcf`;
`/mnt/mibNGSapps/opt/samtools-0.1.19/bcftools/bcftools view $prefix.samtools.var.raw.bcf | /mnt/mibNGSapps/opt/samtools-0.1.19/bcftools/vcfutils.pl varFilter -D 500 > "$prefix"_samtools.vcf`;

