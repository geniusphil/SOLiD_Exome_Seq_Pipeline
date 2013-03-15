#!/usr/bin/perl -w

## 
## SOLiD Exome Seq. SNP calling Pipeline
## 
## Software update:
##                  GATK 2.3-9
##                  Picard 1.85
##                  samtools 0.1.18
##
## reference: /mnt/mibNGS_abi_bio_ExonCaptures/Analysis/Reference/humanDna.fasta
## known site: /mnt/mibNGS_abi_bio_ExonCaptures/Analysis/Reference/hg19.dbsnp.vcf
## server: 72.69
##



use Getopt::Std;

sub Usage{
    print STDERR <<EOF;
    SOLiD Exome Seq. SNP calling Pipeline script
    
    -r reference
    -b bam file prefix
    -h Help
    
EOF
    exit;
}
my %opt;
getopt("r:b:h:", \%opt);
my $ref = $opt{r} or &Usage();
my $prefix = $opt{b} or &Usage();
$opt{h} and &Usage();

## GATK RealignerTargetCreator
`java -jar /opt/GenomeAnalysisTK-2.3-9/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt 4 -R $ref -I $prefix.ma.bam -o $prefix.ma.bam.intervals`;
## GATK IndelRealigner
`java -jar /opt/GenomeAnalysisTK-2.3-9/GenomeAnalysisTK.jar -T IndelRealigner -R $ref -I $prefix.ma.bam -targetIntervals $prefix.ma.bam.intervals -o $prefix.realigned.bam`;

## Picard Make duplicates
`java -jar /opt/picard-tools-1.85/picard-tools-1.85/MarkDuplicates.jar I=$prefix.realigned.bam O=$prefix.realigned.rmdup.bam METRICS_FILE=$prefix.duplicate_report.txt VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true ASSUME_SORTED=true`;

## SAMtools BAM file index
`samtools index $prefix.realigned.rmdup.bam`;
## SAMtools BAM file stat
`samtools flagstat $prefix.realigned.rmdup.bam > $prefix.realigned.rmdup.bam.stat`;

## GATK BaseRecalibrator
`java -jar /opt/GenomeAnalysisTK-2.3-9/GenomeAnalysisTK.jar -T BaseRecalibrator -nt 4 -R $ref -I $prefix.realigned.rmdup.bam -knownSites ../../hg19.dbsnp.vcf -o $prefix.recal_data.grp --covariate QualityScoreCovariate --covariate ReadGroupCovariate --covariate ContextCovariate --covariate CycleCovariate --solid_nocall_strategy PURGE_READ --solid_recal_mode SET_Q_ZERO_BASE_N`;
## GATK PrintReads
`java -jar /opt/GenomeAnalysisTK-2.3-9/GenomeAnalysisTK.jar -T PrintReads -R $ref -I $prefix.realigned.rmdup.bam -BQSR $prefix.recal_data.grp -o $prefix.realigned.rmdup.recali.bam`;
## GATK UinfiedGenotyper
`java -jar /opt/GenomeAnalysisTK-2.3-9/GenomeAnalysisTK.jar -T UnifiedGenotyper -nt 4 -R $ref -I $prefix.realigned.rmdup.recali.bam -o "$prefix"_gatk.vcf -stand_call_conf 30.0 -stand_emit_conf 10.0 -glm both -D ../../hg19.dbsnp.vcf`;

## SAMtools mpileup
`samtools mpileup -ugf $ref $prefix.realigned.rmdup.recali.bam | bcftools view -bcvg - > $prefix.samtools.var.raw.bcf`;
`bcftools view $prefix.samtools.var.raw.bcf | vcfutils.pl varFilter -D 500 > "$prefix"_samtools.vcf`;

