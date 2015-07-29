#!/bin/bash -e
if [ $# -lt 1 ]
then
	echo usage: $0 [bam]
	exit 1
fi

input=$1
output=${input/.bam/.FP.vcf}

/home/D131748/apps/JAVA/jre1.7.0_60/bin/java -jar /home/D131748/apps/ETC/GATK/GenomeAnalysis-3.3.0/GenomeAnalysisTK.jar\
	-T HaplotypeCaller\
	-R /home/D131748/Reference/Human/b37/human_g1k_v37.fasta\
	-I $input\
	--dbsnp /home/D131748/Reference/dbsnp/dbsnp_142.b37.vcf\
	-stand_call_conf 30\
	-stand_emit_conf 10\
	-L /home/D131748/Reference/Target_intervals/OncoPanel_v2/b37/OPv2_FingerPriting.interval_list\
	-o $output
