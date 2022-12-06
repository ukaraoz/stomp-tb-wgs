#! /bin/bash

conda activate java

BASE=/global/homes/u/ukaraoz/cscratch
REFERENCE=${BASE}"/mtb/ref/MTB_ancestor_reference.fna"
MTBLIB=${BASE}"/mtb/mtb_lib"

BINDIR="/global/homes/u/ukaraoz/bin"
TRIMMOMATIC="/global/homes/u/ukaraoz/bin/Trimmomatic-0.39/trimmomatic-0.39.jar"
#KRAKENDB=/global/homes/u/ukaraoz/cscratch/kraken2db

peppe_phage_repeats2discardFILE=${MTBLIB}"/modules/ref/peppe_phage_repeats2discard.bed"
resistance_genes_listFILE=${MTBLIB}"/modules/ref/resistance_genes_list_03092021.bed"
GENOMEFILE=${BASE}"/mtb/mtb_lib/modules/ref/MTB_ancestor_reference.genome"
indelwindow=10
SNPEFFDATADIR=${BASE}"/mtb/ref/snpEff/data"

KRAKEN=${BINDIR}"/kraken2/kraken2"
BWA=${BINDIR}"/bwa"

R="/global/u2/u/ukaraoz/bin/anaconda3/bin/Rscript"
SEQTK=${BINDIR}"/seqtk"
SAMTOOLS=${BINDIR}"/samtools"
PICARD=${BINDIR}"/picard.2.18.29.jar"
GENOMECOVERAGEBED=${BINDIR}"/bedtools2/bin/genomeCoverageBed"
VARSCAN=${BINDIR}"/VarScan.v2.3.7.jar"
GATK=${BINDIR}"/gatk-4.0.2.1/gatk"
BCFTOOLS=${BINDIR}"/bcftools"
VCF2BED=${BINDIR}"/bedops/vcf2bed"
BEDTOOLS=${BINDIR}"/bedtools2/bin/bedtools"
VCFTOOLS=${BINDIR}"/vcftools"
VCF2PHYLIP=${BINDIR}"/vcf2phylip-2.7/vcf2phylip.py"
SNPEFF=${BINDIR}"/snpEff/snpEff.jar"
LINEAGECALLER=${BASE}"/mtb/scripts/fast-lineage-caller-vcf.py"
LINEAGESNPS=${BASE}"/mtb/mtb_lib/modules/ref/lineage_snps/coll.tsv"

MTBLIBSCRIPT=${BASE}"/mtb/mtb_lib/Main.py"

DR_SCRIPT="/Users/ukaraoz/Work/MTB/lib/mtb_lib/modules/scripts/Calling_Subscript.DR.sh"
EPI_SCRIPT="/Users/ukaraoz/Work/MTB/lib/mtb_lib/modules/scripts/Calling_Subscript.EPI.sh"

RSCRIPT="/Users/ukaraoz/Work/MTB/lib/mtb_lib/modules/scripts/2_1_filtering_ThePipeline.R"
WDENSITY=10
WINDEL=0

SAMPLENAMESFILE=/global/homes/u/ukaraoz/cscratch/mtb/stomptb/reads/raw/batch1_ids.txtad
#SAMPLENAMESFILE=$1
echo "##contig=<ID=NC_000962.3,length=4411532>" > vcfheader_contig.txt
IFS=$'\n' read -d '' -r -a samples < ${SAMPLENAMESFILE}

#KRAKENDB_KEY="20200102"
#KRAKENDB_KEY="20201202"
#KRAKENDB_KEY="20210517"
#KRAKENDB_KEY="20220607"
KRAKENDB_KEY="20220926"
KRAKENDB=/global/homes/u/ukaraoz/m3523/kraken/krakendb-${KRAKENDB_KEY}

for SAMPLE in "${samples[@]}"
do
	echo ${SAMPLE}
	java -jar ${TRIMMOMATIC} PE -threads 30 -phred33 \
	${SAMPLE}_R1_001.fastq ${SAMPLE}_R2_001.fastq ${SAMPLE}.P1.clean.fastq ${SAMPLE}.U1.clean.fastq ${SAMPLE}.P2.clean.fastq ${SAMPLE}.U2.clean.fastq \
	TRAILING:10 SLIDINGWINDOW:25:20 MINLEN:50 2> ${SAMPLE}.trimmomatic.log

	${KRAKEN} --threads 40 --db ${KRAKENDB} --use-names --report ${SAMPLE}.${KRAKENDB}.kraken.report --paired ${SAMPLE}.P1.clean.fastq ${SAMPLE}.P2.clean.fastq > ${SAMPLE}.kraken.${KRAKENDB_KEY}
	grep -E "Mycobacterium tuberculosis complex \(taxid 77643\)|Mycobacterium tuberculosis \(taxid 1773\)" ${SAMPLE}.kraken.${KRAKENDB_KEY} | \
	cut -f2 > ${SAMPLE}.${KRAKENDB_KEY}.P12.filtered.readlist
	${SEQTK} subseq ${SAMPLE}.P1.clean.fastq ${SAMPLE}.${KRAKENDB_KEY}.P12.filtered.readlist > ${SAMPLE}.P1.filtered.fastq
	${SEQTK} subseq ${SAMPLE}.P2.clean.fastq ${SAMPLE}.${KRAKENDB_KEY}.P12.filtered.readlist > ${SAMPLE}.P2.filtered.fastq

	rgroup="@RG\\tID:"${SAMPLE}"\\tSM:"${SAMPLE}"_sm\\tPU:"${SAMPLE}"_pu\\tLB:"${SAMPLE}"_lb"
	${BWA} mem -t 40 -R ${rgroup} ${REFERENCE} ${SAMPLE}.P1.filtered.fastq ${SAMPLE}.P2.filtered.fastq | \
	awk '$1 ~ /^@/ || $5 == 60' | \
	samtools view -bt ${REFERENCE} - | \
	samtools sort -o ${SAMPLE}.sort.bam

	java -jar ${PICARD} MarkDuplicates I=${SAMPLE}.sort.bam O=${SAMPLE}.nd.sort.bam M=${SAMPLE}.dup.metrix \
	ASSUME_SORTED=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT
	cp ${SAMPLE}.sort.bam ${SAMPLE}.dup.sort.bam
	cp ${SAMPLE}.nd.sort.bam ${SAMPLE}.sort.bam

	${GENOMECOVERAGEBED} -ibam ${SAMPLE}.sort.bam -d -g ${REFERENCE} > ${SAMPLE}.coverage

	${SAMTOOLS} mpileup -AB -f ${REFERENCE} ${SAMPLE}.sort.bam > ${SAMPLE}.mpileup

	java -jar ${VARSCAN} pileup2snp ${SAMPLE}.mpileup \
	--min-coverage 3 --min-reads2 3 --min-freq-for-hom 0.9 --min-var-freq 0.05 > ${SAMPLE}.snp.pileup2snp

	java -jar ${VARSCAN} filter ${SAMPLE}.snp.pileup2snp \
	--p-value 0.01 --min-coverage 20 --min-reads2 20 --min-avg-qual 25 --min-strands2 2 --min-var-freq 0.90 > ${SAMPLE}.EPI.snp.nodensityfilter
	
	${SAMTOOLS} index ${SAMPLE}.sort.bam
	touch ${SAMPLE}".done"
done
