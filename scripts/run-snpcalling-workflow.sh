#! /bin/bash

conda activate java

BASE=/global/homes/u/ukaraoz/cscratch
REFERENCE=${BASE}"/mtb/stomptb/ref/GCF_000195955.2_ASM19595v2_genomic.fna"
MTBLIB=${BASE}"/mtb/mtb_lib"

BINDIR="/global/homes/u/ukaraoz/bin"
KRAKENDB=/global/homes/u/ukaraoz/cscratch/kraken2db
peppe_phage_repeats2discardFILE=${MTBLIB}"/modules/ref/peppe_phage_repeats2discard.bed"
resistance_genes_listFILE=${MTBLIB}"/modules/ref/resistance_genes_list_03092021.bed"
GENOMEFILE=${BASE}"/mtb/mtb_lib/modules/ref/H37Rv.genome"
indelwindow=10
SNPEFFDATADIR=${BASE}"/mtb/ref/snpEff/data"

KRAKEN=${BINDIR}"/kraken2/kraken2"
BWA=${BINDIR}"/bwa"

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

#SAMPLENAMESFILE=${BASE}/philippines_ids.txt
SAMPLENAMESFILE=$1
echo "##contig=<ID=NC_000962.3,length=4411532>" > vcfheader_contig.txt

IFS=$'\n' read -d '' -r -a samples < ${SAMPLENAMESFILE}
for SAMPLE in "${samples[@]}"
do
	echo ${SAMPLE}
	${KRAKEN} --threads 16 --db ${KRAKENDB} --use-names --report ${SAMPLE}.kraken.report --paired ${SAMPLE}.P1.clean.fastq ${SAMPLE}.P2.clean.fastq > ${SAMPLE}.kraken
	grep -E "Mycobacterium tuberculosis complex \(taxid 77643\)|Mycobacterium tuberculosis \(taxid 1773\)" ${SAMPLE}.kraken | \
	cut -f2 > ${SAMPLE}.P12.filtered.readlist
	${SEQTK} subseq ${SAMPLE}.P1.clean.fastq ${SAMPLE}.P12.filtered.readlist > ${SAMPLE}.P1.filtered.fastq
	${SEQTK} subseq ${SAMPLE}.P2.clean.fastq ${SAMPLE}.P12.filtered.readlist > ${SAMPLE}.P2.filtered.fastq

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

	# Standard
	# 1. SNP
	java -jar ${VARSCAN} mpileup2snp ${SAMPLE}.mpileup \
	--min-coverage 3 --min-reads2 3 --min-freq-for-hom 0.9 --min-var-freq 0.05 --output-vcf 1 > ${SAMPLE}.snp.vcf
	
	## 1.1. filter DR - these are vSNPs
	java -jar ${VARSCAN} filter ${SAMPLE}.snp.vcf \
	--p-value 0.01 --min-coverage 10 --min-reads2 6 --min-avg-qual 25 --min-strands2 2 --min-var-freq 0.10 > ${SAMPLE}.DR.snp.final.vcf
	
	## 1.2. filter EPI - these are fixedSNPs
	java -jar ${VARSCAN} filter ${SAMPLE}.snp.vcf \
	--p-value 0.01 --min-coverage 20 --min-reads2 20 --min-avg-qual 25 --min-strands2 2 --min-var-freq 0.90 > ${SAMPLE}.EPI.snp.nodensityfilter.vcf
	
	# 2. indel1
	java -jar ${VARSCAN} mpileup2indel ${SAMPLE}.mpileup \
	--min-coverage 10 --min-freq-for-hom 0.90 --min-var-freq 0.10 --output-vcf 1 > ${SAMPLE}.varscan.snp.indel.vcf
	
	# 3. indel2
	java -jar ${VARSCAN} mpileup2indel ${SAMPLE}.mpileup \
	--min-coverage 10 --min-freq-for-hom 0.90 --min-var-freq 0.10 --output-vcf 1 > ${SAMPLE}.DR.snp.indel.vcf
	
	# 4. IndexBAM 
	${SAMTOOLS} index ${SAMPLE}.sort.bam

	# 5. IndelCalling
	${GATK} HaplotypeCaller -R ${REFERENCE} -I ${SAMPLE}.sort.bam -O ${SAMPLE}.gatk.vcf -ploidy 1
	${GATK} SelectVariants -R ${REFERENCE} -V ${SAMPLE}.gatk.vcf -select-type INDEL -O ${SAMPLE}.gatk.indel.vcf
	${GATK} VariantFiltration -R ${REFERENCE} -V ${SAMPLE}.gatk.indel.vcf \
	--filter-expression "QD<2.0" --filter-name "LowQD" \
	--filter-expression "FS>200.0" --filter-name "HighFS" \
	--filter-expression "SOR>3.0" --filter-name "HighSOR" \
	--filter-expression "MQ<40.0" --filter-name "Low40MQ" \
	--filter-expression "ReadPosRankSum<-20.0" --filter-name "LowReadPRS" \
	--filter-expression "DP<20" --filter-name "LowDepth" \
	-O ${SAMPLE}.gatk.indel.filter.vcf

	${MTBLIBSCRIPT} coverage -e bam -p ${SAMPLE}

	### change sample name from generic to actual
	echo ${SAMPLE} > sample.txt
	${BCFTOOLS} reheader -s sample.txt ${SAMPLE}.DR.snp.final.vcf | \
	${BCFTOOLS} annotate --header-lines vcfheader_contig.txt -O z -o ${SAMPLE}.DR.snp.final.vcf.gz -
	${BCFTOOLS} index ${SAMPLE}.DR.snp.final.vcf.gz

	# filtering vcfs by region
	## prepare bed files for filtering
	### concat indels from varscan and gatk, create windows around it, and turn to bed
	### change sample name from generic to actual
	${BCFTOOLS} reheader -s sample.txt ${SAMPLE}.DR.snp.indel.vcf -o temp.vcf; cp temp.vcf ${SAMPLE}.DR.snp.indel.vcf; rm temp.vcf
	${BCFTOOLS} reheader -s sample.txt ${SAMPLE}.gatk.indel.filter.vcf -o temp.vcf; cp temp.vcf ${SAMPLE}.gatk.indel.filter.vcf; rm temp.vcf
	${BCFTOOLS} concat ${SAMPLE}.DR.snp.indel.vcf ${SAMPLE}.gatk.indel.filter.vcf -o ${SAMPLE}.all.indel.vcf
	${VCF2BED} < ${SAMPLE}.all.indel.vcf > temp.bed; ${BEDTOOLS} slop -i temp.bed -b ${indelwindow} -g ${GENOMEFILE} > ${SAMPLE}.all.indel.bed; rm temp.bed

	## exclude regions
	[ -e ${SAMPLE}.filteredregions.log.txt ] && rm ${SAMPLE}.filteredregions.log.txt
	### exclude repeats etc.
	${VCFTOOLS} --vcf ${SAMPLE}.EPI.snp.nodensityfilter.vcf --exclude-bed ${peppe_phage_repeats2discardFILE} --recode --out ${SAMPLE}.EPI.snp.nodensityfilter.excluderepeats 2>>${SAMPLE}.filteredregions.log.txt
	### exclude resistance genes
	${VCFTOOLS} --vcf ${SAMPLE}.EPI.snp.nodensityfilter.excluderepeats.recode.vcf --exclude-bed ${resistance_genes_listFILE} --recode --out ${SAMPLE}.EPI.snp.nodensityfilter.excluderepeats_drugres 2>>${SAMPLE}.filteredregions.log.txt
	### exclude around indels
	${VCFTOOLS} --vcf ${SAMPLE}.EPI.snp.nodensityfilter.excluderepeats_drugres.recode.vcf --exclude-bed ${SAMPLE}.all.indel.bed --recode --out ${SAMPLE}.EPI.snp.nodensityfilter.excluderepeats_drugres_indels 2>>${SAMPLE}.filteredregions.log.txt
	### change sample name from generic to actual
	${BCFTOOLS} reheader -s sample.txt ${SAMPLE}.EPI.snp.nodensityfilter.excluderepeats_drugres_indels.recode.vcf | \
	${BCFTOOLS} annotate --header-lines vcfheader_contig.txt -O z -o ${SAMPLE}.EPI.snp.nodensityfilter.final.vcf.gz -
	${BCFTOOLS} index ${SAMPLE}.EPI.snp.nodensityfilter.final.vcf.gz

	## exclude regions
	# echo ${SAMPLE} > sample.txt
	[ -e ${SAMPLE}.filteredregions1.log.txt ] && rm ${SAMPLE}.filteredregions1.log.txt
	### exclude repeats etc.
	${VCFTOOLS} --vcf ${SAMPLE}.EPI.snp.nodensityfilter.vcf --exclude-bed ${peppe_phage_repeats2discardFILE} --recode --out ${SAMPLE}.EPI.snp.nodensityfilter.excluderepeats 2>>${SAMPLE}.filteredregions1.log.txt
	### exclude around indels
	${VCFTOOLS} --vcf ${SAMPLE}.EPI.snp.nodensityfilter.excluderepeats.recode.vcf --exclude-bed ${SAMPLE}.all.indel.bed --recode --out ${SAMPLE}.EPI.snp.nodensityfilter.excluderepeats_indels 2>>${SAMPLE}.filteredregions1.log.txt
	### change sample name from generic to actual
	${BCFTOOLS} reheader -s sample.txt ${SAMPLE}.EPI.snp.nodensityfilter.excluderepeats_indels.recode.vcf | \
	${BCFTOOLS} annotate --header-lines vcfheader_contig.txt -O z -o ${SAMPLE}.EPI.snp.nodensityfilter.excluderepeats_indels.final.vcf.gz -
	${BCFTOOLS} index ${SAMPLE}.EPI.snp.nodensityfilter.excluderepeats_indels.final.vcf.gz


	# lineage=$(${LINEAGECALLER} ${SAMPLE}.EPI.snp.nodensityfilter.final.vcf ${LINEAGESNPS})
	# echo -e ${SAMPLE}'\t'${lineage} > ${SAMPLE}.lineage.txt
	
	# remove no-intergenic
	#java -jar ${SNPEFF} -dataDir ${SNPEFFDATADIR} \
	#-v -nodownload -hgvs1LetterAa -noStats -no-downstream -no-intron -no-upstream -noLof \
	#mtb_h37rv ${SAMPLE}.EPI.snp.nodensityfilter.vcf 2>/dev/null | \
    #grep -o '^[^#]*' | \
    #sed 's/|/\t/g' > ${SAMPLE}.EPI.snp.nodensityfilter.snpeff.vcf

    #java -jar ${SNPEFF} -dataDir ${SNPEFFDATADIR} \
	#-v -nodownload -hgvs1LetterAa -noStats -no-downstream -no-intron -no-upstream -noLof \
	#mtb_h37rv ${SAMPLE}.EPI.snp.nodensityfilter.vcf > \
	#${SAMPLE}.EPI.snp.nodensityfilter.snpeff_nottabbed.vcf 2>/dev/null
done


#! /bin/bash
SAMPLENAMESFILE=$1
IFS=$'\n' read -d '' -r -a samples < ${SAMPLENAMESFILE}
[ -e vcf.stats ] && rm vcf.stats
for SAMPLE in "${samples[@]}"
do
	echo ${SAMPLE}
	[ -e ${SAMPLE}.vcf.stats ] && rm ${SAMPLE}.vcf.stats

	bcftools stats ${SAMPLE}.DR.snp.final.vcf 2>/dev/null | awk -v SAMPLE="$SAMPLE" '/^QUAL/ {print SAMPLE"\t""DR.snp.final.vcf""\t"$0}' >> vcf.stats

	bcftools stats ${SAMPLE}.EPI.snp.nodensityfilter.vcf 2>/dev/null | awk -v SAMPLE="$SAMPLE" '/^QUAL/ {print SAMPLE"\t""EPI.snp.nodensityfilter.vcf""\t"$0}' >> vcf.stats
	bcftools stats ${SAMPLE}.EPI.snp.nodensityfilter.excluderepeats.recode.vcf 2>/dev/null | awk -v SAMPLE="$SAMPLE" '/^QUAL/ {print SAMPLE"\t""EPI.snp.nodensityfilter.excluderepeats.recode.vcf""\t"$0}' >> vcf.stats
	bcftools stats ${SAMPLE}.EPI.snp.nodensityfilter.excluderepeats_drugres.recode.vcf 2>/dev/null | awk -v SAMPLE="$SAMPLE" '/^QUAL/ {print SAMPLE"\t""EPI.snp.nodensityfilter.excluderepeats_drugres.recode.vcf""\t"$0}' >> vcf.stats
	bcftools stats ${SAMPLE}.EPI.snp.nodensityfilter.excluderepeats_drugres_indels.recode.vcf 2>/dev/null | awk -v SAMPLE="$SAMPLE" '/^QUAL/ {print SAMPLE"\t""EPI.snp.nodensityfilter.excluderepeats_drugres_indels.recode.vcf""\t"$0}' >> vcf.stats
	bcftools stats ${SAMPLE}.EPI.snp.nodensityfilter.excluderepeats_indels.final.vcf.gz 2>/dev/null | awk -v SAMPLE="$SAMPLE" '/^QUAL/ {print SAMPLE"\t""EPI.snp.nodensityfilter.excluderepeats_indels.final.vcf.gz""\t"$0}' >> vcf.stats

	#bcftools stats ${SAMPLE}.gatk.indel.vcf 2>/dev/null |grep "^QUAL"  >> ${SAMPLE}.vcf.stats
	#bcftools stats ${SAMPLE}.gatk.indel.filter.vcf 2>/dev/null |grep "^QUAL"  >> ${SAMPLE}.vcf.stats
	#bcftools stats ${SAMPLE}.all.indel.vcf 2>/dev/null |grep "^QUAL"  >> ${SAMPLE}.vcf.stats
done

ls *EPI.snp.nodensityfilter.final.vcf.gz > EPI.snp.nodensityfilter.final.vcf.gz.files
bcftools merge --file-list EPI.snp.nodensityfilter.final.vcf.gz.files -O v --missing-to-ref -o EPI.snp.nodensityfilter.final.merged.vcf

${VCF2PHYLIP} \
--min-samples-locus 0 --write-used-sites \
-i EPI.snp.nodensityfilter.final.merged.vcf -f 

[ -e all.indel.bed ] && rm all.indel.bed
ls *all.indel.bed > all.indel.bed
${VCFTOOLS} --vcf EPI.snp.nodensityfilter.final.merged.vcf --exclude-bed all.indel.bed --recode --out final.merged.vcf

${VCF2PHYLIP} \
--min-samples-locus 0 --write-used-sites \
-i final.merged.vcf.recode.vcf -f


