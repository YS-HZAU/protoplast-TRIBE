#!/bin/sh
#########################################################
#    This BASH script was used to call SNPs using 
#         Cutadapt(v4.9)
#         STAR(v2.5.3a);
#         SAMtools(v1.7-1);
#         Picard(v2.23.1);
#         bcftools(v1.7);
#    Please free to use other methods and softwares 
#    for calling SNPs as substitution
#########################################################
# ---------Update the following varibles------------
ref_genome="/home/ys/Desktop/japonica_genome/ensemble_plant_genome/ensemble_plant_annotation_file/Oryza_sativa.IRGSP_1.0.dna.toplevel.fa"
gtf="/home/ys/Desktop/japonica_genome/ensemble_plant_genome/ensemble_plant_annotation_file/Oryza_sativa.IRGSP_1.0.49.gtf"
star_indices="/home/ys/Desktop/japonica_genome/ensemble_plant_genome/rice_ensemble_plant_STAR_index"
TRIMMOMATIC_JAR="/home/ys/biosoft/trimmomatic/Trimmomatic-0.38/trimmomatic-0.38.jar"
PICARD_JAR="/home/ys/biosoft/picard/picard_2.23.1/picard.jar"
# ---------End Update variable------------

######################################################################
# Filter raw data
# all paired-end raw fq.gz files for TRIBEseq should in this directory
######################################################################

# Use the loop below to run the code for each paired-end sample (*.1.fq.gz, *.2.fq.gz),
# filelist=`ls .1.fq.gz`#
# for file in ${filelist[@]} 
# do
sample=$1
index=$(basename $sample |sed 's/1.fq.gz//')
trim_galore -j 4 -q 30 --phred33 --stringency 3 --length 150 -e 0.1 --paired --fastqc ${index}1.fq.gz ${index}2.fq.gz --gzip -o ./trim_clean_data
# done

## cd ./trim_clean_data

# Now use the clean data for aligning library with STAR
# # Use the loop below to run the code
# for file in $(ls *_1_val_1.fq.gz)
# do
sample=$1
prefix=$(basename $sample |sed 's/_1_val_1.fq.gz//')
echo
echo ==============================
echo  "Now using" $sample "for alignment"
echo ==============================
echo

STAR  --runThreadN 40 --outFilterMismatchNoverLmax 0.07 --outFileNamePrefix ./trim_clean_data/star_alignment/"$prefix"_ --outSAMstrandField intronMotif --outFilterMultimapNmax 1  --genomeDir $star_indices --readFilesIn  ./trim_clean_data/"$prefix"_1_val_1.fq.gz ./trim_clean_data/"$prefix"_2_val_2.fq.gz --sjdbGTFfile $gtf --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat

## cd ./trim_clean_data/star_alignment

# following code is tested with samtools 1.7.1, you might have to tweak it a bit bases your installed verison of samtools (these flags can be problematic for older version of samtools: -@, -o)
# remove low quality alignment and convert sam to bam
samtools view -@ 10 -bSh -q 30 $prefix"_Aligned.sortedByCoord.out.bam" > $prefix"_Aligned.sortedByCoord.out.HQ30.bam"

# summary the main aligement information
if [ ! -d align_summary ]; then
	mkdir align_summary
else
	echo "============================"
	echo "align_summary exists"
	echo "============================"
fi

grep -E "Number of input reads|Average input read length|Uniquely mapped reads number|Uniquely mapped reads %|Average mapped length" $prefix"_Log.final.out" > $prefix"_align_summary.txt"
sed -i "1 i $prefix"_ailgn_infomation"\tresults" $prefix"_align_summary.txt"

# sorted bam file is used for assessing gene expression through featurecounts
count_out=$prefix"_featurecount.txt"
samtools index $prefix"_Aligned.sortedByCoord.out.HQ30.bam"
featureCounts -T 8 -p -t exon -g gene_id -a $gtf -o ./$count_out $prefix"_Aligned.sortedByCoord.out.HQ30.bam"

# Use the loop to run Picard to remove duplicates
sample=$prefix"_Aligned.sortedByCoord.out.HQ30.bam"
index=$(basename $sample |sed 's/.bam//')
picard MarkDuplicates INPUT=$sample OUTPUT=${index}.nodup.bam METRICS_FILE= ${index}.dup.txt VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true TMP_DIR=tmp ASSUME_SORTED=true

# Use the loop to extract vcf file from nodup.bam file
# for file in $(ls *.nodup.bam)
# do 
sample=$1
index=$(basename $sample |sed 's/.bam//')
samtools mpileup -ugf $ref_genome $sample > ./${index}.bcf &&  bcftools call -vm -Oz -o ./${index}_against_REF.vcf.gz ./${index}.bcf"

# Use the loop to summary the vcf file
# for file in $(ls *_against_REF.vcf.gz)
# do
sample=$1
bcftools stats $sample > ${sample}.stats && cat ${sample}.stats grep "^ST" ${sample}.stats | awk '{print $3"\t"$4}'| sed "1 i type_of_mutation\t${sample%*stats}" > ${sample}.stats"_temp"
# After all vcf file were get stats-ed, then run follow code
paste *_temp > nucleotide_mutation_summary.txt

