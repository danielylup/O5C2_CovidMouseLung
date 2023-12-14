#!/usr/bin/bash

INPUT_PATH = "~/Project/nx107"
SAMPLES = "GG01_S1 GG03_S3 GG05_S5 GG06_S6"

for SAMPLE in ${SAMPLES}; do

fastp -w 16 -i ${INPUT_PATH}/Sample_${SAMPLE}/LTT23_${SAMPLE}_L001_R1_001.fastq.gz \
      -o ${INPUT_PATH}/QC_file/LTT23_${SAMPLE}_L001_R1.filtered.fastq.gz \
      -h ${INPUT_PATH}/QC_file/QC_report/LTT23_${SAMPLE}.html \
      -j ${INPUT_PATH}/QC_file/QC_report/LTT23_${SAMPLE}.json

STAR --runThreadN 56 \
     --readFilesCommand zcat
     --genomeDir ~/Project/Reference/Mmus/Ref/STAR/Mmus_GRCm39/ \
     --outFileNamePrefix ${INPUT_PATH}/STAR_output/${SAMPLE}_STAR \
     --readFilesIn ${INPUT_PATH}/QC_file/LTT23_${SAMPLE}_L001_R1.filtered.fastq.gz

samtools view -@ 56 -b ${INPUT_PATH}/STAR_output/${SAMPLE}_STAR/${SAMPLE}_STARAligned.out.sam | samtools sort -@ 56 -o ${INPUT_PATH}/STAR_output/${SAMPLE}_STAR/${SAMPLE}_sorted.bam

qualimap rnaseq -bam ${INPUT_PATH}/STAR_output/${SAMPLE}_STAR/LTT23_GG01_sorted.bam \
                -gtf ~/Project/Reference/Mmus/data/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.gtf \
                -outdir ${INPUT_PATH}/QC_file/qualimap_output/${SAMPLE}/ \
                -s --java-mem-size=56G

done

featureCounts -T 56 
              -a ~/Project/Reference/Mmus/data/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.gtf 
              -o ~/Project/nx107/featureCounts_output/bulkCML_GG_count_C2.tsv \
		  ~/Project/nx107/STAR_output/LTT23_GG01_STARAligned.out.sam \
		  ~/Project/nx107/STAR_output/LTT23_GG03_STARAligned.out.sam \
		  ~/Project/nx107/STAR_output/LTT23_GG05_STARAligned.out.sam \
		  ~/Project/nx107/STAR_output/LTT23_GG06_STARAligned.out.sam
