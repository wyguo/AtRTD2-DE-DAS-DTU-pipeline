---
title: "BaRTv2.0: Short read assembly pipeline"
author: "Wenbin Guo"
date: "2021-06-28"
output: html_document
---


## Data pre-processing
### Raw RNA-seq read trimming

```shell
# Trimmomatic version 0.39 
trimmomatic PE $input_1 $input_2 $output_1P $output_1U $output_2P $output_2U ILLUMINACLIP:adapters.fasta:2:30:7 LEADING:20 TRAILING:20 MINLEN:20
```

### Fastqc read quality report

```shell
fastqc $fastq_file -o $report_dir
```

## RNA-seq read mapping
### STAR mapping index

```shell
### index star mapping pass1
STAR \
--runMode genomeGenerate \
--genomeDir $output_dir \
--genomeFastaFiles $genome_fasta \
--outFileNamePrefix $output_dir \
--limitGenomeGenerateRAM 240000000000

### merge sj files from star mapping pass1
cat star_result_pass1/*/SJ.out.tab | awk '($5 > 0 && $7 > 2 && $6==0)' | cut -f1-6 | sort | uniq > $sj_file

### index star mapping pass2
STAR \
--runMode genomeGenerate \
--genomeDir $output_dir \
--genomeFastaFiles $genome_fasta \
--outFileNamePrefix $output_dir \
--limitGenomeGenerateRAM 240000000000 \
--sjdbFileChrStartEnd $sj_file

```

### STAR mapping

```shell
### mapping scripts for pass1 and pass2
STAR \
--genomeDir $Index_dir \
--readFilesIn $read1 $read1 \
--sjdbOverhang 100 \
--alignIntronMin 60 \
--alignIntronMax 15000 \
--alignMatesGapMax 2000 \
--alignEndsType Local \
--alignSoftClipAtReferenceEnds Yes \
--outSAMprimaryFlag AllBestScore \
--outFilterType BySJout \
--outFilterMismatchNmax 0 \
--outFilterMismatchNoverLmax 0.3 \
--outFilterScoreMinOverLread 0.66 \
--outFilterMatchNmin 0 \
--outFilterScoreMin 0 \
--outFilterMultimapNmax 15 \
--outFilterIntronMotifs RemoveNoncanonical \
--outFilterIntronStrands RemoveInconsistentStrands \
--outSJfilterReads All \
--outSJfilterCountUniqueMin -1 5 5 5 \
--outSJfilterCountTotalMin -1 5 5 5 \
--outSAMstrandField intronMotif \
--outSAMtype BAM SortedByCoordinate \
--alignTranscriptsPerReadNmax 30000 \
--twopassMode None \
--readFilesCommand zcat \
--outReadsUnmapped Fastx \
--outFileNamePrefix $outfolder \
--outTmpDir $outTmpDir \
--alignSJoverhangMin 5 \
--alignSJDBoverhangMin 3 \
--outSJfilterOverhangMin -1 12 12 12 \
--outFilterMatchNminOverLread 0.66 \
--outFilterMismatchNoverReadLmax 1 \
--alignSJstitchMismatchNmax 0 0 0 0

```

## Transcript assemnly
### Cufflinks

```shell
cufflinks \
-u \
-F 0.01 \
--min-intron-length 60 \
--max-intron-length 15000 \
--library-type fr-firststrand \
-o $save_dir \
$bam_file
```

### Stringtie

```shell
stringtie $bam_file \
-o ${save_dir}trans.gtf \
--rf \
-A gene_abund.tab \
-a 5 \
-j 0.1 \
-f 0.3 \
-g 50 \
-M 1 \
-c 2.5
```

### Scallop

```shell
scallop \
-i $bam_file \
-o ${save_dir}Scallop_${sample_id}.gtf \
--library_type first
```
