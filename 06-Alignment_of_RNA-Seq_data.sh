# Quality and adapter trimming:
for f1 in *R1.fq.gz; do f2=${f1/_R1/_R2} && echo $f1 $f2 && trim_galore --paired --illumina $f1 $f2; done

for file in *val*fq.gz; do mv $file ${file/val_?/trimmed}; done

# Alignment to TAIR10:
for f1 in *R1_trimmed.fq.gz; do f2=${f1/_R1/_R2} && echo $f1 $f2 && STAR --genomeDir tair10 --readFilesIn $f1 $f2 --readFilesCommand zcat --runThreadN 4 --outFileNamePrefix ${f1/R1_trimmed.fq.gz/} --outSAMmultNmax 1 --alignEndsType Local; done
rm *out *tab; rmdir *STARtmp

for file in *sam; do mv $file ${file/_Aligned.out/}; done

# Sort SAM files, filter for MAPQ >= 10 and convert to BAM
# (these BAM files are used as input for calling DE genes and introns in 07-Differential_expression.R):
for file in *sam; do echo $file && samtools view -huq 10 $file | samtools sort - -o ${file/.sam/_sorted.bam} && rm $file; done

# Merge replicates:
for f1 in *biorep1*sorted.bam; do f2=${file/rep1/rep2} && echo $f1 $f2 && samtools merge ${f1/biorep1/merged/} $f1 $f2; done

# Make unstranded Bedgraph files (for visualization in genomic browsers):
for file in *bam; do echo $file && bedtools genomecov -ibam $file -bg -split | sort -k1,1 -k2,2n | sed "1i track type=bedGraph" | gzip > ${file/bam/bedgraph.gz}; done

