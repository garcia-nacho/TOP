sample=test
index=../Hi17103bt2

r1=$(ls *R1*.fastq.gz)
r2=$(ls *R2*.fastq.gz)
(bowtie2 -p 1 -x ${index} -1 ${r1} -2 ${r2} -S ${sample}.sam) 2> ${sample}_Bowtie2summary.txt
samtools view -b -o ${sample}.bam  ${sample}.sam
samtools sort ${sample}.bam -o ${sample}.sorted.bam
samtools index ${sample}.sorted.bam
samtools depth -a ${sample}.sorted.bam > ${sample}_depth.tsv
rm *sam
