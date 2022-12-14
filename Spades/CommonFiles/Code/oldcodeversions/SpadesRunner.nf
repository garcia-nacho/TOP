// Declare syntax version
nextflow.enable.dsl=2

params.reads = "./RunX/*/*_{R1,R2}*.fastq.gz"
params.publishDir = './results'
params.threads = 10
params.toolcores = 2

process Trimming {

    maxForks = params.threads
    tag { sample }

    input:
    tuple val(sample), path(reads)
    
    output:
    path ("*P.fastq.gz")
     
    script:
    def (fq1, fq2) = reads

    """

    trimmomatic PE -basein ${fq1} -baseout ${sample}.fastq.gz  ILLUMINACLIP:/home/docker/CommonFiles/adapters/Kapa-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:3:15 MINLEN:36

    """
}

process TrimmingOutput {
    cpus = 1
    maxForks = 1

    publishDir(
        path: "${params.publishDir}/",
        mode: 'copy',
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    path(x)
    
    output:
    path ("*")

    script:

    """
    mkdir fastq
    cp *.fastq.gz fastq
    """

}


process Spades {
    cpus = params.toolcores
    maxForks = params.threads

    input:
    path(input)

    output:
    

    script:

    """
    """
}




process SpadesRun {
    cpus = params.toolcores
    maxForks = params.threads
    tag { sample }

   // publishDir(
   //     path: "${params.publishDir}/",
   //     mode: 'copy',
   //     saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
   // )

    input:
    tuple val(sample), path(reads) 

    
    output:

    path ("*{.tsv,.csv,.bam,.bai,Bowtie2summary.txt,contigs.fasta,.zip}")
    

    script:
    def (fq1, fq2) = reads

    """
    fastqc ${fq1}
    fastqc ${fq2}
    trimmomatic PE -basein ${fq1} -baseout ${sample}.fastq.gz  ILLUMINACLIP:/home/docker/CommonFiles/adapters/Kapa-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:3:15 MINLEN:36

    fastqc ${sample}_1P.fastq.gz 
    fastqc ${sample}_2P.fastq.gz
    spades.py -o . --careful --cov-cutoff auto -t ${params.toolcores} -1 ${sample}_1P.fastq.gz -2 ${sample}_2P.fastq.gz
    mv contigs.fasta ${sample}_contigs.fasta
    Rscript /home/docker/CommonFiles/Code/ContigCleaner.R
    mv clean_contigs.fasta ${sample}_clean_contigs.fasta
    mv clean_contigs.stats.csv ${sample}_contigs.stats.csv
    bowtie2-build ${sample}_clean_contigs.fasta ${sample}_bt2
    (bowtie2 -p ${params.toolcores} -x ${sample}_bt2 -1 ${sample}_1P.fastq.gz -2 ${sample}_2P.fastq.gz -S ${sample}.sam) 2> ${sample}_Bowtie2summary.txt
    samtools view -b -o ${sample}.bam  ${sample}.sam
    samtools sort ${sample}.bam -o ${sample}.sorted.bam
    samtools index ${sample}.sorted.bam
    
    rm ${sample}.sam
    rm ${sample}.bam
    

    /kraken2-2.1.2/kraken2 --use-names --report  ${sample}_R1_Raw_kraken_summaries.tsv --db /Kraken2DB/ ${fq1} > ${sample}_R1_Raw.kraken.tsv
    Rscript /home/docker/CommonFiles/Code/KrakenParser.R

    /kraken2-2.1.2/kraken2 --use-names --report  ${sample}_R2_Raw_kraken_summaries.tsv --db /Kraken2DB/ ${fq2} > ${sample}_R2_Raw.kraken.tsv
    Rscript /home/docker/CommonFiles/Code/KrakenParser.R

    /kraken2-2.1.2/kraken2 --use-names --report  ${sample}_R1_Trimmed_kraken_summaries.tsv --db /Kraken2DB/ ${sample}_1P.fastq.gz  > ${sample}_R1_Trimmed.kraken.tsv
    Rscript /home/docker/CommonFiles/Code/KrakenParser.R

    /kraken2-2.1.2/kraken2 --use-names --report  ${sample}_R2_Trimmed_summaries.tsv --db /Kraken2DB/ ${sample}_2P.fastq.gz  > ${sample}_R2_Trimmed.kraken.tsv
    Rscript /home/docker/CommonFiles/Code/KrakenParser.R

    /kraken2-2.1.2/kraken2 --use-names --report  ${sample}_Contigs_kraken_summaries.tsv --db /Kraken2DB/ ${sample}_clean_contigs.fasta  > ${sample}_cleancontigs.kraken.tsv

    rm *.fastq.gz
    Rscript /home/docker/CommonFiles/Code/KrakenParser.R

    """
}

process Integration {
    cpus 1
    maxForks = 1

    publishDir(
        path: "${params.publishDir}/",
        mode: 'copy',
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    path(x)
    
    output:
    path ("*")

    script:

    """
    multiqc ./
    Rscript /home/docker/CommonFiles/Code/Summarizer.R

    """

}
//New process to parse all summaries


workflow {
   sample_reads = Channel.fromFilePairs( params.reads )
   trimmingreads=Trimming(sample_reads)
   organizingtrim=TrimmingOutput(trimmingreads.collect())
   ouputspades=SpadesRun(sample_reads)
   outputintegration=Integration(ouputspades.collect())
}