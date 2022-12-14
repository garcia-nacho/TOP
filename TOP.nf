// Declare syntax version
nextflow.enable.dsl=2

params.reads = "./RunX/*/*_{R1,R2}*.fastq.gz"
params.publishDir = './results'
params.threads = 10
params.toolcores = 2

process Trimming {

    maxForks = params.threads - 1
    tag { sample }

    publishDir(
    path: "${params.publishDir}/fastq",
    mode: 'copy',
    saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    tuple val(sample), path(reads)
    
    output:
    tuple val(sample), path ("*1P.fastq.gz"), path("*2P.fastq.gz") 
     
    script:
    def (fq1, fq2) = reads

    """

    trimmomatic PE -basein ${fq1} -baseout ${sample}.fastq.gz  ILLUMINACLIP:/home/docker/CommonFiles/adapters/Kapa-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:3:15 MINLEN:36

    """
}

process KrakenRaw {
    maxForks = 1
    
    input:
    tuple val(sample), path(reads)

    output:
    path ("*{.tsv,.csv,.zip}")

    script:
    def (fq1, fq2) = reads

    """
    fastqc ${fq1}
    fastqc ${fq2}

    /kraken2-2.1.2/kraken2 --use-names --report  ${sample}_R1_Raw_kraken_summaries.tsv --db /Kraken2DB/ ${fq1} > ${sample}_R1_Raw.kraken.tsv
    Rscript /home/docker/CommonFiles/Code/KrakenParser.R

    /kraken2-2.1.2/kraken2 --use-names --report  ${sample}_R2_Raw_kraken_summaries.tsv --db /Kraken2DB/ ${fq2} > ${sample}_R2_Raw.kraken.tsv
    Rscript /home/docker/CommonFiles/Code/KrakenParser.R
    """ 

}


process Spades {

    maxForks = params.threads - 2

    input:
    tuple val(sample), path (trimmedR1), path(trimmedR2) 

    output:
    tuple val(sample), path ("*raw_contigs.fasta"),  path ("*clean_contigs.fasta"), path (trimmedR1), path(trimmedR2), emit: spadesraw 
    path("*_contigs.stats.csv"), emit: spadessum 
    path ("*clean_contigs.fasta"), emit: fastasclean
    path ("*raw_contigs.fasta"), emit: fastasraw
    val(sample), emit: sample_name

    script:

    """
    spades.py -o . --careful --cov-cutoff auto -t 1 -1 ${trimmedR1} -2 ${trimmedR2} 
    mv contigs.fasta ${sample}_raw_contigs.fasta
    Rscript /home/docker/CommonFiles/Code/ContigCleaner.R
    mv clean_contigs.fasta ${sample}_clean_contigs.fasta
    mv clean_contigs.stats.csv ${sample}_contigs.stats.csv
    
    """
}


process Rmlst {

    input:
    path(input)
    val(sample)

    output:
    path("*mlst{.json,.csv}"), emit: mlstfiles

    script:
    """
    (echo -n '{"base64":true,"details":true,"sequence": "'; base64 ${input}; echo '"}') | \
    curl -s -H "Content-Type: application/json" -X POST "http://rest.pubmlst.org/db/pubmlst_rmlst_seqdef_kiosk/schemes/1/sequence" -d @- > \
    ${sample}_rmlst.json
    Rscript /home/docker/CommonFiles/Code/rmlst_parser.R
    Rscript /home/docker/CommonFiles/Code/R_RunnerMLST.R

    """

}


process Prokka {
    container = 'garcianacho/prokka'

    input:
    path(input)
    val(sample)

    output:
    path("*")

    script:
    """
    prokka --outdir ${sample} --prefix ${sample} ${input} 
    """

}


process KrakenClean {
    maxForks = 2
    
    input:
    tuple val(sample), path(raw), path(clean) , path(trimmedR1), path(trimmedR2)

    output:
    path ("*{.tsv,.csv}")

    script:
  
    """
    /kraken2-2.1.2/kraken2 --use-names --report  ${sample}_Contigs_kraken_summaries.tsv --db /Kraken2DB/ ${clean} > ${sample}_cleancontigs.kraken.tsv
    Rscript /home/docker/CommonFiles/Code/KrakenParser.R
    """ 

}


process KrakenTrimmed {
    maxForks = 1
    
    input:
    tuple val(sample), path(trimmedR1), path(trimmedR2) 

    output:
    path ("*{.tsv,.csv,.zip}")

    script:
  

    """
    fastqc ${trimmedR1}
    fastqc ${trimmedR2}
    /kraken2-2.1.2/kraken2 --use-names --report  ${sample}_R1_Trimmed_kraken_summaries.tsv --db /Kraken2DB/ ${trimmedR1}  > ${sample}_R1_Trimmed.kraken.tsv
    Rscript /home/docker/CommonFiles/Code/KrakenParser.R

    /kraken2-2.1.2/kraken2 --use-names --report  ${sample}_R2_Trimmed_kraken_summaries.tsv --db /Kraken2DB/ ${trimmedR2}  > ${sample}_R2_Trimmed.kraken.tsv
    Rscript /home/docker/CommonFiles/Code/KrakenParser.R
    """ 

}

process Mapping {

    maxForks = 2
    
    input:
    tuple val(sample), path(raw),  path(clean), path(trimmedR1), path(trimmedR2) 

    output:
    path("*sorted.bam"), emit: bam
    path("*.bai"), emit: bai
    path ("*Bowtie2summary.txt"), emit: bt2sum

    script:

    """
    bowtie2-build ${clean} ${sample}_bt2
    (bowtie2 -p 1 -x ${sample}_bt2 -1 ${trimmedR1} -2 ${trimmedR2} -S ${sample}.sam) 2> ${sample}_Bowtie2summary.txt
    samtools view -b -o ${sample}.bam  ${sample}.sam
    samtools sort ${sample}.bam -o ${sample}.sorted.bam
    samtools index ${sample}.sorted.bam
    
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
    path(kkraw)
    path(kkclean)
    path(kktri)
    path(bam)
    path(bai)
    path(btsum)
    path(spsum)
    path(fastaclean)
    path(fastaraw)
    path(mlst_in)
  
    output:
    path ("*")

    script:

    """
    multiqc ./
    Rscript /home/docker/CommonFiles/Code/Summarizer.R
    mkdir bam
    mv *.bam ./bam
    mv *.bai ./bam
    mkdir fasta
    mv *_clean_contigs.fasta ./fasta
    mkdir QC
    mv *_raw_contigs.fasta ./QC
    mv *kraken* ./QC
    mv *contigs.stats.csv ./QC
    mv *.txt ./QC
    mv *.zip ./QC
    mv *.json ./QC
    mv *mlst.csv ./QC

    """
}


workflow {
   sample_reads = Channel.fromFilePairs( params.reads )
   trimmed=Trimming(sample_reads)
   ktrim=KrakenTrimmed(trimmed)
   kkraw=KrakenRaw(sample_reads)
   ouputspades=Spades(trimmed)
   kkcon=KrakenClean(ouputspades.spadesraw)
   mapped=Mapping(ouputspades.spadesraw)
   mlst=Rmlst(ouputspades.fastasclean, ouputspades.sample_name)
   integ=Integration(kkraw.collect(),
                     kkcon.collect(),
                     ktrim.collect(),
                     mapped.bam.collect(),
                     mapped.bai.collect(),
                     mapped.bt2sum.collect(),
                     ouputspades.spadessum.collect(),
                     ouputspades.fastasclean.collect(),
                     ouputspades.fastasraw.collect(),
                     mlst.mlstfiles.collect())
}