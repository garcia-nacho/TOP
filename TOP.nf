#!/usr/bin/env nextflow

// Declare syntax version
nextflow.enable.dsl=2

params.readsfolder = "."
params.publishDir = params.readsfolder+"/TOPresults"
params.threads = 10
params.toolcores = 2
params.forceSp="none"

params.reads=params.readsfolder+"/*/*_{R1,R2}*.fastq.gz"

process Trimming {
 
    container 'garcianacho/top:spades'
    //containerOptions '--volume /media/nacho/Data/kraken2_standard_20220926/:/Kraken2DB'

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
    container 'garcianacho/top:spades'
    containerOptions '--volume /media/nacho/Data/kraken2_standard_20220926/:/Kraken2DB'
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

    container 'garcianacho/top:spades'
    //containerOptions '--volume /media/nacho/Data/kraken2_standard_20220926/:/Kraken2DB'

    maxForks = params.threads - 2

    input:
    tuple val(sample), path (trimmedR1), path(trimmedR2) 

    output:
    tuple val(sample), path ("*raw_contigs.fasta"),  path ("*clean_contigs.fasta"), path (trimmedR1), path(trimmedR2), emit: spadesraw 
    path("*_contigs.stats.csv"), emit: spadessum 
    path ("*clean_contigs.fasta"), emit: fastasclean
    path ("*raw_contigs.fasta"), emit: fastasraw
    path (trimmedR1), emit: r1spades
    path (trimmedR2), emit: r2spades

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
    container 'garcianacho/top:spades'
    //containerOptions '--volume /media/nacho/Data/kraken2_standard_20220926/:/Kraken2DB'
    input:
    path(input)
    val(sample)
    path(r1)
    path(r2)

    output:
    path("*mlst{.json,.csv}"), emit: mlstresults
    path("*.agent"), emit: agent
    path(input), emit: clean_contigs_frommlst
    val(sample), emit: sample_frommlst
    path (r1), emit: r1mlst
    path (r2), emit: r2mlst

    shell:
    """
    (echo -n '{"base64":true,"details":true,"sequence": "'; base64 ${input}; echo '"}') | \
    curl -s -H "Content-Type: application/json" -X POST "http://rest.pubmlst.org/db/pubmlst_rmlst_seqdef_kiosk/schemes/1/sequence" -d @- > \
    ${sample}_rmlst.json
    Rscript /home/docker/CommonFiles/Code/rmlst_parser.R
    Rscript /home/docker/CommonFiles/Code/seqmlst_parser.R
    #Missing genus 

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
    container 'garcianacho/top:spades'
    containerOptions '--volume /media/nacho/Data/kraken2_standard_20220926/:/Kraken2DB'
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
    container 'garcianacho/top:spades'
    containerOptions '--volume /media/nacho/Data/kraken2_standard_20220926/:/Kraken2DB'
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
    container 'garcianacho/top:spades'
    containerOptions '--volume /media/nacho/Data/kraken2_standard_20220926/:/Kraken2DB'

    maxForks = 2
    
    input:
    tuple val(sample), path(raw),  path(clean), path(trimmedR1), path(trimmedR2) 

    output:
    path("*sorted.bam"), emit: bam
    path("*.bai"), emit: bai
    path ("*Bowtie2summary.txt"), emit: bt2sum
    path ("*depth.tsv"), emit: bt2depth

    script:

    """
    bowtie2-build ${clean} ${sample}_bt2
    (bowtie2 -p 1 -x ${sample}_bt2 -1 ${trimmedR1} -2 ${trimmedR2} -S ${sample}.sam) 2> ${sample}_Bowtie2summary.txt
    samtools view -b -o ${sample}.bam  ${sample}.sam
    samtools sort ${sample}.bam -o ${sample}.sorted.bam
    samtools index ${sample}.sorted.bam
    samtools depth -a ${sample}.sorted.bam > ${sample}_depth.tsv
    
    """ 

}

process Integration {
    container 'garcianacho/top:spades'
    //containerOptions '--volume /media/nacho/Data/kraken2_standard_20220926/:/Kraken2DB'
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
    path(abri)
    path(hicap)
    path(seroba)
    path(depth)
  
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
  
    mv *_Abricate.csv ./QC

    if test -f "*.tsv"; then mv *.tsv ./QC; fi
    if test -f "*.gbk"; then mv *.gbk ./QC; fi
    if test -f "*.log"; then mv *.log ./QC; fi
    if test -f "*.svg"; then mv *.svg ./QC; fi
    """
}

process Abricate { 
    container 'garcianacho/top:abricate'
    cpus 1
    maxForks = 1

    input:
    path(fastaclean)
    val(sample)
    path(agent)

    output:
    path("*{.tsv,.csv}"), emit: abricate_results

    script:

    """
    if test -f "Hinf.agent"; 
    then
        abricate-get_db --db ncbi --force
        abricate-get_db --db vfdb --force
        abricate --db vfdb --quiet *.fasta > ${sample}_vfdb.tsv 
        abricate --db HinfFtsI --quiet *.fasta > ${sample}_HinfFtsI.tsv
        abricate --db HinfGyrSubA --quiet *.fasta > ${sample}_HinfGyrSubA.tsv
        abricate --db HinfTopoIVsubA --quiet *.fasta > ${sample}_HinfTopoIVSubA.tsv
        abricate --db ncbi --quiet *.fasta > ${sample}_ncbi.tsv
        #Integration Abricate
    else
        abricate-get_db --db ncbi --force
        abricate-get_db --db vfdb --force
        abricate --db vfdb --quiet *.fasta > ${sample}_vfdb.tsv
        abricate --db ncbi --quiet *.fasta > ${sample}_ncbi.tsv
        #Dummy file
    fi

    #Test if the files can be empty or they just dont exist
    
    Rscript /home/docker/Code/AbricateParser.R
    mv Abricate.csv ${sample}_Abricate.csv

    """

}


process Hicap { 
    container 'garcianacho/top:hicap'
    cpus 1
    maxForks = 1

    input:
    path(fastaclean)
    val(sample)
    path(agent)

    output:
    path("*"), emit: hicap_results

    script:

    """
    if test -f "Hinf.agent"; 
    then
        /home/docker/Code/hicapwrapper.sh

    else
        cat "NoHi" > ${sample}_HiCap.tsv
        #Dummy Hicap file

    fi

    """

}

process Seroba { 
    container 'garcianacho/top:seroba'
    cpus 1
    maxForks = 1

    input:
    path(r1)
    path(r2)
    val(sample)
    path(agent)

    output:
    path("*.tsv"), emit: seroba_results

    script:

    """
    source activate seroba
    
    seroba runSerotyping /home/docker/seroba/database/ ${r1} ${r2} ${sample}
    mv ${sample}/*.tsv ./${sample}_seroba.tsv
    rm -rf dummy

    conda deactivate

    """

}


workflow {
   sample_reads = Channel.fromFilePairs( params.reads )
   trimmed=Trimming(sample_reads)
   ktrim=KrakenTrimmed(trimmed)
   kkraw=KrakenRaw(sample_reads)
   outputspades=Spades(trimmed)
   kkcon=KrakenClean(outputspades.spadesraw)
   mapped=Mapping(outputspades.spadesraw)
   mlst=Rmlst(outputspades.fastasclean, outputspades.sample_name,outputspades.r1spades, outputspades.r2spades)
   abri=Abricate(mlst.clean_contigs_frommlst, mlst.sample_frommlst, mlst.agent)
   hicap=Hicap(mlst.clean_contigs_frommlst, mlst.sample_frommlst, mlst.agent)
   seroba=Seroba(mlst.r1mlst, mlst.r2mlst, mlst.sample_frommlst, mlst.agent)
   integ=Integration(kkraw.collect(),
                     kkcon.collect(),
                     ktrim.collect(),
                     mapped.bam.collect(),
                     mapped.bai.collect(),
                     mapped.bt2sum.collect(),
                     outputspades.spadessum.collect(),
                     outputspades.fastasclean.collect(),
                     outputspades.fastasraw.collect(),
                     mlst.mlstresults.collect(),
                     abri.abricate_results.collect(),
                     hicap.hicap_results.collect(),
                     seroba.seroba_results.collect()
                     mapped.bt2depth.collect())
}