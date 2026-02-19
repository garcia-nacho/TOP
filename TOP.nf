#!/usr/bin/env nextflow

// Declare syntax version
nextflow.enable.dsl=2
println ">>> softtrimming param received: ${params.softtrimming}"

params.readsfolder = "."
params.publishDir = params.readsfolder+"/TOPresults"
params.threads = 10
params.spadescores = 8
params.forceSp="none"
params.krakenDB="/media/nacho/Data/kraken2_standard_20220926/"
params.TBDB="/mnt/N/NGS/TB_pipeline/TB_pipeline_database/DB/"
params.tempfolder="/media/nacho/Data/temp/toptest/tempdb/"
params.devrun="No"
params.min_reads = 40000 
params.adapters="TruSeq"
params.progress_dir = "/media/nacho/Data/temp/toptest/TOPv1.1test/test"
params.softtrimming = "No"

params.reads=params.readsfolder+"/*/*_{R1,R2}*.fastq.gz"

process PreFilter {

    container 'ghcr.io/garcia-nacho/top_spades:v1.1'
    containerOptions '--volume '+params.progress_dir+':/logs'
    maxForks = params.threads - 1
    errorStrategy 'ignore'
    tag { sample }

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path(reads)

    script:
    """
    r1=\$(ls ${sample}_R1*.fastq.gz)
    r2=\$(ls ${sample}_R2*.fastq.gz)

    num_reads_r1=\$(zcat \$r1 | wc -l)
    num_reads_r1=\$((num_reads_r1 / 4))

    num_reads_r2=\$(zcat \$r2 | wc -l)
    num_reads_r2=\$((num_reads_r2 / 4))
    echo $task.attempt > attemp.txt
    if [ \$num_reads_r1 -ge ${params.min_reads} ] && [ \$num_reads_r2 -ge ${params.min_reads} ]; then
        echo "$sample passed filtering with \$num_reads_r1 and \$num_reads_r2 reads."
        echo "$sample" > pass_filter.txt
        

        logstring="${sample},Prefilter,Completed"
        logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
        echo \${logstring} >> /logs/\${logkey}_logs.tsv

        
    else


        logstring="${sample},Prefilter,Failed"
        logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
        echo \${logstring} >> /logs/\${logkey}_logs.tsv
        exit 1

    fi
    """
}

process Trimming {
 
    container 'ghcr.io/garcia-nacho/top_spades:v1.1'
    containerOptions '--volume '+params.progress_dir+':/logs'

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
    ln -s *_R1* ${sample}_ln_R1_001.fastq.gz
    ln -s *_R2* ${sample}_ln_R2_001.fastq.gz
    echo ">>> adapters = ${params.adapters}"
    echo ">>> softtrimming = ${params.softtrimming}"
    
    if [ "${params.adapters}" == "Kapa" ]
    then
        if [ "${params.softtrimming}" == "Yes" ]
        then
            trimmomatic PE -phred33 -basein ${sample}_ln_R1_001.fastq.gz -baseout ${sample}.fastq.gz  ILLUMINACLIP:/home/docker/CommonFiles/adapters/Kapa-PE.fa:2:30:10
        else
            trimmomatic PE -phred33 -basein ${sample}_ln_R1_001.fastq.gz -baseout ${sample}.fastq.gz  ILLUMINACLIP:/home/docker/CommonFiles/adapters/Kapa-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:3:15 MINLEN:36
        fi
        
    else
        if [ "${params.softtrimming}" == "Yes" ]
        then
            trimmomatic PE -phred33 -basein ${sample}_ln_R1_001.fastq.gz -baseout ${sample}.fastq.gz  ILLUMINACLIP:/home/docker/CommonFiles/adapters/TruSeq3-PE-2.fa:2:30:10
        else
            trimmomatic PE -phred33 -basein ${sample}_ln_R1_001.fastq.gz -baseout ${sample}.fastq.gz  ILLUMINACLIP:/home/docker/CommonFiles/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:3:15 MINLEN:36
        fi
    fi

    rm ${sample}_ln_R1_001.fastq.gz
    rm ${sample}_ln_R2_001.fastq.gz

    logstring="${sample},Trimming,Completed"
    logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
    echo \${logstring} >> /logs/\${logkey}_logs.tsv


    """
}

process KrakenRaw {
    container 'ghcr.io/garcia-nacho/top_spades:v1.1'
    containerOptions '--volume '+params.krakenDB+':/Kraken2DB --volume '+params.progress_dir+':/logs'
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

    # /kraken2-2.1.2/kraken2 --use-names --memory-mapping --report  ${sample}_R1_Raw_kraken_summaries.tsv --db /Kraken2DB/ ${fq1} > ${sample}_R1_Raw.kraken.tsv
    # Rscript /home/docker/CommonFiles/Code/KrakenParser.R

    /kraken2-2.1.2/kraken2 --use-names --memory-mapping --report  ${sample}_R2_Raw_kraken_summaries.tsv --db /Kraken2DB/ --paired ${fq1} ${fq2} > ${sample}_Raw.kraken.tsv
    Rscript /home/docker/CommonFiles/Code/KrakenParser.R



    logstring="${sample},KrakenRaw,Completed"
    logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
    echo \${logstring} >> /logs/\${logkey}_logs.tsv

    """ 

}


process CPUcounter {

    input:
    tuple val(sample), path (trimmedR1), path(trimmedR2) 

    output:
    tuple val(sample), path (trimmedR1), path(trimmedR2), path ("cpus.txt"), emit: cpu_counts   

    script: 
    """
    cpu_spades=\$(ps aux | grep '[s]pades.py' | wc -l)
    cpu_trimmomatic=\$(ps aux | grep '[t]rimmomatic' | wc -l)
    threads_spades=\$(ps -eo cmd | grep '[s]pades.py' | grep -Eo '(-t|--threads) [0-9]+' | awk '{sum += \$2} END {print sum}')

    echo "\$cpu_spades \$cpu_trimmomatic" > cpus.txt 
    """
 
}

process Spades {

    container 'ghcr.io/garcia-nacho/top_spades:v1.1'
    containerOptions '--volume '+params.progress_dir+':/logs'

    maxForks = params.spadescores

    input:
    tuple val(sample), path (trimmedR1), path(trimmedR2), path(cpu_file) 

    output:
    tuple val(sample), path ("*raw_contigs.fasta"),  path ("*clean_contigs.fasta"), path (trimmedR1), path(trimmedR2), emit: spadesraw 
    path("*_contigs.stats.csv"), emit: spadessum 
    path ("*clean_contigs.fasta"), emit: fastasclean
    path ("*raw_contigs.fasta"), emit: fastasraw
    path (trimmedR1), emit: r1spades
    path (trimmedR2), emit: r2spades
    path ("*sequencerID.tsv"), emit: sqID

    val(sample), emit: sample_name

    script:

    """
    read cpu_spades cpu_trimmomatic < ${cpu_file}

    total_used=\$(( cpu_spades  +  cpu_trimmomatic +1 )) 
    

    threadssp=\$(echo "${params.spadescores} / \$total_used" | bc)

 
    if [ \$threadssp -lt 1 ]; then
        threadssp=1
    fi



    sqid=\$(gzip -cd ${trimmedR1} | head -n 1)   
    echo \${sqid}  > ${sample}_sequencerID.tsv
    spades.py -o . --careful --cov-cutoff auto -t \${threadssp} -1 ${trimmedR1} -2 ${trimmedR2} 
    mv contigs.fasta ${sample}_raw_contigs.fasta
    Rscript /home/docker/CommonFiles/Code/ContigCleaner.R
    mv clean_contigs.fasta ${sample}_clean_contigs.fasta
    source activate coverm
    coverm genome -1 ${trimmedR1} -2 ${trimmedR2} -r ${sample}_clean_contigs.fasta -t \${threadssp} --single-genome -m mean > read_coverage.tsv
    conda deactivate
    Rscript /home/docker/CommonFiles/Code/readcovadd.R

    mv clean_contigs.stats.csv ${sample}_contigs.stats.csv
    rm -rf ./corrected
    rm read_coverage.tsv

    logstring="${sample},Assembly,Completed"
    logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
    echo \${logstring} >> /logs/\${logkey}_logs.tsv

    
    """
}

process Rmlst {
    container 'ghcr.io/garcia-nacho/top_spades:v1.1'
    containerOptions '--volume '+params.progress_dir+':/logs'
    errorStrategy 'retry'
    //errorStrategy { task.attempt < 10 ? 'retry' : 'ignore' }
    maxRetries 10

    
    maxForks = 1
    

    input:
    path(input)
    val(sample)
    path(r1)
    path(r2)

    output:
    path("*mlst{.json,.csv}"), emit: mlstresults
    path("*.agent"), emit: agent
    path("*localmlst.tsv"), emit: localmlst
    path(input), emit: clean_contigs_frommlst
    val(sample), emit: sample_frommlst
    path (r1), emit: r1mlst
    path (r2), emit: r2mlst

    shell:
    """
    (echo -n '{"base64":true,"details":true,"sequence": "'; base64 ${input}; echo '"}') | \
    curl -s -H "Content-Type: application/json" -X POST "https://rest.pubmlst.org/db/pubmlst_rmlst_seqdef_kiosk/schemes/1/sequence" -d @- > \
    ${sample}_rmlst.json
    Rscript /home/docker/CommonFiles/Code/rmlst_parser.R
    Rscript /home/docker/CommonFiles/Code/seqmlst_parser.R
    #Missing genus 
    #Produce error if nothing happened
    source activate mlst
    mlst --blastdb /home/docker/CommonFiles/blast/mlst.fa ${sample}_clean_contigs.fasta > ${sample}_localmlst.tsv
    conda deactivate

    
    if [[ "$task.attempt" -lt 9 ]]; then

    Rscript /home/docker/CommonFiles/Code/ErrorRaiserMLST.R
    
    fi

    logstring="${sample},MLST,Completed"
    logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
    echo \${logstring} >> /logs/\${logkey}_logs.tsv

    """

}

process KrakenClean {
    container 'ghcr.io/garcia-nacho/top_spades:v1.1'
    containerOptions '--volume '+params.krakenDB+':/Kraken2DB --volume '+params.progress_dir+':/logs'
    maxForks = 1
    
    input:
    tuple val(sample), path(raw), path(clean) , path(trimmedR1), path(trimmedR2)

    output:
    path ("*{.tsv,.csv}")

    script:
  
    """
    /kraken2-2.1.2/kraken2 --use-names --memory-mapping --report  ${sample}_Contigs_kraken_summaries.tsv --db /Kraken2DB/ ${clean} > ${sample}_cleancontigs.kraken.tsv
    Rscript /home/docker/CommonFiles/Code/KrakenParser.R
    echo "${sample},KrakenClean,Completed" >> /logs/run_logs.tsv
    """ 

}

process KrakenTrimmed {
    container 'ghcr.io/garcia-nacho/top_spades:v1.1'
    containerOptions '--volume '+params.krakenDB+':/Kraken2DB'
    maxForks = 1
    
    input:
    tuple val(sample), path(trimmedR1), path(trimmedR2) 

    output:
    path ("*{.tsv,.csv,.zip}")

    script:
  

    """
    fastqc ${trimmedR1}
    fastqc ${trimmedR2}
    # /kraken2-2.1.2/kraken2 --use-names --memory-mapping --report  ${sample}_R1_Trimmed_kraken_summaries.tsv --db /Kraken2DB/ ${trimmedR1}  > ${sample}_R1_Trimmed.kraken.tsv
    # Rscript /home/docker/CommonFiles/Code/KrakenParser.R

    /kraken2-2.1.2/kraken2 --use-names --memory-mapping --report  ${sample}_R2_Trimmed_kraken_summaries.tsv --db /Kraken2DB/ --paired ${trimmedR1} ${trimmedR2}  > ${sample}_Trimmed.kraken.tsv
    Rscript /home/docker/CommonFiles/Code/KrakenParser.R
    """ 

}

process CPUcounterBT2 {

    input:
    tuple val(sample), path(raw),  path(clean), path(trimmedR1), path(trimmedR2) 

    output:
    tuple val(sample), path(raw),  path(clean), path(trimmedR1), path(trimmedR2) , path ("cpus_bt2.txt"), emit: cpu_bt2  

    script: 
    """
    cpu_bt2=\$(ps aux | grep '[s]bowtie' | wc -l) 
    cpu_spades=\$(ps aux | grep '[s]pades.py' | wc -l)
    echo "\$cpu_bt2 \$cpu_bt2" > cpus_bt2.txt 
    """
 
}

process Mapping {
    container 'ghcr.io/garcia-nacho/top_spades:v1.1'
    containerOptions '--volume '+params.progress_dir+':/logs'

    maxForks = 1
    
    input:
    tuple val(sample), path(raw),  path(clean), path(trimmedR1), path(trimmedR2),path (cpu_file2) 
    

    output:
    path("*sorted.bam"), emit: bam
    path("*.bai"), emit: bai
    path ("*Bowtie2summary.txt"), emit: bt2sum
    path ("*depth.tsv"), emit: bt2depth

    script:

    """
    read cpu_bt2 cpu_spades < ${cpu_file2}

    total_used=\$(( cpu_spades  +  cpu_bt2 +1 )) 

    threadssp=1

    if [ \$total_used -lt 3 ]; then
        threadssp=${params.spadescores}
    fi



    bowtie2-build ${clean} ${sample}_bt2
    (bowtie2 -p \$threadssp -x ${sample}_bt2 -1 ${trimmedR1} -2 ${trimmedR2} -S ${sample}.sam) 2> ${sample}_Bowtie2summary.txt
    samtools view -b -o ${sample}.bam  ${sample}.sam
    samtools sort ${sample}.bam -o ${sample}.sorted.bam
    samtools index ${sample}.sorted.bam
    samtools depth -a ${sample}.sorted.bam > ${sample}_depth.tsv
    rm ${sample}.sam
    rm ${sample}.bam

    logstring="${sample},Mapping,Completed"
    logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
    echo \${logstring} >> /logs/\${logkey}_logs.tsv


    """ 

}

process Abricate { 
    container 'ghcr.io/garcia-nacho/top_abricate:v1.1'
    containerOptions '--volume '+params.progress_dir+':/logs'
    //cpus 1
    //maxForks = 1

    input:
    path(fastaclean)
    val(sample)
    path(agent)

    output:
    path("*{.tsv,.csv}"), emit: abricate_results

    script:

    """
    abricate --db vfdb2 --quiet *.fasta > ${sample}_vfdb.tsv
    abricate --db ncbi --quiet *.fasta > ${sample}_ncbi.tsv
    abricate --db plasmidfinder --quiet *.fasta > ${sample}_plasmidfinder.tsv
    Rscript /home/docker/Code/AbricateParser.R
    mv Abricate.csv ${sample}_Abricate.csv

    logstring="${sample},Abricate,Completed"
    logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
    echo \${logstring} >> /logs/\${logkey}_logs.tsv


    """

}

process NGstar { 
    container 'ghcr.io/garcia-nacho/top_ngstar:v1.1'
    containerOptions '--volume '+params.progress_dir+':/logs'
    cpus 1
    maxForks = 1

    input:
    path(fastaclean)
    val(sample)
    path(agent)

    output:
    path("*ngstar_results.tsv"), emit: ngstar_results

    script:

    """
    if test -f "Ngon.agent"; 
    then
      cp -L ${sample}_clean_contigs.fasta ${sample}_clean_contigs_nolink.fasta
      cat ${sample}_clean_contigs_nolink.fasta | tr 'a-z' 'A-Z' > ${sample}_clean_contigs_nocap.fasta
      #python3 /home/docker/pyngSTar/pyngSTar.py -f -i ${sample}_clean_contigs_nocap.fasta -p /home/docker/pyngSTar/pyngSTarDB_02012024/ -o ${sample}_ngstar_results.tsv
      
      #URL="https://ngstar.canada.ca/alleles/download?lang=en&loci_name=23S"
      
      mkdir /home/docker/pyngSTar/pyngSTarDB_rolling/
      Rscript /home/docker/CreateDB.R
      rm /home/docker/pyngSTar/pyngSTarDB_rolling/profiles.xlsx
      python3 /home/docker/pyngSTar/pyngSTar.py -p /home/docker/pyngSTar/pyngSTarDB_rolling/ -u
       
      #python3 /home/docker/pyngSTar/pyngSTar.py -f -i ${sample}_clean_contigs_nocap.fasta -p /home/docker/pyngSTar/pyngSTarDB_03032025/ -o ${sample}_ngstar_results.tsv
      python3 /home/docker/pyngSTar/pyngSTar.py -f -i ${sample}_clean_contigs_nocap.fasta -p /home/docker/pyngSTar/pyngSTarDB_rolling/ -o ${sample}_ngstar_results.tsv  
      

        logstring="${sample},NGStar,Completed"
        logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
        echo \${logstring} >> /logs/\${logkey}_logs.tsv
      
    else
      echo "NoNgon" > ${sample}_ngstar_results.tsv

      logstring="${sample},NGStar,Skipped"
      logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
      echo \${logstring} >> /logs/\${logkey}_logs.tsv

    fi

    """

}


process HiCgmlst{
    container 'ghcr.io/garcia-nacho/top_spades:v1.1'
    containerOptions '--volume '+params.progress_dir+':/logs'

    input:
    path(fastaclean)
    val(sample)
    path(agent)

    output:
    path("*Hicgmlst.csv"), emit: hicgmlst_results

    
    script:

    """
    if [[ -f "Hinf.agent" ]]
    then

        cp /home/docker/CommonFiles/profiles_cgmlst.tsv.tar.gz ./profiles_cgmlst.tsv.tar.gz && tar -xf profiles_cgmlst.tsv.tar.gz
        #GET https://rest.pubmlst.org/db/pubmlst_hinfluenzae_seqdef/schemes/56/profiles_csv > profiles_cgmlst.tsv
    
        /home/docker/CommonFiles/Code/REST_Runner.sh ${sample}_clean_contigs.fasta https://rest.pubmlst.org/db/pubmlst_hinfluenzae_seqdef/schemes/56/sequence ${sample}_Hicgmlst.json
        Rscript /home/docker/CommonFiles/Code/Hicgmlst.R
        mv Hicgmlst.csv ${sample}_Hicgmlst.csv

        logstring="${sample},Hicgmlst,Completed"
        logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
        echo \${logstring} >> /logs/\${logkey}_logs.tsv

    else
        echo "NoHi" > ${sample}_Hicgmlst.csv
        logstring="${sample},Hicgmlst,Skipped"
        logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
        echo \${logstring} >> /logs/\${logkey}_logs.tsv
    fi

    """
}


process Hicap { 
    container 'ghcr.io/garcia-nacho/top_hicap:v1.1'
    containerOptions '--volume '+params.progress_dir+':/logs'
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
    if [[ -f "Hinf.agent" ]] || [[ -f "Hpar.agent" ]] ; 
    then
        /home/docker/Code/hicapwrapper.sh

        logstring="${sample},Hicap,Completed"
        logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
        echo \${logstring} >> /logs/\${logkey}_logs.tsv

    else
        echo "NoHi" > ${sample}_HiCap.tsv

        logstring="${sample},Hicap,Skipped"
        logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
        echo \${logstring} >> /logs/\${logkey}_logs.tsv

    fi

    """

}

process HinfPBP3 { 
    container 'ghcr.io/garcia-nacho/top_pbp3'
    containerOptions '--volume '+params.progress_dir+':/logs'
    cpus 1
    maxForks = 1

    input:
    path(fastaclean)
    val(sample)
    path(agent)

    output:
    path("*PBP3Mutations.csv"), emit: hipbpresults

    script:

    """
    if test -f "Hinf.agent"; 
    then
       Rscript /home/docker/pbp3/PBP_scanner.R
       mv PBP3Mutations.csv ${sample}_PBP3Mutations.csv


        logstring="${sample},HinfPBP3,Completed"
        logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
        echo \${logstring} >> /logs/\${logkey}_logs.tsv

       
    else
        echo "NoHi" > ${sample}_PBP3Mutations.csv

        logstring="${sample},HinfPBP3,Skipped"
        logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
        echo \${logstring} >> /logs/\${logkey}_logs.tsv

    fi

    """

}

process Seroba { 
    container 'ghcr.io/garcia-nacho/top_seroba:v1.1'
    containerOptions '--volume '+params.progress_dir+':/logs'
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
    if test -f "Spne.agent"; 
    then
        source activate seroba
        seroba runSerotyping /home/docker/seroba/database/ ${r1} ${r2} ${sample}
        mv ${sample}/*.tsv ./${sample}_seroba.tsv
        rm -rf dummy
        conda deactivate

        logstring="${sample},Seroba,Completed"
        logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
        echo \${logstring} >> /logs/\${logkey}_logs.tsv

    else
        echo "NoSpne" > ${sample}_seroba.tsv

        logstring="${sample},Seroba,Skipped"
        logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
        echo \${logstring} >> /logs/\${logkey}_logs.tsv
    fi

    """
}

process STX { 
    container 'ghcr.io/garcia-nacho/top_virfinder:v1.1'
    containerOptions '--volume '+params.progress_dir+':/logs'
    cpus 2
    maxForks = 1
    time '15m'
    errorStrategy 'ignore'
    //errorStrategy { task.exitStatus == 143 ? 'ignore' : 'terminate' }

    input:
    //path(r1)
    //path(r2)
    val(sample)
    path(agent)
    path(rawreads)

    output:
    path("*_fastq_virfinder.json"), emit: stx_results

    script:

    """
    if [[ -f "Ecol.agent" ]] || [[ -f "Shige.agent" ]] ; 
    then
    r1=\$(ls ${sample}*1P.fastq.gz)
    r2=\$(ls ${sample}*2P.fastq.gz)
    r1_count=\$(ls -1 ${sample}*1P.fastq.gz | wc -l) 

    if [ \${r1_count} == 1 ];
    then

    virulencefinder.py -i \${r1} \${r2} -o .

    mv data.json ${sample}_fastq_virfinder.json
    echo "${sample},STX_Fastq,Completed" >> /logs/run_logs.tsv

    else
    
    error "Invalid Sample Name: ${sample}"


    logstring="${sample},STX_Fastq,Failed"
    logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
    echo \${logstring} >> /logs/\${logkey}_logs.tsv

    fi

    else
        echo "NoEcol" > ${sample}_fastq_virfinder.json

        logstring="${sample},STX_Fastq,Skipped"
        logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
        echo \${logstring} >> /logs/\${logkey}_logs.tsv
    fi

    """
}

process STX_Contigs { 
    container 'ghcr.io/garcia-nacho/top_virfinder:v1.1'
    containerOptions '--volume '+params.progress_dir+':/logs'
    cpus 1
    maxForks = 1

    input:
    path(fastaclean)
    val(sample)
    path(agent)

    output:
    path("*_contigs_virfinder.json"), emit: stx_contigs_results

    script:

    """
    if [[ -f "Ecol.agent" ]] || [[ -f "Shige.agent" ]] ; 
    then
    
    virulencefinder.py -i ${sample}_clean_contigs.fasta -o .
    mv data.json ${sample}_contigs_virfinder.json


    logstring="${sample},STX_Contigs,Completed"
    logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
    echo \${logstring} >> /logs/\${logkey}_logs.tsv

    else
        echo "NoEcol" > ${sample}_contigs_virfinder.json

        logstring="${sample},STX_Contigs,Skipped"
        logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
        echo \${logstring} >> /logs/\${logkey}_logs.tsv
    fi

    """
}

process EMMtyper { 
    container 'ghcr.io/garcia-nacho/top_emmtyper:v1.1'
    containerOptions '--volume '+params.progress_dir+':/logs'
    cpus 1
    maxForks = 1

    input:
    path(fastaclean)
    val(sample)
    path(agent)

    output:
    path("*{.tsv,EMM_Allele.fa}"), emit: emm_results

    script:

    """
    if test -f "Spyo.agent"; 
    then

      emmtyper --blast_db /emmdb27022025/emmdb.tfa ${sample}_clean_contigs.fasta > ${sample}_emmtyper.tsv  
      Rscript /home/docker/EMM_Extraction.R
      mv EMM_seqs_extended.fa ${sample}_EMM_Allele.fa


        logstring="${sample},EMMTypper,Completed"
        logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
        echo \${logstring} >> /logs/\${logkey}_logs.ts

    else
      echo "NoSpy" > ${sample}_emmtyper.tsv

        logstring="${sample},EMMTypper,Skipped"
        logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
        echo \${logstring} >> /logs/\${logkey}_logs.ts

    fi

    """

}

process MeningoTyper { 
    container 'ghcr.io/garcia-nacho/top_meningotype:v1.1'
    containerOptions '--volume '+params.progress_dir+':/logs'

    cpus 1
    maxForks = 1

    input:
    path(fastaclean)
    val(sample)
    path(agent)

    output:
    path("*meningotype.txt"), emit: meningotype_results

    script:

    """
    if test -f "Nmen.agent"; 
    then
      meningotype --all ${sample}_clean_contigs.fasta >> ${sample}_meningotype.txt

        logstring="${sample},MeningoTyper,Completed"
        logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
        echo \${logstring} >> /logs/\${logkey}_logs.ts


    else
      echo "NoNmen" > ${sample}_meningotype.txt

        logstring="${sample},MeningoTyper,Skipped"
        logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
        echo \${logstring} >> /logs/\${logkey}_logs.ts

    fi

    """

}

process NGmaster { 
    container 'ghcr.io/garcia-nacho/top_ngmaster:v1.1'
    containerOptions '--volume '+params.progress_dir+':/logs'
    cpus 1
    maxForks = 1

    input:
    path(fastaclean)
    val(sample)
    path(agent)

    output:
    path("*ngmast_results.txt"), emit: ngmaster_results

    script:

    """
    if test -f "Ngon.agent"; 
    then
      ngmaster ${sample}_clean_contigs.fasta >> ${sample}_ngmast_results.txt


        logstring="${sample},NGmaster,Completed"
        logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
        echo \${logstring} >> /logs/\${logkey}_logs.ts

    else
      echo "NoNgon" > ${sample}_ngmast_results.txt

        logstring="${sample},NGmaster,Skipped"
        logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
        echo \${logstring} >> /logs/\${logkey}_logs.ts

    fi

    """

}

process Seqsero { 
    container 'ghcr.io/garcia-nacho/top_seqsero:v1.1'
    containerOptions '--volume '+params.progress_dir+':/logs'
    cpus 2
    maxForks = 1

    input:
    val(sample)
    path(agent)
    path(rawreads)
    //tuple val(sampledummy), path (trimmedR1), path(trimmedR2) 

    output:
    path("*seqsero_results.tsv"), emit: seqsero_results
    path("*seqsero.tar.gz"), emit: seqsero_gzip

    script:

    """
    if test -f "Salmo.agent"; 
    then
    R1=\$(ls ${sample}*1P.fastq.gz)
    R2=\$(ls ${sample}*2P.fastq.gz)
      SeqSero2_package.py -p 2 -t 2 -n ${sample} -i \${R1} \${R2}
      mv \$(ls -d */) SeqSeroResults_allele
      mv SeqSeroResults_allele/SeqSero_result.tsv ./${sample}_seqsero_results.tsv
      tar -zcvf ${sample}_seqsero.tar.gz SeqSeroResults_allele 
      rm -rf SeqSeroResults_allele

      SeqSero2_package.py -m k -p 2 -t 2 -n ${sample} -i \${R1} \${R2}
      mv \$(ls -d */) SeqSeroResults_kmer
      mv SeqSeroResults_kmer/SeqSero_result.tsv ./${sample}_kmer_seqsero_results.tsv
      tar -zcvf ${sample}_kmer_seqsero.tar.gz SeqSeroResults_kmer
      rm -rf SeqSeroResults_kmer
   

        logstring="${sample},Seqsero,Completed"
        logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
        echo \${logstring} >> /logs/\${logkey}_logs.ts
      
    else
      echo "NoSalmo" > ${sample}_seqsero_results.tsv
      echo "NoSalmo" > ${sample}_dummy_seqsero.tar.gz
      echo "NoSalmo" > ${sample}_dummy_kmer_seqsero.tar.gz
      echo "NoSalmo" > ${sample}_kmer_seqsero_results.tsv

        logstring="${sample},Seqsero,Skipped"
        logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
        echo \${logstring} >> /logs/\${logkey}_logs.ts
    fi

    """

}

process Tartrate { 
    container 'ghcr.io/garcia-nacho/top_tartrate:v1.1'
    containerOptions '--volume '+params.progress_dir+':/logs'
    cpus 1
    maxForks = 1

    input:
    path(fastaclean)
    val(sample)
    path(agent)

    output:
    path("*tartrate.txt"), emit: tartrate_results

    script:

    """
    if test -f "Salmo.agent"; 
    then
      tartrate ${sample}_clean_contigs.fasta > ${sample}_tartrate.txt  


        logstring="${sample},Tartrate,Completed"
        logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
        echo \${logstring} >> /logs/\${logkey}_logs.ts

    else
      echo "NoSalmo" > ${sample}_tartrate.txt 

        logstring="${sample},Tartrate,Skipped"
        logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
        echo \${logstring} >> /logs/\${logkey}_logs.ts

    fi

    """

}

process TBpipelineP1{
    container 'ghcr.io/garcia-nacho/top_tbpipeline:v1.1'
    containerOptions '--volume '+params.TBDB+':/mnt/global_collection --volume '+params.progress_dir+':/logs'
    
    cpus 2
    maxForks = 2

    input:
    val(sample)
    path(agent)
    path(rawreads)

    output:
    path("*.tar.gz"), emit: tbpipeline_p1_results

    script:

    """
    if test -f "Myco.agent"; 
    then
      if test -d "/mnt/local_collection/localdb_dummy"; then rm -rf /mnt/local_collection/localdb_dummy; fi
      mkdir ${sample}
      r1=\$(ls ${sample}*1P.fastq.gz)
      r2=\$(ls ${sample}*2P.fastq.gz)
      #trimmomatic PE -basein \${r1} -baseout ${sample}.fastq.gz  ILLUMINACLIP:/home/docker/CommonFiles/adapters/Kapa-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:3:15 MINLEN:36

      mv \${r1} ${sample}/${sample}_R1.fastq.gz
      mv \${r2} ${sample}/${sample}_R2.fastq.gz
      rm ./*.fastq.gz

      niph_tb_pipeline1
      rm -rf COPY_TO_REPORTS
      rm -rf COPY_TO_TB_PIPELINE_DATABASE
      rm -rf ${sample}/*.fastq.gz
      tar -zcvf ${sample}.tar.gz ${sample}
      rm -rf ${sample}

        logstring="${sample},Tartrate,Skipped"
        logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
        echo \${logstring} >> /logs/\${logkey}_logs.ts

    else
      
      mkdir ${sample}_nonTB
      echo "dummy" > ${sample}_nonTB/${sample}_nonTB.txt 
      echo "dummy" > ${sample}_dummy_tbp.tar.gz

        logstring="${sample},TBPipeline,Skipped"
        logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
        echo \${logstring} >> /logs/\${logkey}_logs.ts

    fi

    """
}

process BPEprofiler{
    container 'ghcr.io/garcia-nacho/top_bpprofiler'
    containerOptions '--volume '+params.progress_dir+':/logs'
    cpus 1
    maxForks = 1
    errorStrategy 'retry'
    maxRetries 10

    input:
    path(fastaclean)
    val(sample)
    path(agent)

    output:
    path("*bpe_mlst.csv"), emit: bpeprofiler_results
    path("*.json"), emit: bpeprofiler_json

    script:

    """
    if test -f "Bper.agent"; 
    then
        /home/docker/code/bpe_mlst.sh
        mv BPE_MLST.csv ${sample}_bpe_mlst.csv

        logstring="${sample},BPEprofiler,Completed"
        logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
        echo \${logstring} >> /logs/\${logkey}_logs.ts

    else
 
      echo "NoBper" > ${sample}_bpe_mlst.csv
      echo "NoBper" > ${sample}_dummy.json


        logstring="${sample},BPEprofiler,Skipped"
        logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
        echo \${logstring} >> /logs/\${logkey}_logs.ts

    fi

    """
}

process Diphtoscan{
    
    container 'ghcr.io/garcia-nacho/top_diphtoscan'
    containerOptions '--volume '+params.progress_dir+':/logs'
    
    cpus 1
    maxForks = 1

    input:
    path(fastaclean)
    val(sample)
    path(agent)

    output:
    path("*diphtoscan.csv"), emit: diphto_res

    script:

    """
    if test -f "Cory.agent"; 
    then

        /home/docker/diphtoscan/Diphtorunner.sh

        mv diphtoresults/diphtoscan_results.csv ${sample}_diphtoscan.csv

        logstring="${sample},Diphtoscan,Completed"
        logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
        echo \${logstring} >> /logs/\${logkey}_logs.ts

    else
        echo "NoDiphto" > ${sample}_diphtoscan.csv

        logstring="${sample},Diphtoscan,Skipped"
        logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
        echo \${logstring} >> /logs/\${logkey}_logs.ts
 
    fi

    """
}

process TBprofiler{
    container 'ghcr.io/garcia-nacho/top_tbprofiler'
    containerOptions '--volume '+params.progress_dir+':/logs'
    
    cpus 2
    maxForks = 2

    input:
    val(sample)
    path(agent)
    path(rawreads)

    output:
    path("*tb_profiler.*"), emit: tbprofiler_results

    script:

    """
    if test -f "Myco.agent"; 
    then
      r1=\$(ls ${sample}*1P.fastq.gz)
      r2=\$(ls ${sample}*2P.fastq.gz)
      source activate tb-profiler
      tb-profiler profile -1 \${r1} -2 \${r2} -p ${sample} --csv
      mv results/${sample}.results.csv ${sample}.tb_profiler.csv
      mv results/${sample}.results.json ${sample}.tb_profiler.json
      rm -rf bam
      rm -rf vcf
      rm -rf results  
      conda deactivate
      Rscript /home/docker/Code/TBprofilerparser.R

        logstring="${sample},"TBProfiler",Completed"
        logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
        echo \${logstring} >> /logs/\${logkey}_logs.ts

    else
 
      echo "dummy" > ${sample}_tb_profiler.tsv

        logstring="${sample},TBProfiler,Skipped"
        logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
        echo \${logstring} >> /logs/\${logkey}_logs.ts


    fi

    """
}

process TBpipelineP2{
    container 'ghcr.io/garcia-nacho/top_tbpipeline:v1.1'
    containerOptions '--volume '+params.TBDB+':/mnt/global_collection'
    
    cpus 2
    maxForks = 2

    input:
    path(tbpiperes1)

    output:
    path("*"), emit: tbpipeline_p2_results

    publishDir(
    path: "${params.publishDir}/TB_Pipeline",
    mode: 'copy',
    saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    script:

    """
    Rscript /home/tbuser/Code/TBCleaner.R 
    mkdir topdummy_nonTB
    counttb=\$(ls -dl */ | wc -l)
    countnontb=\$(ls -dl *_nonTB/ | wc -l)

    rm -rf *_nonTB 
 

    if [ \${counttb} -gt \${countnontb} ]
    then

        #parsing
        mkdir COPY_TO_REPORTS
        mkdir COPY_TO_TB_PIPELINE_DATABASE
        niph_tb_pipeline2
        Rscript /home/tbuser/Code/TBParser.R
        
    else
        echo "Non_MTBC_samples_in_the_run" > Non_MTBC_samples_in_the_run
    fi

    """
}

process JonEcoPipe { 
    container 'ghcr.io/garcia-nacho/top_ecoli:v1.1'
    containerOptions '--volume '+params.progress_dir+':/logs'
    cpus 2
    maxForks = 1

    input:
    //path(r1)
    //path(r2)
    val(sample)
    path(agent)
    path(rawreads)

    output:
    path("*_ecopipeline*"), emit: eco_results

    script:

    """
    if [[ -f "Ecol.agent" ]] || [[ -f "Shige.agent" ]] ; 
    then
        r1=\$(ls ${sample}*1P.fastq.gz)
        r2=\$(ls ${sample}*2P.fastq.gz)
        mkdir Fasta Forward Reverse
        mv \${r1} Forward
        mv \${r2} Reverse

        EcoliPipelineTOP.sh Fasta/ Forward/ Reverse/
        Rscript /home/docker/Ecoparsing.R

        mv ecopipeline.csv ${sample}_raw_ecopipeline.csv
        mv EcoliPipelineReceiptFile* ${sample}_ecopipeline_ReceiptFileRaw.txt

        logstring="${sample},BohlinsECPipeline,Completed"
        logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
        echo \${logstring} >> /logs/\${logkey}_logs.ts

    else
        echo "NoEcol" > ${sample}_raw_ecopipeline.csv
        echo "NoEcol" > ${sample}_ecopipeline_ReceiptFileRaw.txt

        logstring="${sample},BohlinsECPipeline,Skipped"
        logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
        echo \${logstring} >> /logs/\${logkey}_logs.ts
    fi

    """
}

process JonEcoPipeFasta { 
    container 'ghcr.io/garcia-nacho/top_ecoli:v1.1'
    containerOptions '--volume '+params.progress_dir+':/logs'
    cpus 1
    maxForks = 1

    input:
    path(fastaclean)
    val(sample)
    path(agent)

    output:
    path("*_ecopipeline*"), emit: eco_results_fasta

    script:

    """
    if [[ -f "Ecol.agent" ]] || [[ -f "Shige.agent" ]] ; 
    then

        mkdir Fasta Forward Reverse
        mv ${sample}_clean_contigs.fasta Fasta/${sample}_clean_contigs.fasta

        EcoliPipelineTOP.sh Fasta/ Forward/ Reverse/
        Rscript /home/docker/Ecoparsing.R

        mv ecopipeline.csv ${sample}_fasta_ecopipeline.csv
        mv EcoliPipelineReceiptFile* ${sample}_ecopipeline_ReceiptFileFasta.txt

        logstring="${sample},BohlinsECPipelineContigs,Completed"
        logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
        echo \${logstring} >> /logs/\${logkey}_logs.ts

    else
        echo "NoEcol" > ${sample}_fasta_ecopipeline.csv
        echo "NoEcol" > ${sample}_ecopipeline_ReceiptFileFasta.txt

        logstring="${sample},BohlinsECPipelineContigs,Skipped"
        logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
        echo \${logstring} >> /logs/\${logkey}_logs.ts
    fi

    """
}

process AmrFinderPlus{
    container 'ghcr.io/garcia-nacho/top_amrfinderplus'
    containerOptions '--volume '+params.progress_dir+':/logs'

    //cpus 1
    //maxForks = 1

    input:
    path(fastaclean)
    val(sample)
    path(agent)

    output:
    path("*_amrfinderplus.tsv"), emit: amrfinder_results

    script:

    """
    if test -f "Ecol.agent"; 
    then
        /home/docker/amrfinder/amrfinder -n *.fasta --plus --organism Escherichia -o ${sample}_amrfinderplus.tsv
    
    elif test -f "Ngon.agent"; 
    then
        /home/docker/amrfinder/amrfinder -n *.fasta --plus --organism Neisseria_gonorrhoeae -o ${sample}_amrfinderplus.tsv

    elif test -f "Nmen.agent"; 
    then
        /home/docker/amrfinder/amrfinder -n *.fasta --plus --organism Neisseria_meningitidis -o ${sample}_amrfinderplus.tsv

    elif test -f "Salmo.agent"; 
    then
        /home/docker/amrfinder/amrfinder -n *.fasta --plus --organism Salmonella -o ${sample}_amrfinderplus.tsv

    elif test -f "Spyo.agent"; 
    then
        /home/docker/amrfinder/amrfinder -n *.fasta --plus --organism Streptococcus_pyogenes -o ${sample}_amrfinderplus.tsv

    elif test -f "Spne.agent"; 
    then
        /home/docker/amrfinder/amrfinder -n *.fasta --plus --organism Streptococcus_pneumoniae -o ${sample}_amrfinderplus.tsv

    elif test -f "Vcol.agent"; 
    then
        /home/docker/amrfinder/amrfinder -n *.fasta --plus --organism Vibrio_cholerae -o ${sample}_amrfinderplus.tsv

    elif test -f "Abau.agent"; 
    then
        /home/docker/amrfinder/amrfinder -n *.fasta --plus --organism Acinetobacter_baumannii -o ${sample}_amrfinderplus.tsv

    elif test -f "Bcep.agent"; 
    then
        /home/docker/amrfinder/amrfinder -n *.fasta --plus --organism Burkholderia_cepacia -o ${sample}_amrfinderplus.tsv
    
    elif test -f "Bmal.agent"; 
    then
        /home/docker/amrfinder/amrfinder -n *.fasta --plus --organism Burkholderia_mallei -o ${sample}_amrfinderplus.tsv

    elif test -f "Bpse.agent"; 
    then
        /home/docker/amrfinder/amrfinder -n *.fasta --plus --organism Burkholderia_pseudomallei -o ${sample}_amrfinderplus.tsv
    
    elif test -f "Cfre.agent"; 
    then
        /home/docker/amrfinder/amrfinder -n *.fasta --plus --organism Citrobacter_freundii -o ${sample}_amrfinderplus.tsv

    elif test -f "Cdif.agent"; 
    then
        /home/docker/amrfinder/amrfinder -n *.fasta --plus --organism Clostridioides_difficile -o ${sample}_amrfinderplus.tsv
    
    elif test -f "Cory.agent"; 
    then
        /home/docker/amrfinder/amrfinder -n *.fasta --plus --organism Corynebacterium_diphtheriae -o ${sample}_amrfinderplus.tsv

    elif test -f "Kpne.agent"; 
    then
        /home/docker/amrfinder/amrfinder -n *.fasta --plus --organism Klebsiella_pneumoniae -o ${sample}_amrfinderplus.tsv

    elif test -f "Paur.agent"; 
    then
        /home/docker/amrfinder/amrfinder -n *.fasta --plus --organism Pseudomonas_aeruginosa -o ${sample}_amrfinderplus.tsv

    elif test -f "Smar.agent"; 
    then
        /home/docker/amrfinder/amrfinder -n *.fasta --plus --organism Serratia_marcescens -o ${sample}_amrfinderplus.tsv

    elif test -f "Saur.agent"; 
    then
        /home/docker/amrfinder/amrfinder -n *.fasta --plus --organism Staphylococcus_aureus -o ${sample}_amrfinderplus.tsv
    
    elif test -f "Saga.agent"; 
    then
        /home/docker/amrfinder/amrfinder -n *.fasta --plus --organism Streptococcus_agalactiae -o ${sample}_amrfinderplus.tsv

    elif test -f "Vpar.agent"; 
    then
        /home/docker/amrfinder/amrfinder -n *.fasta --plus --organism Vibrio_parahaemolyticus -o ${sample}_amrfinderplus.tsv

    elif test -f "Vvul.agent"; 
    then
        /home/docker/amrfinder/amrfinder -n *.fasta --plus --organism Vibrio_vulnificus -o ${sample}_amrfinderplus.tsv

    else
        /home/docker/amrfinder/amrfinder -n *.fasta --plus -o ${sample}_amrfinderplus.tsv
    fi
        logstring="${sample},AMRFinderPlus,Completed"
        logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
        echo \${logstring} >> /logs/\${logkey}_logs.tsv


    """

}

process Integration {
    container 'ghcr.io/garcia-nacho/top_spades:v1.1'
    containerOptions '--volume '+params.progress_dir+':/logs'


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
    path(amrfinderplus_in)
    path(hicap)
    path(hicgmlst)
    path(seroba)
    path(depth)
    path(emmtyp)
    path(meningotype)
    path(ngmast)
    path(ngstar)
    path(stxtp)
    path(stxtp_contig)
    path(seqsero)
    path(seqsero_tar_gz)
    path(tartrate_res)
    path(tbres)
    path(tbprofiler)
    path(ecopipelinefiles)
    path(ecopipelinefilesfasta)
    path(sequencerid)
    path(bpeprofres)
    path(bpeprofjson)
    path(diphtores)
    path(localmlist)
    path(pbpresults)
    
  
    output:
    path ("*")

    script:

    """
    #breakpoint
    multiqc ./
    Rscript /home/docker/CommonFiles/Code/Summarizer.R
    cat Status.tsv >> /logs/run_logs.tsv 
    rm Status.tsv
    mkdir bam
    mv *.bam ./bam
    mv *.bai ./bam
    mkdir fasta
    mv *_clean_contigs.fasta ./fasta

    EMMS=\$(find -iname "*EMM_Allele.fa")
    if [ "\${EMMS}"  !=  ""  ]; then    
        mkdir EMM_Alleles
        cp *EMM_Allele.fa EMM_Alleles
    fi

    mkdir QC
    mv *_raw_contigs.fasta ./QC
    mv *kraken* ./QC
    mv *contigs.stats.csv ./QC
    mv *.txt ./QC
    mv *.zip ./QC
    mv *.json ./QC
    mv *vfdb.tsv ./QC
    mv *ncbi.tsv ./QC
    mv *depth.tsv ./QC
    mv *_seqmlst.csv ./QC
    mv *_rmlst.csv ./QC
    mv *localmlst.tsv ./QC
    mv *_Virulencefactors.csv ./QC
    mv *_STXType.csv ./QC
    mv *_amrfinderplus.tsv ./QC
    mv *_Abricate.csv ./QC
    mv *bpe_mlst.csv ./QC
    mv *tb_profiler* ./QC
    mv *diphtoscan.csv ./QC


    if test -f "*dummy_seqsero_tar.gz"; then rm *dummy_seqsero_tar.gz; fi
    if test -f "*seqsero_tar.gz"; then mv seqsero_tar.gz ./QC; fi

    if test -f "*.tsv"; then mv *.tsv ./QC; fi
    if test -f "*.gbk"; then mv *.gbk ./QC; fi
    if test -f "*.log"; then mv *.log ./QC; fi
    if test -f "*.svg"; then mv *.svg ./QC; fi
    if test -f "*.csv"; then mv *.csv ./QC; fi


    #Rscript /home/docker/CommonFiles/Code/FileParserSummary.R
    
    """
}

workflow {
   sample_reads = Channel.fromFilePairs( params.reads )
   all_raw_reads = Channel.fromPath(params.reads)
   filtered_samples = PreFilter(sample_reads) 
   trimmed=Trimming(filtered_samples)
   ktrim=KrakenTrimmed(trimmed)
   kkraw=KrakenRaw(filtered_samples)
   statuscpu=CPUcounter(trimmed)
   outputspades=Spades(statuscpu.cpu_counts)
   kkcon=KrakenClean(outputspades.spadesraw)
   statuscpu2=CPUcounterBT2(outputspades.spadesraw)
   mapped=Mapping(statuscpu2.cpu_bt2)
   mlst=Rmlst(outputspades.fastasclean, outputspades.sample_name,outputspades.r1spades, outputspades.r2spades)
   abri=Abricate(mlst.clean_contigs_frommlst, mlst.sample_frommlst, mlst.agent)
   hicap=Hicap(mlst.clean_contigs_frommlst, mlst.sample_frommlst, mlst.agent)
   hicgmlst=HiCgmlst(mlst.clean_contigs_frommlst, mlst.sample_frommlst, mlst.agent)

   seroba=Seroba(mlst.r1mlst, mlst.r2mlst, mlst.sample_frommlst, mlst.agent)
   emmtyp=EMMtyper(mlst.clean_contigs_frommlst, mlst.sample_frommlst, mlst.agent)
   meningotype=MeningoTyper(mlst.clean_contigs_frommlst, mlst.sample_frommlst, mlst.agent)
   ngmast=NGmaster(mlst.clean_contigs_frommlst, mlst.sample_frommlst, mlst.agent) 
   ngstar=NGstar(mlst.clean_contigs_frommlst, mlst.sample_frommlst, mlst.agent) 
   //stxtyp=STX(mlst.sample_frommlst, mlst.agent, all_raw_reads.collect())
   stxtyp=STX(mlst.sample_frommlst, mlst.agent, trimmed.map{ it.drop(1) }.flatten().collect())
   stxtypcontig=STX_Contigs(mlst.clean_contigs_frommlst, mlst.sample_frommlst, mlst.agent)
   //seqsero=Seqsero(mlst.sample_frommlst, mlst.agent, all_raw_reads.collect())
   seqsero=Seqsero(mlst.sample_frommlst, mlst.agent, trimmed.map{ it.drop(1) }.flatten().collect())
   tartrate=Tartrate(mlst.clean_contigs_frommlst, mlst.sample_frommlst, mlst.agent)
   //tbpipe1=TBpipelineP1(mlst.sample_frommlst, mlst.agent, all_raw_reads.collect())
   tbpipe1=TBpipelineP1(mlst.sample_frommlst, mlst.agent, trimmed.map{ it.drop(1) }.flatten().collect())
   //tbprof=TBprofiler(mlst.sample_frommlst, mlst.agent, all_raw_reads.collect())
   tbprof=TBprofiler(mlst.sample_frommlst, mlst.agent, trimmed.map{ it.drop(1) }.flatten().collect())
   tbpipe2=TBpipelineP2(tbpipe1.tbpipeline_p1_results.collect())
   //ecopipe=JonEcoPipe(mlst.sample_frommlst, mlst.agent, all_raw_reads.collect())
   ecopipe=JonEcoPipe(mlst.sample_frommlst, mlst.agent, trimmed.map{ it.drop(1) }.flatten().collect())
   ecopipefasta=JonEcoPipeFasta(mlst.clean_contigs_frommlst, mlst.sample_frommlst, mlst.agent)
   amrfindplus=AmrFinderPlus(mlst.clean_contigs_frommlst, mlst.sample_frommlst, mlst.agent)
   bpe=BPEprofiler(mlst.clean_contigs_frommlst, mlst.sample_frommlst, mlst.agent)
   diphto=Diphtoscan(mlst.clean_contigs_frommlst, mlst.sample_frommlst, mlst.agent)
   pbp3hinf=HinfPBP3(mlst.clean_contigs_frommlst, mlst.sample_frommlst, mlst.agent)
   

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
                     amrfindplus.amrfinder_results.collect(),
                     hicap.hicap_results.collect(),
                     hicgmlst.hicgmlst_results.collect(),
                     seroba.seroba_results.collect(),
                     mapped.bt2depth.collect(),
                     emmtyp.emm_results.collect(),
                     meningotype.meningotype_results.collect(),
                     ngmast.ngmaster_results.collect(),
                     ngstar.ngstar_results.collect(),
                     stxtyp.stx_results.collect(),
                     stxtypcontig.stx_contigs_results.collect(),
                     seqsero.seqsero_results.collect(),
                     seqsero.seqsero_gzip.collect(),  
                     tartrate.tartrate_results.collect(), 
                     tbpipe2.tbpipeline_p2_results.collect(),
                     tbprof.tbprofiler_results.collect(),
                     ecopipe.eco_results.collect(),
                     outputspades.sqID.collect(),
                     ecopipefasta.eco_results_fasta.collect(),
                     bpe.bpeprofiler_results.collect(),
                     bpe.bpeprofiler_json.collect(),
                     diphto.diphto_res.collect(),
                     mlst.localmlst.collect(),
                     pbp3hinf.hipbpresults.collect())
}