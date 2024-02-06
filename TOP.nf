#!/usr/bin/env nextflow

// Declare syntax version
nextflow.enable.dsl=2

params.readsfolder = "."
params.publishDir = params.readsfolder+"/TOPresults"
params.threads = 10
params.spadescores = 8
params.forceSp="none"
params.krakenDB="/media/nacho/Data/kraken2_standard_20220926/"
params.TBDB="/mnt/N/NGS/TB_pipeline/TB_pipeline_database/DB/"
params.tempfolder="/media/nacho/Data/temp/toptest/tempdb/"

params.reads=params.readsfolder+"/*/*_{R1,R2}*.fastq.gz"

process Trimming {
 
    container 'ghcr.io/garcia-nacho/top_spades'
    //containerOptions '--volume /media/nacho/Data/kraken2_standard_20220926/:/Kraken2DB'
    errorStrategy 'ignore'
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
    container 'ghcr.io/garcia-nacho/top_spades'
    containerOptions '--volume '+params.krakenDB+':/Kraken2DB'
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

    container 'ghcr.io/garcia-nacho/top_spades'
    //containerOptions '--volume /media/nacho/Data/kraken2_standard_20220926/:/Kraken2DB'

    maxForks = params.spadescores

    input:
    tuple val(sample), path (trimmedR1), path(trimmedR2) 

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
    sqid=\$(gzip -cd ${trimmedR1} | head -n 1)   
    echo \${sqid}  > ${sample}_sequencerID.tsv
    spades.py -o . --careful --cov-cutoff auto -t 1 -1 ${trimmedR1} -2 ${trimmedR2} 
    mv contigs.fasta ${sample}_raw_contigs.fasta
    Rscript /home/docker/CommonFiles/Code/ContigCleaner.R
    mv clean_contigs.fasta ${sample}_clean_contigs.fasta
    mv clean_contigs.stats.csv ${sample}_contigs.stats.csv
    rm -rf ./corrected
    
    
    """
}

process Rmlst {
    container 'ghcr.io/garcia-nacho/top_spades'
    //containerOptions '--volume /media/nacho/Data/kraken2_standard_20220926/:/Kraken2DB'
    errorStrategy 'retry'
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
    curl -s -H "Content-Type: application/json" -X POST "http://rest.pubmlst.org/db/pubmlst_rmlst_seqdef_kiosk/schemes/1/sequence" -d @- > \
    ${sample}_rmlst.json
    Rscript /home/docker/CommonFiles/Code/rmlst_parser.R
    Rscript /home/docker/CommonFiles/Code/seqmlst_parser.R
    #Missing genus 
    source activate mlst
    mlst ${sample}_clean_contigs.fasta > ${sample}_localmlst.tsv
    conda deactivate

    """

}

process Prokka {
    container = 'ghcr.io/garcia-nacho/top_prokka'

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
    container 'ghcr.io/garcia-nacho/top_spades'
    containerOptions '--volume '+params.krakenDB+':/Kraken2DB'
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
    /kraken2-2.1.2/kraken2 --use-names --report  ${sample}_R1_Trimmed_kraken_summaries.tsv --db /Kraken2DB/ ${trimmedR1}  > ${sample}_R1_Trimmed.kraken.tsv
    Rscript /home/docker/CommonFiles/Code/KrakenParser.R

    /kraken2-2.1.2/kraken2 --use-names --report  ${sample}_R2_Trimmed_kraken_summaries.tsv --db /Kraken2DB/ ${trimmedR2}  > ${sample}_R2_Trimmed.kraken.tsv
    Rscript /home/docker/CommonFiles/Code/KrakenParser.R
    """ 

}

process Mapping {
    container 'ghcr.io/garcia-nacho/top_spades'
    //containerOptions '--volume /media/nacho/Data/kraken2_standard_20220926/:/Kraken2DB'

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
    rm ${sample}.sam
    rm ${sample}.bam
    """ 

}

process Integration {
    container 'ghcr.io/garcia-nacho/top_spades'
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
    path(ecopipelinefiles)
    path(ecopipelinefilesfasta)
    path(sequencerid)
    path(localmlist)
  
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
    mv *vfdb.tsv ./QC
    mv *ncbi.tsv ./QC
    mv *depth.tsv ./QC
    mv *_seqmlst.csv ./QC
    mv *_rmlst.csv ./QC
    mv *_Virulencefactors.csv ./QC
    mv *_STXType.csv ./QC
    mv *_Abricate.csv ./QC

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

process Abricate { 
    container 'ghcr.io/garcia-nacho/top_abricate'
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
        #abricate-get_db --db ncbi --force
        #abricate-get_db --db vfdb --force
        #abricate-get_db --db plasmidfinder --force
        abricate --db vfdb --quiet *.fasta > ${sample}_vfdb.tsv 
        abricate --db HinfFtsI --quiet *.fasta > ${sample}_HinfFtsI.tsv
        abricate --db HinfGyrSubA --quiet *.fasta > ${sample}_HinfGyrSubA.tsv
        abricate --db HinfTopoIVsubA --quiet *.fasta > ${sample}_HinfTopoIVSubA.tsv
        abricate --db ncbi --quiet *.fasta > ${sample}_ncbi.tsv
        abricate --db plasmidfinder --quiet *.fasta > ${sample}_plasmidfinder.tsv
        #Integration Abricate
    else
        #abricate-get_db --db ncbi --force
        #abricate-get_db --db vfdb --force
        #abricate-get_db --db plasmidfinder --force
        abricate --db vfdb --quiet *.fasta > ${sample}_vfdb.tsv
        abricate --db ncbi --quiet *.fasta > ${sample}_ncbi.tsv
        abricate --db plasmidfinder --quiet *.fasta > ${sample}_plasmidfinder.tsv
        #Dummy file
    fi

    #Test if the files can be empty or they just dont exist
    
    Rscript /home/docker/Code/AbricateParser.R
    mv Abricate.csv ${sample}_Abricate.csv

    """

}

process NGstar { 
    container 'ghcr.io/garcia-nacho/top_ngstar'
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
      python3 /home/docker/pyngSTar/pyngSTar.py -f -i ${sample}_clean_contigs_nocap.fasta -p /home/docker/pyngSTar/pyngSTarDB_02012024/ -o ${sample}_ngstar_results.tsv 
      
    else
      echo "NoNgon" > ${sample}_ngstar_results.tsv

    fi

    """

}

process Hicap { 
    container 'ghcr.io/garcia-nacho/top_hicap'
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

    else
        echo "NoHi" > ${sample}_HiCap.tsv
        #Dummy Hicap file

    fi

    """

}

process Seroba { 
    container 'ghcr.io/garcia-nacho/top_seroba'
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

    else
        echo "NoSpne" > ${sample}_seroba.tsv
    fi

    """
}

process STX { 
    container 'ghcr.io/garcia-nacho/top_virfinder'
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
    if test -f "Ecol.agent"; 
    then
    r1=\$(ls ${sample}_R1*.fastq.gz)
    r2=\$(ls ${sample}_R2*.fastq.gz)
    r1_count=\$(ls -1 ${sample}_R1*.fastq.gz | wc -l) 

    if [ \${r1_count} == 1 ];
    then

    virulencefinder.py -i \${r1} \${r2} -o .

    mv data.json ${sample}_fastq_virfinder.json

    else
    
    error "Invalid Sample Name: ${sample}"

    fi

    else
        echo "NoEcol" > ${sample}_fastq_virfinder.json
    fi

    """
}

process STX_Contigs { 
    container 'ghcr.io/garcia-nacho/top_virfinder'
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
    if test -f "Ecol.agent"; 
    then
    
    virulencefinder.py -i ${sample}_clean_contigs.fasta -o .
    mv data.json ${sample}_contigs_virfinder.json

    else
        echo "NoEcol" > ${sample}_contigs_virfinder.json
    fi

    """
}

process EMMtyper { 
    container 'ghcr.io/garcia-nacho/top_emmtyper'
    cpus 1
    maxForks = 1

    input:
    path(fastaclean)
    val(sample)
    path(agent)

    output:
    path("*.tsv"), emit: emm_results

    script:

    """
    if test -f "Spyo.agent"; 
    then
      emmtyper ${sample}_clean_contigs.fasta > ${sample}_emmtyper.tsv  

    else
      echo "NoSpy" > ${sample}_emmtyper.tsv

    fi

    """

}

process MeningoTyper { 
    container 'ghcr.io/garcia-nacho/top_meningotype'
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

    else
      echo "NoNmen" > ${sample}_meningotype.txt

    fi

    """

}

process NGmaster { 
    container 'ghcr.io/garcia-nacho/top_ngmaster'
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

    else
      echo "NoNgon" > ${sample}_ngmast_results.txt

    fi

    """

}

process Seqsero { 
    container 'ghcr.io/garcia-nacho/top_seqsero'
    cpus 2
    maxForks = 1

    input:
    val(sample)
    path(agent)
    path(rawreads)

    output:
    path("*seqsero_results.tsv"), emit: seqsero_results
    path("*seqsero.tar.gz"), emit: seqsero_gzip

    script:

    """
    if test -f "Salmo.agent"; 
    then
    R1=\$(ls ${sample}*R1*.fastq.gz)
    R2=\$(ls ${sample}*R2*.fastq.gz)
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
      
    else
      echo "NoSalmo" > ${sample}_seqsero_results.tsv
      echo "NoSalmo" > ${sample}_dummy_seqsero.tar.gz
      echo "NoSalmo" > ${sample}_dummy_kmer_seqsero.tar.gz
      echo "NoSalmo" > ${sample}_kmer_seqsero_results.tsv
    fi

    """

}

process Tartrate { 
    container 'ghcr.io/garcia-nacho/top_tartrate'
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

    else
      echo "NoSalmo" > ${sample}_tartrate.txt 

    fi

    """

}

process TBpipeline{
    container 'ghcr.io/garcia-nacho/top_tbpipeline'
    containerOptions '--volume '+params.TBDB+':/mnt/global_collection --volume ' +params.tempfolder+':/mnt/local_collection'
    
    cpus 2
    maxForks = 2

    publishDir "${params.tempfolder}", pattern: "localdb_*", mode: 'copy'

    input:
    val(sample)
    path(agent)
    path(rawreads)

    output:
    path("*tbp.csv"), emit: tbpipeline_results
    path("*tbp_bulk.tar.gz"), emit: tbpipeline_bulk
    path("localdb_*"), emit: tbpipeline_localdb

    script:

    """
    if test -f "Myco.agent"; 
    then
      if test -d "/mnt/local_collection/localdb_dummy"; then rm -rf /mnt/local_collection/localdb_dummy; fi
      mkdir ${sample}
    
      r1=\$(ls ${sample}_R1*.fastq.gz)
      r2=\$(ls ${sample}_R2*.fastq.gz)
       
      mv \${r1} ${sample}
      mv \${r2} ${sample}
      rm ./*.fastq.gz

      niph_tb_pipeline
      Rscript /home/tbuser/Code/TBParser.R
      mkdir localdb_${sample}
      mkdir tb_bulk
      mkdir tb_bulk/${sample}

      rm ${sample}/*.fastq.gz
      mv ${sample}/* tb_bulk/${sample}/
      
      cp -R COPY_TO_TB_PIPELINE_DATABASE/*/* localdb_${sample}  
      mv COPY_TO_TB_PIPELINE_DATABASE/ tb_bulk
      mv COPY_TO_REPORTS/ tb_bulk
      mv *.log tb_bulk
      mv *.nwk tb_bulk
      mv TB* tb_bulk
      tar -czvf ${sample}_tbp_bulk.tar.gz tb_bulk
      rm -rf ${sample}

    else
      if test -d "/mnt/local_collection/localdb_dummy"; then rm -rf /mnt/local_collection/localdb_dummy; fi
      echo "dummy" > ${sample}_dummy_tbp_bulk.tar.gz
      echo "NoMyco" > ${sample}_tbp.csv
      mkdir localdb_dummy
      echo "dummy" > localdb_dummy/dummy.txt   

    fi

    """
}

process TBpipelineP1{
    container 'ghcr.io/garcia-nacho/top_tbpipeline'
    containerOptions '--volume '+params.TBDB+':/mnt/global_collection'
    
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
    
      r1=\$(ls ${sample}_R1*.fastq.gz)
      r2=\$(ls ${sample}_R2*.fastq.gz)
       
      mv \${r1} ${sample}
      mv \${r2} ${sample}
      rm ./*.fastq.gz

      niph_tb_pipeline1
      rm -rf COPY_TO_REPORTS
      rm -rf COPY_TO_TB_PIPELINE_DATABASE
      rm -rf ${sample}/*.fastq.gz
      tar -zcvf ${sample}.tar.gz ${sample}
      rm -rf ${sample}

    else
      
      mkdir ${sample}_nonTB
      echo "dummy" > ${sample}_nonTB/${sample}_nonTB.txt 
      echo "dummy" > dummy_tbp.tar.gz

    fi

    """
}

process TBpipelineP2{
    container 'ghcr.io/garcia-nacho/top_tbpipeline'
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
    rm -rf topdummy_nonTB
    countnontb=\$((\$countnontb-1))
    
    if [ \${counttb} -gt \${countnontb} ]
    then
        if [ \${countnontb} -gt 0 ]
        then
            rm -rf *_nonTB 
        fi
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
    container 'ghcr.io/garcia-nacho/top_ecoli'
    cpus 2
    maxForks = 1

    input:
    //path(r1)
    //path(r2)
    val(sample)
    path(agent)
    path(rawreads)

    output:
    path("*_ecopipeline.csv"), emit: eco_results

    script:

    """
    if [[ -f "Ecol.agent" ]] || [[ -f "Shige.agent" ]] ; 
    then
        r1=\$(ls ${sample}_R1*.fastq.gz)
        r2=\$(ls ${sample}_R2*.fastq.gz)
        mkdir Fasta Forward Reverse
        mv \${r1} Forward
        mv \${r2} Reverse

        EcoliPipelineTOP.sh Fasta/ Forward/ Reverse/
        Rscript /home/docker/Ecoparsing.R

        mv ecopipeline.csv ${sample}_raw_ecopipeline.csv
        mv EcoliPipelineReceiptFile* ${sample}_ecopipeline_ReceiptFileRaw.txt

    else
        echo "NoEcol" > ${sample}_raw_ecopipeline.csv
        echo "NoEcol" > ${sample}_ecopipeline_ReceiptFileRaw.txt
    fi

    """
}

process JonEcoPipeFasta { 
    container 'ghcr.io/garcia-nacho/top_ecoli'
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

    else
        echo "NoEcol" > ${sample}_fasta_ecopipeline.csv
        echo "NoEcol" > ${sample}_ecopipeline_ReceiptFileFasta.txt
    fi

    """
}

workflow {
   sample_reads = Channel.fromFilePairs( params.reads )
   all_raw_reads = Channel.fromPath(params.reads)
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
   emmtyp=EMMtyper(mlst.clean_contigs_frommlst, mlst.sample_frommlst, mlst.agent)
   meningotype=MeningoTyper(mlst.clean_contigs_frommlst, mlst.sample_frommlst, mlst.agent)
   ngmast=NGmaster(mlst.clean_contigs_frommlst, mlst.sample_frommlst, mlst.agent) 
   ngstar=NGstar(mlst.clean_contigs_frommlst, mlst.sample_frommlst, mlst.agent) 
   //stxtyp=STX(mlst.r1mlst, mlst.r2mlst, mlst.sample_frommlst, mlst.agent, all_raw_reads)
   stxtyp=STX(mlst.sample_frommlst, mlst.agent, all_raw_reads.collect())
   stxtypcontig=STX_Contigs(mlst.clean_contigs_frommlst, mlst.sample_frommlst, mlst.agent)
   seqsero=Seqsero(mlst.sample_frommlst, mlst.agent, all_raw_reads.collect())
   tartrate=Tartrate(mlst.clean_contigs_frommlst, mlst.sample_frommlst, mlst.agent)
   //tbpipe=TBpipeline(mlst.sample_frommlst, mlst.agent, all_raw_reads.collect())
   tbpipe1=TBpipelineP1(mlst.sample_frommlst, mlst.agent, all_raw_reads.collect())
   tbpipe2=TBpipelineP2(tbpipe1.tbpipeline_p1_results.collect())
   ecopipe=JonEcoPipe(mlst.sample_frommlst, mlst.agent, all_raw_reads.collect())
   ecopipefasta=JonEcoPipeFasta(mlst.clean_contigs_frommlst, mlst.sample_frommlst, mlst.agent)

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
                     ecopipe.eco_results.collect(),
                     outputspades.sqID.collect(),
                     ecopipefasta.eco_results_fasta.collect(),
                     mlst.localmlst.collect() )
}