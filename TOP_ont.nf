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
params.devrun="No"
params.progress_dir = "/media/nacho/Data/temp/toptest/testont/logs"
params.reads=params.readsfolder+"/*"
params.softtrimming = "No"

include {  KrakenClean; Rmlst; Abricate; Hicap;
           HiCgmlst; EMMtyper; BPEprofiler;  MeningoTyper;STX_Contigs;Tartrate;
           JonEcoPipeFasta;
           AmrFinderPlus;
           Diphtoscan; HinfPBP3; NGmaster;NGstar; 
           Seqsero; STX; Seroba; TBpipelineP1; TBprofiler; TBpipelineP2; JonEcoPipe } from './TOP.nf'


workflow {


   root_ch = Channel.value( file(params.readsfolder) )

   prep = Renaming(root_ch)

   merged_with_id_ch = (
   prep.merged_fastqs
    .flatten()
    .map { f ->
      def sample_id = f.name.replaceFirst(/_merged\.fastq\.gz$/, '')
      tuple(sample_id, f)
    }
   )

   assembly = ONTassembly(merged_with_id_ch)
   nanoplot = Nanoplot(merged_with_id_ch)
   cleaning = Cleaning(assembly.withid)
   kkraw = KrakenONT(merged_with_id_ch)
   kkcon = KrakenClean(cleaning.dummyparsed)
   mlst=Rmlst(cleaning.fastasclean, cleaning.sample_name,cleaning.dummyr1, cleaning.dummyr2)
   abri=Abricate(mlst.clean_contigs_frommlst, mlst.sample_frommlst, mlst.agent)
   hicap=Hicap(mlst.clean_contigs_frommlst, mlst.sample_frommlst, mlst.agent)
   hicgmlst=HiCgmlst(mlst.clean_contigs_frommlst, mlst.sample_frommlst, mlst.agent)
   emmtyp=EMMtyper(mlst.clean_contigs_frommlst, mlst.sample_frommlst, mlst.agent)
   bpe=BPEprofiler(mlst.clean_contigs_frommlst, mlst.sample_frommlst, mlst.agent)
   ecopipefasta=JonEcoPipeFasta(mlst.clean_contigs_frommlst, mlst.sample_frommlst, mlst.agent)
   amrfindplus=AmrFinderPlus(mlst.clean_contigs_frommlst, mlst.sample_frommlst, mlst.agent)
   diphto=Diphtoscan(mlst.clean_contigs_frommlst, mlst.sample_frommlst, mlst.agent)
   pbp3hinf=HinfPBP3(mlst.clean_contigs_frommlst, mlst.sample_frommlst, mlst.agent)   
   meningotype=MeningoTyper(mlst.clean_contigs_frommlst, mlst.sample_frommlst, mlst.agent)
   ngmast=NGmaster(mlst.clean_contigs_frommlst, mlst.sample_frommlst, mlst.agent) 
   ngstar=NGstar(mlst.clean_contigs_frommlst, mlst.sample_frommlst, mlst.agent) 

   stxtypcontig=STX_Contigs(mlst.clean_contigs_frommlst, mlst.sample_frommlst, mlst.agent)
   tartrate=Tartrate(mlst.clean_contigs_frommlst, mlst.sample_frommlst, mlst.agent)

   dumfq=DummyFastq(mlst.sample_frommlst, mlst.clean_contigs_frommlst, mlst.agent)

   seqsero=Seqsero(dumfq.sample, dumfq.agent, dumfq.dummyreads)
   stxtyp=STX(dumfq.sample, dumfq.agent, dumfq.dummyreads)
   seroba=Seroba(dumfq.r1, dumfq.r2, dumfq.sample, dumfq.agent)
   tbpipe1=TBpipelineP1(dumfq.sample, dumfq.agent, dumfq.dummyreads)
   tbprof=TBprofiler(dumfq.sample, dumfq.agent, dumfq.dummyreads)
   tbpipe2=TBpipelineP2(tbpipe1.tbpipeline_p1_results.collect())
   ecopipe=JonEcoPipe(dumfq.sample, dumfq.agent, dumfq.dummyreads)

   integ=Integration(
                     nanoplot.tosum.collect(),
                     cleaning.fastasclean.collect(),
                     cleaning.contigstats.collect(),
                     cleaning.b2dum.collect(),
                     kkraw.collect(),
                     kkcon.collect(),
                     nanoplot.dummyfastqc.collect(),
                     assembly.raw.collect(),
                     mlst.mlstresults.collect(),
                     abri.abricate_results.collect(),
                     amrfindplus.amrfinder_results.collect(),
                     hicap.hicap_results.collect(),
                     hicgmlst.hicgmlst_results.collect(),
                     seroba.seroba_results.collect(),
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
                     ecopipefasta.eco_results_fasta.collect(),
                     bpe.bpeprofiler_results.collect(),
                     bpe.bpeprofiler_json.collect(),
                     diphto.diphto_res.collect(),
                     mlst.localmlst.collect(),
                     dumfq.r1.collect(),
                     dumfq.r2.collect(),
                     assembly.seqncerid.collect(),
                     assembly.medakabam.collect(),
                     pbp3hinf.hipbpresults.collect())


}



process Renaming {

   container 'ghcr.io/garcia-nacho/top_ont'
   containerOptions '--volume '+params.progress_dir+':/logs'
   maxForks = 1
   
   publishDir "${params.publishDir}", mode: 'copy'

  input:
    path reads_root

  output:
    path ("merged/*_merged.fastq.gz") , emit: merged_fastqs 
   
   script:
   """
   home=\$(pwd)
   cd ${params.readsfolder}
   Rscript /home/docker/ont_renamer.R
   mv merged \${home}/merged 

   """
    
}

process ONTassembly {

   container 'ghcr.io/garcia-nacho/top_ont'
   containerOptions '--volume '+params.progress_dir+':/logs'
   maxForks = params.spadescores
   tag "$sid"

   input:
      tuple val(sid), path(reads)
   
   output:
      tuple val(sid), path ("*raw_contigs.fasta"), emit: withid
      path ("*raw_contigs.fasta"), emit: raw
      path ("*sequencerID.tsv"),emit: seqncerid
      path ("*_medaka.bam*"), emit: medakabam
   
   script:
   """

      gunzip -c "${reads}" > "${sid}_merged.fastq"
      sqid=\$(gzip -cd ${reads} | head -n 1)   
      echo \${sqid}  > ${sid}_sequencerID.tsv

      cat ${sid}_merged.fastq | NanoFilt -q 10 -l 3000 > ${sid}_filtered.fastq
      cat ${sid}_merged.fastq | NanoFilt -q 8 -l 300 > ${sid}_raw.fastq
      rm ${sid}_merged.fastq
      mkdir fasta

      flye --nano-raw ${sid}_filtered.fastq --out-dir fasta/${sid} --threads 1
      mv fasta/${sid}/assembly.fasta ${sid}_flye.fasta
      minimap2 -x map-ont ${sid}_flye.fasta ${sid}_raw.fastq > ${sid}_overlaps.paf
      racon -t 1 ${sid}_raw.fastq ${sid}_overlaps.paf ${sid}_flye.fasta > ${sid}_racon.fasta
      rm ${sid}_overlaps.paf

      source activate medaka
      medaka_consensus -i ${sid}_raw.fastq -d ${sid}_racon.fasta -o ${sid}_medaka -t 1
      mv ${sid}_medaka/calls_to_draft.bam ${sid}_medaka.bam
      mv ${sid}_medaka/calls_to_draft.bam.bai ${sid}_medaka.bam.bai
      mv ${sid}_medaka/consensus.fasta ${sid}_raw_contigs.fasta

      conda deactivate

   """
    
}


process Nanoplot {

   container 'ghcr.io/garcia-nacho/top_ont'
   containerOptions '--volume '+params.progress_dir+':/logs'
   maxForks = params.threads
   errorStrategy 'ignore'
    
   tag "$sid"
   
   input:
      tuple val(sid), path(reads)
   
   output:
      tuple val(sid), path ("*NanoStats.txt"), emit: nanostats
      tuple val(sid), path ("*.html"), emit: nanoplots
      path ("*{NanoStats.txt, .html}"), emit: tosum
      path("*fastqc.zip"), emit: dummyfastqc
   
   script:
   """
   set +e
   NanoPlot --fastq ${reads} --outdir nanoplot_${sid}
   mv nanoplot_${sid}/NanoStats.txt ./${sid}_NanoStats.txt
   mv nanoplot_${sid}/WeightedLogTransformed_HistogramReadlength.html ./${sid}_WeightedLogTransformed_HistogramReadlength.html
   mv nanoplot_${sid}/Non_weightedLogTransformed_HistogramReadlength.html ./${sid}_Non_weightedLogTransformed_HistogramReadlength.html

   Rscript /home/docker/fastqc_dummyfier.R

   """
}

process Cleaning {

   container 'ghcr.io/garcia-nacho/top_ont'
   containerOptions '--volume '+params.progress_dir+':/logs'
   maxForks = 1

   tag "$sid"

  input:
    tuple val(sid), path (fastas)

  output:
    path ("*clean_contigs.fasta"), emit: fastasclean
    val(sid), emit: sample_name
    path ("dummy_R1.fastq.gz"), emit: dummyr1
    path ("dummy_R2.fastq.gz"), emit: dummyr2
    tuple val(sid), path (fastas),  path ("*clean_contigs.fasta"), path ("dummy_R1.fastq.gz"), path("dummy_R2.fastq.gz"), emit: dummyparsed 
    path ("*clean_contigs.stats.csv"), emit: contigstats
    path ("*Bowtie2summary.txt"), emit: b2dum
   
   script:
   """
   Rscript /home/docker/ont_cleaner.R
   mv clean_contigs.fasta ${sid}_clean_contigs.fasta
   mv clean_contigs.stats.csv ${sid}_clean_contigs.stats.csv

   #Create dummy R1 and R2
   echo 'NA' > dummy_R1.fastq.gz
   echo 'NA' > dummy_R2.fastq.gz 
   echo 'NA' > ${sid}_Bowtie2summary.txt


   """
    
}

process KrakenONT {
    container 'ghcr.io/garcia-nacho/top_spades:v1.1'
    containerOptions '--volume '+params.krakenDB+':/Kraken2DB --volume '+params.progress_dir+':/logs'
    maxForks = 1
    
    tag "$sample"

    input:
      tuple val(sample), path(reads)

    output:
      path ("*{.tsv,.csv,.zip}")

    script:

    """

    /kraken2-2.1.2/kraken2 --use-names --memory-mapping --report  ${sample}_R2_Raw_kraken_summaries.tsv --db /Kraken2DB/ ${reads}  > ${sample}_Raw.kraken.tsv
    Rscript /home/docker/CommonFiles/Code/KrakenParser.R

    cp ${sample}_Raw.resultskraken.csv ${sample}_Trimmed.resultskraken.csv

    logstring="${sample},KrakenONT,Completed"
    logkey=\$(echo -n \${logstring} | md5sum | cut -d ' ' -f1 | cut -c1-8)
    echo \${logstring} >> /logs/\${logkey}_logs.tsv

    """ 

}

process DummyFastq {

   container 'ghcr.io/garcia-nacho/top_ont'
   containerOptions '--volume '+params.progress_dir+':/logs'
   maxForks = 1

   tag "$sid"

   input:
      val(sid)
      path (fastas)
      path (agent)

   output:
      path ("*.fastq.gz"), emit: dummyreads
      path ("*R1_001.fastq.gz"), emit: r1
      path ("*R2_001.fastq.gz"), emit: r2
      val (sid), emit: sample
      path (agent), emit:agent


   script:

   """
   wgsim -N 1000000 -1 150 -2 150 -d 350 -S 42 -e 0.0 ${fastas} ${sid}_R1_001.fastq ${sid}_R2_001.fastq
   gzip ${sid}_R1_001.fastq
   gzip ${sid}_R2_001.fastq


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
      path(nanoplot)
      path(fastaclean)
      path(stats)
      path(b2dummy)
      path(kkn)
      path(kkc)
      path(fasqc)
      path(fastaraw)
      path(mlst)
      path(abricate)
      path(amrfinderplus)
      path(hicap)
      path(hicgmlis)
      path(seroba)
      path(emmtype)
      path(menigotype)
      path(ngmast)
      path(ngstar)
      path(stxtyp)
      path(stxtypecontig)
      path(seqsero)
      path(seqserozip)
      path(tartrat)
      path(tbpide2)
      path(tbprof)
      path(ecopipe1)
      path(ecopipe2)
      path(bperes)
      path(bepjson)
      path(diptho)
      path(mlst)
      path(dummyr1)
      path(dummyr2)
      path(sqncerid)
      path(medakabam)
      path(pbp3hif)    
  
    output:
       path ("*")

    script:

    """
    #breakpoint

    
    Rscript /home/docker/CommonFiles/Code/Summarizer.R
    cat Status.tsv >> /logs/run_logs.tsv 
    rm Status.tsv

    mkdir illumina_fastq
    mv *.fastq.gz ./illumina_fastq

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

    mkdir bam
    mv *.bam ./bam
    mv *.bai ./bam



    #Rscript /home/docker/CommonFiles/Code/FileParserSummary.R
    
    """
}