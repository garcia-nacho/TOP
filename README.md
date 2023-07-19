# The One Pipeline   
#### *One Pipeline to rule them all, One Pipeline to find them, One Pipeline to bring them all, and in the darkness bind them*   
   
 [![Docker](https://badgen.net/badge/icon/docker?icon=docker&label)](https://https://docker.com/) ![NF](https://badgen.net/badge/_/Nextflow/green?icon=terminal)   


The One Pipeline (aka TOP) uses Illumina paired-reads from many different bacterial species to generate a report with:

-Species   
-MLST   
-AMR genes   
-Virulence factors   
-Specific analysis for *E. coli* (Stx genes), *S. pyogenes* (EMM type), *S. pneumonia* (Serotype) and *H. influenze* (Serotype, AMR-related-genes alleles)     
-A lot of quality parameters (N50, L50, Q30, depth, contaminants, etc)
    
 It generates fasta files containing the genomes and bam files containing the reads aligned against these genomes. It saves all the information regarding the tools included in the pipeline.

## Installing TOP   
To install TOP you need a computer with docker and conda installed.   
<code>git clone https://github.com/garcia-nacho/TOP</code>   
<code>cd TOP</code>   
<code>./Install.sh</code>   
The <code>Install.sh</code> script will take care of the rest of dependencies. 
You can now link the <code>TOP.sh</code> script to any folder on your path or run it from there.   
   
## Running TOP   
The command <code>TOP.sh FolderWithFastq</code> will run the pipeline using the fastq files present in the *./FolderWithFastq* folder.   
The pipeline is expecting a set of subfolders with two fastq files per folder. The fastq files must contain the strings R1/R2 to identify the FW and RV reads. 

## Updating TOP main script   
To update TOP script and docker images, run <code>TOP.sh --update</code>
      
## Updating TOP databases and internal tools    
TO BE DONE
   
## Uninstalling TOP
Remove the TOP folder and run the following command <code>conda remove -n top_nf --all</code>

## Under the hood   
Under development
