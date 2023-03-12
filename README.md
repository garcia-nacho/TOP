# The One Pipeline: *One pipeline to rule them all*
### *One Pipeline to rule them all, One Pipeline to find them, One Pipeline to bring them all, and in the darkness bind them*   
   
[![Linux](https://svgshare.com/i/Zhy.svg)](https://www.linux.org/)   [![Docker](https://badgen.net/badge/icon/docker?icon=docker&label)](https://https://docker.com/) ![NF](https://badgen.net/badge/_/Nextflow/green?icon=terminal)   


The One Pipeline (aka TOP) uses Illumina paired-reads from many different bacterial species to generate a report with:

-Species   
-MLST   
-AMR genes   
-Virulence factors   
-Specific analysis for *E. coli* (Stx genes), *S. pyogenes* (EMM type), *S. pneumonia* (Serotype) and *H. influenze* (Serotype, AMR-related-genes alleles)     
-A lot of quality parameters (N50, L50, Q30, depth, contaminants, etc)
    
 It generates fasta files containing the genomes and bam files containing the reads aligned against these genomes. It saves all the information regarding the tools included in the pipeline.

## Installation   
To install TOP you need a computer with docker and conda installed.    
The Install.sh script will take care of the rest of dependencies. 
Note that the this script is under development at the moment and it will require for you to point to the kraken database. In the future the script will also install the database.     
   
## Running TOP   
The command <code>TOP.sh FolderWithFastq</code> will run the pipeline using the fastq files present in the *./FolderWithFastq* folder.   
The pipeline is expecting a set of subfolders with two fastq files per folder. The fastq files must contain the strings R1/R2 to identify the FW and RV reads. 

## Updating TOP   
The pipeline gets the last nextflow script from this github repository and therefore it is always updated. Same with the docker images, they will we pulled on-the-fly. 
If any major update is eventually implemented you should uninstall TOP and install it again using the information on this repository.
   
## Uninstalling TOP
remove the TOP folder and run the following command <code>conda remove -n top_nf --all</code>

## Under the hood   
Under development
