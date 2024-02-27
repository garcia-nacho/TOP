# The One Pipeline   
#### *"One Pipeline to rule them all. One Pipeline to find them. One Pipeline to bring them all, and in the darkness bind them"*   

<img src="dalletop.webp"
      type="image/webp"
       width="500" height="500" />   
DallE representation of TOP
      
         
 [![Docker](https://badgen.net/badge/icon/docker?icon=docker&label)](https://https://docker.com/) ![NF](https://badgen.net/badge/_/Nextflow/green?icon=terminal)   


## TOP v1.0
The One Pipeline (aka TOP) is a *multi-agent-multi-analysis* NGS pipeline developed at the Norwegian Institute of Public Health to generate data and reports than can be used in epidemiological analysis.    
TOP uses Illumina paired-reads from different bacterial species and generates a report that includes: 

1. Species identification   
2. MLST (From PUB-MLST) 
3. AMR genes   
4. Virulence factors   
5. Specific analysis for 
* *Escherichia coli / Shigella sp.*: Stx genes, serotypes and specific virulence factors
* *Streptococcus pyogenes*: EMM types
* *Streptococcus pneumoniae*: Serotypes
* *Salmonella sp.*: Serotypes and tartrate fermentation hability
* *Neisseria meningitidis*: Meningotypes
* *Neisseria gonorrhoeae*: NGstar-types and NGmaster-types
* *Haemophilus influenze*: Capsule-type and AMR-related-mutations
* *Mycobacterium tuberculosis*: Lineage and AMR        
6. Quality parameters (N50, L50, Q30, depth, contaminants, etc, Kraken2 reports, etc)
    
TOP generates **fasta** files containing the assembled genomes and bam files containing the reads aligned against these genomes. It also saves all the information regarding the tools included in the pipeline.

TOP is a tool under continuous development. All stable versions are stored as releases in GitHub. We encourage you to use the stable releases since they have been extensiverly tested and the results have been validated.   

## Installing TOP   
To install TOP you need a computer with docker and conda installed.   
<code>git clone https://github.com/garcia-nacho/TOP --branch 1.0</code>   
<code>cd TOP</code>   
<code>./install.sh -p /path/to/TOP</code>

This commands will install TOP in your system and a link to the pipeline will be created in the path defined in /path/to/TOP. Be sure that the path is in your PATH variable, otherwise you will have to call TOP manually (i.e./path/to/TOP/TOP.sh).   
If no path is provided a link to the pipeline will be create in the TOP folder    

During the installation of the pipeline you can specify certain parameters:
      
**Install.sh -c (or --cores)** will define the maximum number of cores available to TOP. If not provided TOP will use up to 10 cores.   
**Install.sh -t (or --tbdb)** will define the location of the *Mycobacterium tuberculosis* database (MTDB).  If not provided, TOP will run without MTDB

The installation script will take care of the rest of dependencies.    
   
TOP have been thoroughly tested in Linux environments, but it should work out the box in MacOS and Windows containing WSL systems.
     
## Running TOP   
To run TOP you just need to execute the TOP.sh command.

TOP.sh can take several arguments:

**TOP.sh -r** is used to define the path the the paired reads. If no path is provided, TOP will try to find paired reads in the current directory.   
**TOP.sh -c** is used to define the number of cores used by the pipeline. It overrides the number of cores defined during the installation.    
**TOP.sh -t** is used to will define the path to the *Mycobacterium tuberculosis* database. It overrides the MTDB defined during the installation.   

The names of the paired reads must be in the following format: *Sample1_R1_001.fastq.gz / Sample1_R2_001.fastq.gz*   

## Updating TOP   
To update TOP (including the main script and all docker images), run <code>TOP.sh --update</code>
      
## Uninstalling TOP
To remove TOP you must run <code>TOP.sh --uninstall</code>

## Future improvements:   
* Speed-up the STX process
* Implement TB-WHOCatalog https://github.com/GTB-tbsequencing/mutation-catalogue-2023/tree/main/Final%20Result%20Files
* Implement an offline version of PathogenWatch
* cgMLST for *H. influenaze*
* Implement TBTyper
* Iclude fastas (instead of fastq) as entry point for the pipeline
* Create a plot to visualize Kraken results
* Dynamic allocation of resources
* Running TB-pipeline from trimmed fastq
   
## Under the hood   
TOP is a Nextflow-based modular pipeline that used docker to isolate and run the different tools included in the analysis. The main TOP.sh is just a wrapper for the nextflow script (TOP.nf), with augmented functionalities such us update, and arguments passed to the nextflow script. The highly modular nature of TOP makes it easy to quickly add new analysis. You just need to add a process in the TOP.nf, include it in the workflow section of the TOP.nf and implement the parsing of the results into the output. The most straight forward to do it is by implementing it in the Summarizer.R script inside the Spades' docker image.     