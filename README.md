# About
Recent sequencing technologies like Pacbio and Oxford Nanopore have made significant progress in sequencing accuracy and sequencing output per run. Thus, there is a need to develop methods that can use long read RNA sequencing data for tasks, where long reads overpower short reads, such as transcript quantification and gene fusion detection. 
AERON is an alignment based pipeline for quantification and detection of gene-fusion events using only long RNA-reads. It uses a state-of-the-art sequence-to-graph aligner to align reads generated from long read sequencing technologies to a reference transcriptome. 
It makes use of a novel way to assign reads to transcripts, based on the position of the mapping of the read on the transcript and the fraction of the read contained in a transcript. AERON also introduces the first long read specific gene-fusion detection algorithm.
AERON was tested on different datasets of varying length and coverage and was found to provide accurate transcript quantification. Also, AERON was able to detect experimentally validated fusion events and novel fusion events which could be validated further.

## Prerequisites
* snakemake (The experiments were run using snakemake version 5.0.0)  
* Python (version >= 2.7)  
* variation graphs tools also known as vg tools. Can be installed from https://github.com/vgteam/vg (experiments were run using vgtools v1.5.0-499-ge8a9bcb)  

## Download

The software can be downloaded by using the following command
```
	git clone https://github.com/SchulzLab/Aeron.git 
```

## Pipeline
The downloaded folder should contain a "snakemake-pipeline" folder which contains the following files and folders:

* Binaries: Folder consisting of all the binaries required by the aligner.
* aeronscripts: Folder consisting of additional scripts required by the aligner
* Graphbuilder: Folder consisting of scripts to generate the graph file namely GraphBuilder.py and ParseGTF.py
* input: Folder containing a tar file consisting of a sample dataset, a sample human transcriptome graph file and a human specie annotation file in gtf format
* config.yaml: Sample config file consisting of all the parameters required by AERON
* Snakemake: Pipeline required to run quantification step of AERON
* Snakefile_fusion: Pipeline required to run the fusion-detection step of AERON

## Version
0.01 - first complete version with quantification and fusion-gene detection ability

## Running
### Graph building
To generate a graph file from a reference sequence, run the following command from the GraphBuilder folder:
'''
	python GraphBuilder -e Path_to_the_genome_sequence -g Path_to_the_genome_sequence -o Output_File  
'''
The above command would generate an "gfa" file which can be used for transcript quantification and gene-fusion detection step

Things to remember:
- The genome sequence file should be in fasta format and should contain the chromosome sequences of the specie.
- The number of chromosomes in the sequence fasta file should match the number of chromosomes in the gtf file.
- The chromosome ids in the sequence fasta file should match the chromosome ids in the gtf file.

### Quantification and gene-fusion detection
1. In the folder of snakemake_pipeline, make a directory titled input  
2. Copy the input files to the "input" folder
3. Input files should include:
	* The input file in fasta or fastq format. 
	* A graph file in .gfa format. 
	* Annotation file of the species in gtf format. 
	* Reference sequence in fasta or fastq format.
4. Make sure that there is no underscores in the file names, graph .gfa, reads .fq, transcripts .fa  
5. Edit config.yaml, add input file names.
6. Example config file is provided in the repository

For quantification:

- Run the following command

```
Snakemake --cores=no_of_cores all  (experiments were run using 10 cores)  
```
- The quantification results will be in a folder named output  


For gene-fusion detection

- Run the following command
```
snakemake --cores 40 all -s Snakefile_fusion
```
- The quantification results will be in a folder named output  

## Things to take care of
- AERON run will create a folder called tmp to store temporary files. The files are named according to the input file given to the program. For instance, a file containing node name to integer mapping for a graph file MyGraph.gfa will be named as MyGraph_nodemapping.txt. If, in another run, the program is using a graph file with the same name (MyGraph), then AERON wont create a new file. Instead it will use the old mapping file 


