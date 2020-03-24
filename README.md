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

* AeronScripts: Folder consisting of scripts to generate the graph file namely GraphBuilder.py and ParseGTF.py and additional scripts required by the pipeline
* Binaries: Folder consisting of all the binaries required by the pipeline
* input: Folder containing a sample graph file (in .gfa format)
* config.yaml: Sample config file consisting of all the parameters required by AERON
* Snakefile: Pipeline required to run quantification step of AERON
* Snakefile_fusion: Pipeline required to run the fusion-detection step of AERON

## Version
0.01 - first complete version with quantification and fusion-gene detection ability

## Running
### Graph building
To generate a graph file from a reference sequence, run the following command from the AeronScripts folder:
```
	python GraphBuilder -e Path_to_the_genome_sequence -g Path_to_the_gtf_file -o Output_File  
```
The above command will generate a "gfa" file which can be used for the transcript quantification and gene-fusion detection step. 

A sample graph file generated from annotated transcripts of human (ENSEMBL v92, hg38) is provided in the input file. 

Things to remember:
- The genome sequence file should be in fasta format with each sequence representing a chromosome.
- The number of chromosomes in the sequence fasta file should match the number of chromosomes in the gtf file.
- The chromosome ids in the sequence fasta file should match the chromosome ids in the gtf file.

### Quantification and gene-fusion event detection
1. In the folder Aeron, make a directory titled input  
2. Copy the input files to the "input" folder
3. Input files should include:
	* The input read file(s) in fasta or fastq format. 
	* A graph file in .gfa format.   
	* Annotation file of the species in gtf format. 
	* Reference sequence in fasta or fastq format.
4. Make sure that there is no underscores in the file names, graph .gfa, reads .fq, transcripts .fa  
5. Edit config.yaml, add input file names. An example config file is provided in the repository


parameter | default | explanation
----- | ----- | -----
graph | - | Relative or absolute path of the input graph generated with the graph building script
transcripts | - | Relative or absolute path of the reference transcripts in fasta/q format
reads | - | Name of the input long read fasta/q file(s). Multiple files can be included by placing each of them in its own line. the files should be in the input folder. Do not include "input/". 
gtffile | - | Relative or absolute path of the gtf file used for building the graph
vgpath | - | Relative or absolute path of the vg toolkit binary (https://github.com/vgteam/vg)
fusion_max_error_rate | 0.2 | Maximum allowed error rate for a read to support a fusion. If a read aligns to a predicted fusion transcript with a higher error rate, it is considered to not support the fusion.
fusion_min_score_difference | 200 | Minimum score difference for a read to support a fusion. If a read aligns to both a reference transcript and a predicted fusion transcript, and the score difference is less than the parameter, it is considered to not support the fusion. Higher values lead to a higher precision (higher fraction of predicted fusions are real) at the cost of lower sensitivity (smaller number of real fusions are detected).
seedsize | 17 | Minimum size of an exact match between a read and a transcript to be used for alignment in the quantification pipeline. Higher values lead to faster runtime for quantification but potentially lower accuracy.
maxseeds | 20 | Number of exact matches used for alignment in the quantification pipeline. Higher values lead to a slower runtime for quantification but potentially higher accuracy.

For quantification:

- Run the following command

```
snakemake --cores=no_of_cores all  (experiments were run using 10 cores)  
```
- The quantification results will be in a folder named output  

For gene-fusion detection

- First run the quantification (see above)
- Then run the following command (adjusting the number of cores as needed)
```
snakemake --cores 40 all -s Snakefile_fusion
```
- The fusion detection results will be in a folder named fusionoutput  

## Things to remember:
1.	The graph file given with the repository has been generated using annotated transcripts of human (ENSEMBL v92, hg38). Hence, the file can be used as an input. There is no need to generate a new graph file.

2.	You only need to create the input folder once. Multiple files can be included in the input folder. Also, multiple snakemake runs can be executed using the same input folder.

3.	 An AERON run will create a folder called tmp to store temporary files. The files are named according to the input file given to the program. For instance, a file containing node name to integer mapping for a graph file MyGraph.gfa will be named as MyGraph_nodemapping.txt. If, in another run, the program is using a graph file with the same name (MyGraph), then AERON wont create a new file. Instead it will use the old mapping file. This tmp folder can be removed after every snakemake run 


