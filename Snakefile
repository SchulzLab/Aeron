configfile: "config.yaml"
VGPATH = config["vgpath"]
SCRIPTPATH=config["scripts"]
ALIGNERBINPATH = config["binaries"]
READFILE_FULLNAME = config["reads"]
TRANSCRIPTFILE_FULLNAME = config["transcripts"]
GRAPHFILE_FULLNAME = config["graph"]
SEEDSIZE = config["seedsize"]
MAXSEEDHITS = config["maxseeds"]
ALIGNERBANDWIDTH = config["aligner_bandwidth"]
GTFPATH = config["gtffile"]
ALIGNMENTSELECTION = config["alignment_selection"]
ECUTOFF = config["alignment_E_cutoff"]

if isinstance(READFILE_FULLNAME, str): READFILE_FULLNAME = [READFILE_FULLNAME]

READFILE_NAME = [path.split('.')[0] for path in READFILE_FULLNAME]
TRANSCRIPTFILE_NAME = TRANSCRIPTFILE_FULLNAME.split('.')[0]
GRAPHFILE_NAME = GRAPHFILE_FULLNAME.split('.')[0]

def readfile_withending(wildcards):
	for r in READFILE_FULLNAME:
		if r.split('.')[0] == wildcards.reads:
			return "input/" + r
	if wildcards.reads == TRANSCRIPTFILE_NAME: return "input/" + TRANSCRIPTFILE_FULLNAME
	assert False

rule all:
	input:
		expand("output/aln_{reads}_{graph}_all.gam", reads=READFILE_NAME, graph=GRAPHFILE_NAME),
		expand("output/aln_{reads}_{graph}_selected.gam", reads=READFILE_NAME, graph=GRAPHFILE_NAME),
		expand("output/aln_{reads}_{graph}_full_length.gam", reads=READFILE_NAME, graph=GRAPHFILE_NAME),
		#expand("output/alignmentstats_{reads}_{graph}.txt", reads=READFILE_NAME, graph=GRAPHFILE_NAME),
		expand("output/aln_{transcript}_{graph}_full_length.gam", transcript=TRANSCRIPTFILE_NAME, graph=GRAPHFILE_NAME),
		#expand("output/alignmentstats_{transcript}_{graph}.txt", transcript=TRANSCRIPTFILE_NAME, graph=GRAPHFILE_NAME),
		expand("output/matrix_{reads}_{transcript}_{graph}_all.txt", reads=READFILE_NAME, transcript=TRANSCRIPTFILE_NAME, graph=GRAPHFILE_NAME),
		#expand("output/matrix_{reads}_{transcript}_{graph}_bestmatch.txt", reads=READFILE_NAME, transcript=TRANSCRIPTFILE_NAME, graph=GRAPHFILE_NAME),
		expand("output/matrixstats_{reads}_{transcript}_{graph}.txt", reads=READFILE_NAME, transcript=TRANSCRIPTFILE_NAME, graph=GRAPHFILE_NAME),
		expand("output/CountMatrix_{reads}_{transcript}_{graph}.txt", reads=READFILE_NAME, transcript=TRANSCRIPTFILE_NAME, graph=GRAPHFILE_NAME)

rule align:
	input:
		graph = "input/{graph}.gfa",
		reads = readfile_withending
	output:
		"output/aln_{reads}_{graph}_all.gam"
	benchmark:
		"benchmark/aln_{reads}_{graph}_all.txt"
	threads: 15
	shell:
		"/usr/bin/time -v {ALIGNERBINPATH}/GraphAligner -g {input.graph} -f {input.reads} --try-all-seeds --seeds-mxm-length {SEEDSIZE} --seeds-mem-count {MAXSEEDHITS} --seeds-mxm-cache-prefix tmp/seedcache -a {output} -t {threads} -b {ALIGNERBANDWIDTH} {ALIGNMENTSELECTION} 1> tmp/aligner_stdout.txt 2> tmp/aligner_stderr.txt"

rule postprocess:
	input:
		all_alns = "output/aln_{reads}_{graph}_all.gam",
		reads = readfile_withending
	output:
		selected_alns = "output/aln_{reads}_{graph}_selected.gam",
		full_len_alns = "output/aln_{reads}_{graph}_full_length.gam",
		summary = "tmp/run_{reads}_{graph}_summary.txt",
	benchmark:
		"benchmark/postprocess_{reads}_{graph}.txt"
	shell:
		"{ALIGNERBINPATH}/Postprocess {input.all_alns} {input.reads} {output.selected_alns} {output.full_len_alns} {output.summary}"

rule pick_longest:
	input:
		"output/aln_{reads}_{graph}_selected.gam"
	output:
		"output/aln_{reads}_{graph}_longest.gam"
	benchmark:
		"benchmark/pick_longest_{reads}_{graph}.txt"
	shell:
		"{ALIGNERBINPATH}/SelectLongestAlignment {input} {output}"

rule assign_reads_to_transcripts:
	input:
		readalns = "output/aln_{reads}_{graph}_selected.gam",
		transcripts = "output/aln_{transcripts}_{graph}_full_length.gam",
		readfa = readfile_withending
	output:
		"output/matrix_{reads}_{transcripts}_{graph}_all.txt"
	benchmark:
		"benchmark/assign_reads_to_transcript_{reads}_{transcripts}_{graph}.txt"
	shell:
		"{ALIGNERBINPATH}/AlignmentSubsequenceIdentity {input.transcripts} {input.readalns} {input.readfa} > {output}"

#rule find_best_assignments:
#	input:
#		"output/matrix_{runid}_all.txt"
#	output:
#		"output/matrix_{runid}_bestmatch.txt"
#	benchmark:
#		"benchmark/bestmatch_{runid}.txt"
#	shell:
#		"{SCRIPTPATH}/find_matrix_bestmatch.py {input} 0.4 0.95 > {output}"

rule output_assignment_statistics:
	input:
		"output/matrix_{reads}_{transcripts}_{graph}_all.txt",
		#"output/matrix_{reads}_{transcripts}_{graph}_bestmatch.txt",
		#"output/matrix_{reads}_{transcripts}_{graph}_unambiguous.txt",
		#read_stdout = "tmp/aligner_stdout_{reads}_{graph}.txt",
		#transcript_stdout = "tmp/aligner_stdout.txt"
	output:
		"output/matrixstats_{reads}_{transcripts}_{graph}.txt"
	params:
		allmatrix = "output/matrix_{reads}_{transcripts}_{graph}_all.txt",
		#bestmatchmatrix = "output/matrix_{reads}_{transcripts}_{graph}_bestmatch.txt",
		#unambiguousmatrix = "output/matrix_{reads}_{transcripts}_{graph}_unambiguous.txt"
	run:
		shell("echo 'Number of reads considered:' >> {output}")
		#shell("grep 'Reads with an alignment:' < {input.read_stdout} | cut -d ' ' -f 5 >> {output}")
		shell("echo 'Number of transcripts considered:' >> {output}")
		#shell("grep 'Output end-to-end alignments:' < {input.transcript_stdout} | cut -d ' ' -f 3 >> {output}")
		shell("echo 'Number of reads which overlap a transcript:' >> {output}"),
		shell("cut -f 1 {params.allmatrix} | sort | uniq | wc -l >> {output}"),
		shell("echo 'Number of transcripts which overlap a read:' >> {output}"),
		shell("cut -f 2 {params.allmatrix} | sort | uniq | wc -l >> {output}"),
		shell("echo 'Number of total read-transcript overlaps:' >> {output}"),
		shell("wc -l < {params.allmatrix} >> {output}"),
		#shell("echo 'Number of bestmatch read-transcript overlaps:' >> {output}"),
		#shell("wc -l < {params.bestmatchmatrix} >> {output}"),
		#shell("echo 'Number of unambiguous overlaps:' >> {output}")
		#shell("wc -l < {params.unambiguousmatrix} >> {output}")

rule generateCountMatrix:
	input:
		matrix = "output/matrix_{reads}_{transcripts}_{graph}_all.txt"
	output:
		"output/CountMatrix_{reads}_{transcripts}_{graph}.txt"
	benchmark:
		"benchmark/generateCount_{reads}_{transcripts}_{graph}.txt"
	shell:
		"python {SCRIPTPATH}/ThreePrime.py -m {input.matrix} >> {output}"
