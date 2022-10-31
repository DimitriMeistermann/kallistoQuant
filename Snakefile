import glob, os, sys, json, pysam, shutil


WORKING_DIR = os.path.dirname(workflow.snakefile)
if config=={}:
	print("Default config file loaded, from " + WORKING_DIR + "/config.json")
	configfile: WORKING_DIR+"/config.json"

## creation of the logs subdirectory
if not os.path.exists(WORKING_DIR+"/log"):
	os.mkdir(WORKING_DIR+"/log")

#put all config variable as variable in the snakefile
for configVar in config:
	if isinstance(config[configVar], str): exec(configVar+"= '"+config[configVar]+"'")
	else: exec(configVar+"="+str(config[configVar]))


## test of the path provided in the config.json file
if not os.path.exists(FASTQ_PATH):
	print("The directory " + FASTQ_PATH + " doesn't exist. Check the field FASTQ_PATH into the config.json file.")
	sys.exit(0)
else:
	## If the path ends by /, the / is suppressed
	if ( FASTQ_PATH[-1:] == "/" ):
		FASTQ_PATH = FASTQ_PATH[:-1]

INPUT_FASTQS = glob.glob(FASTQ_PATH+'/*.fastq.gz')

SAMPLESsplitted = [os.path.basename(f).split(".") for f in INPUT_FASTQS]
SAMPLES=[]

#remove .fastq.gz to get sample names
for s in SAMPLESsplitted:
	SAMPLES.append(".".join(s[0:-2]))

if(OUTPUT_PATH[-1] == "/") : OUTPUT_PATH = OUTPUT_PATH[:-1]

## suppress the .R1. and .R2. elements for paired-end fastq files for the alignement processus in SAMPLES
if IS_PAIRED_END:
	SAMPLES = [itemR2 for itemR2 in SAMPLES if (PAIR_END_FILE_PATTERN+"2") not in itemR2]	
	SAMPLES = [itemR1.replace((PAIR_END_FILE_PATTERN+"1"),'') for itemR1 in SAMPLES]

#OVERHANG calculation
if READ_LENGTH<1:
	with os.popen("gunzip -c " + INPUT_FASTQS[0] + " | sed -n '2p'","r") as inputRead:
		OVERHANG = len(inputRead.readline()) - 2
else:
	OVERHANG = int(READ_LENGTH) - 1

#forward/reverse adaptator accounting
#if no adptator removal ("FORWARD_ADAPTATOR" and "REVERSE_ADAPTATOR" == ""), input of KALLISTO is directly "FASTQ_PATH"
CUTADAPT_PARAMETER="" # -a forward -A reverse for pairedEnd or just -a forward for singleEnd

if FORWARD_ADAPTATOR == "" and  REVERSE_ADAPTATOR =="":
	KALLISTO_FASTQ_FOLDER = FASTQ_PATH
	USE_CUTADAPT = False
else:
	KALLISTO_FASTQ_FOLDER = OUTPUT_PATH + "/FASTQ_CUTADAPT"
	USE_CUTADAPT = True

if IS_PAIRED_END:
	PAIR_SUFFIX = [PAIR_END_FILE_PATTERN+"1",PAIR_END_FILE_PATTERN+"2"]
	if USE_CUTADAPT: CUTADAPT_PARAMETER = "-a " + FORWARD_ADAPTATOR + " -A " + REVERSE_ADAPTATOR
else:
	PAIR_SUFFIX = [""]
	if USE_CUTADAPT: CUTADAPT_PARAMETER = "-a " + FORWARD_ADAPTATOR

if IS_PAIRED_END:  print("Workflow set on paired end mode")
else : print("Workflow set on single end mode")
if USE_CUTADAPT: print("cutadapt will be used on fastq. To disable it, set FORWARD_ADAPTATOR and REVERSE_ADAPTATOR as an empty string in config file")

if os.path.isfile(OUTPUT_PATH+"/log/genomeIsLoaded"): os.remove(OUTPUT_PATH+"/log/genomeIsLoaded")


##############
rule all: 
	input: OUTPUT_PATH+"/results/multiqc_report.html"

rule CUTADAPT:
	input: expand(FASTQ_PATH+"/{{sample}}{pair}.fastq.gz",pair=PAIR_SUFFIX)
	output: expand(OUTPUT_PATH+"/FASTQ_CUTADAPT/{{sample}}{pair}.fastq.gz",pair=PAIR_SUFFIX)
	params:
		minLen = CUTADAPT_MIN_READ_LEN,
		cpu = THREAD_PER_SAMPLE
	log:
		err=OUTPUT_PATH+"/log/CUTADAPT_{sample}.err"
	shell: """
	cutadapt {CUTADAPT_PARAMETER} -j {params.cpu} --minimum-length {params.minLen}  {input} 2> {log.err} | gzip -c > {output}
	"""	

rule FASTQC:
	input: FASTQ_PATH+"/{sample}{pair}.fastq.gz"
	output: multiext(OUTPUT_PATH+"/fastQC/{sample}{pair}_fastqc",".zip",".html")
	params:
		outpath=OUTPUT_PATH+"/fastQC",
		cpu = 1
	log:
		out=OUTPUT_PATH+"/log/FASTQC_{sample}{pair}.out",
		err=OUTPUT_PATH+"/log/FASTQC_{sample}{pair}.err"
	shell: """
	fastqc -o {params.outpath} {input} 1> {log.out} 2> {log.err}
	"""

rule KALLISTO_QUANT:
	input: 
		fastq=expand(KALLISTO_FASTQ_FOLDER+"/{{sample}}{pair}.fastq.gz",pair=PAIR_SUFFIX),
	output:
		directory(OUTPUT_PATH+"/KALLISTO/{sample}")
	params: cpu = THREAD_PER_SAMPLE
	shell: """
	kallisto quant {OTHER_KALLISTO_ARGS} --threads={THREAD_PER_SAMPLE} --index={KALLISTO_INDEX} --output-dir={output} {input}
	"""

rule COUNTS_TABLE:
	input: expand(OUTPUT_PATH+"/KALLISTO/{sample}",sample=SAMPLES)
	output:
		countTable = OUTPUT_PATH+"/results/rawCountsTable.tsv",
		tpmTable = OUTPUT_PATH+"/results/TranscriptPerMillion.tsv"
	params: cpu = 1
	log:
		out=OUTPUT_PATH+"/log/COUNTS_TABLE.out",
		err=OUTPUT_PATH+"/log/COUNTS_TABLE.err"
	shell: "Rscript {WORKING_DIR}/countsTable.R {OUTPUT_PATH}  1> {log.out} 2> {log.err}"


rule MULTIQC:
	input:
		fastqc=expand(OUTPUT_PATH+"/fastQC/{sample}{pair}_fastqc{ext}", sample=SAMPLES,pair=PAIR_SUFFIX,ext=[".zip",".html"]),
		countTable = OUTPUT_PATH+"/results/rawCountsTable.tsv",
		tpmTable = OUTPUT_PATH+"/results/TranscriptPerMillion.tsv"
	output: OUTPUT_PATH+"/results/multiqc_report.html"
	params:
		outpath = OUTPUT_PATH + "/results",
		cpu = 1
	shell: """
	multiqc -f -e general_stats -e tophat -e bowtie2 {OUTPUT_PATH} -o {params.outpath}
	"""
