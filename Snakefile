#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# GDSC-MGX Pipeline (BioBakery Implementation)
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# SET GLOBAL SCOPE PYTHON VARIABLES (EXECUTED BEFORE SNAKEMAKE)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
import pandas as pd

#----- Set path to config file
configfile: "config.yaml"

#----- read in sample data and extract sample names
samples_df = pd.read_csv(config["sample_csv"]).set_index("sample_id", drop=False)
sample_list = list(samples_df['sample_id'])

#----- Set directories
projDir = config["proj_dir"] + '/'
raw = projDir + "raw_data/"
filtered = projDir + "filtered_fq/"
taxonomy = projDir + "taxonomy/"
function = projDir + "function/"

#----- Make necessary directories
dir_log = projDir + "log/"
if not os.path.isdir(dir_log):
    os.mkdir(dir_log)
tmp = config["tmp_dir"]
if not os.path.isdir(tmp):
    os.mkdir(tmp)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# SNAKEMAKE RULES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#    
rule all:
    input:
        expand(filtered + "{sample}.fastq.gz", sample = sample_list),
        expand(taxonomy + "{sample}-profiled_metagenome.txt", sample = sample_list),
        expand(function + "{sample}-humann.log", sample = sample_list),
        expand(function + "{sample}_genefamilies.tsv", sample = sample_list),
        expand(function + "{sample}_pathabundance.tsv", sample = sample_list)


#----- Filter host reads with Kneaddata
rule filter:
    output:
        filtered + "{sample}.fastq.gz"
    params:
        sample = lambda wildcards: wildcards.sample,
        R1 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_1"],
        R2 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_2"] if config["layout"]=="paired" else "None",
        layout = config["layout"],
        kneaddata = config["kneaddata_path"],
        bowtie_ref = config["bowtie_ref"],
        trimmomatic_path = config["trimmomatic_path"]
    conda:
        "kneaddata",
    resources: cpus="10", maxtime="7:00:00", mem_mb="60gb",
    threads: 16
    shell: """
        echo {params.sample}
        if [ {params.layout} == "single" ]
        then
            #== Run Kneaddata
            python --version
                {params.kneaddata} --unpaired {params.R1} \
                    -db {params.bowtie_ref} \
                    -t 16 --bypass-trf \
                    --remove-intermediate-output \
                    --output {filtered}/{params.sample} \
                    --cat-final-output \
                    --output-prefix {params.sample} \
                    --trimmomatic {params.trimmomatic_path} \
                
            #----- gzip kneaddata outputs
            pigz {filtered}{params.sample}/{params.sample}.fastq
            mv {filtered}{params.sample}/{params.sample}.fastq.gz {filtered}
            mv {filtered}{params.sample}/{params.sample}*contam*.fastq {tmp}
        
        else

            #== Run Kneaddata
            python --version
                {params.kneaddata} --input1 {params.R1} --input2 {params.R2} \
                    -db {params.bowtie_ref} \
                    -t 16 --bypass-trf \
                    --remove-intermediate-output \
                    --output {filtered}/{params.sample} \
                    --cat-final-output \
                    --output-prefix {params.sample} \
                    --trimmomatic {params.trimmomatic_path} \
                
            #----- gzip kneaddata outputs
            pigz {filtered}{params.sample}/{params.sample}.fastq
            mv {filtered}{params.sample}/{params.sample}.fastq.gz {filtered}
            mv {filtered}{params.sample}/{params.sample}*contam*.fastq {tmp}
        fi

"""

rule taxonomy:
    input:
        filtered + "{sample}.fastq.gz"
    output: 
        taxonomy + "{sample}-profiled_metagenome.txt"
    params:
        metaphlan_path = config["metaphlan_path"],
        metaphlan_db = config["metaphlan_db"]
    conda:
        "metaphlan",
    resources: cpus="40", maxtime="20:00:00", mem_mb="20gb",
    shell: """
        {params.metaphlan_path} \
            {input[0]} \
            --input_type fastq \
            -o {output[0]} \
            --bowtie2db {params.metaphlan_db}
      
"""

rule function:
    input: 
        filtered + "{sample}.fastq.gz",
        taxonomy + "{sample}-profiled_metagenome.txt"  
    output: 
        function + "{sample}-humann.log",
        function + "{sample}_genefamilies.tsv",
        function + "{sample}_pathabundance.tsv"
    params:
        humann_path = config["humann_path"],
        sample = lambda wildcards: wildcards.sample,
        mode = config["mode"],
        nt_db = config["nt_db"],
        aa_db = config["aa_db"]
    conda:
        "humann3",
    resources: cpus="40", maxtime="25:00:00", mem_mb="40gb",
    shell: """
		{params.humann_path} -i {input[0]} \
			-o {function} \
            --output-basename {params.sample} \
			--threads 16 --o-log {output[0]} \
			--taxonomic-profile {input[1]} \
			--search-mode {params.mode} \
			--nucleotide-database {params.nt_db} \
			--remove-temp-output \
			--protein-database {params.aa_db}
   """


