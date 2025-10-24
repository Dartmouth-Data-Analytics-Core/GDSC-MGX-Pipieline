#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# GDSC-MGX Pipeline (BioBakery Implementation)
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#≠≠≠≠≠ To Do ≠≠≠≠≠#
# Method for extracting uniprot IDs for groups of interest and then searching in Centrifuger output
# Need some way of memory management for temporary humann bowtie tsvs once queried (delete?)


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
centrifuger = projDir + "centrifuger/"
function = projDir + "function/"
results = projDir + "results/"
figures = projDir + "Figures/"

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
        #----- KneadData outputs
        expand(filtered + "{sample}.fastq.gz", sample = sample_list),
        figures + "KneadData_Diagnostics.png",

        #----- Taxonomy outputs
        expand(taxonomy + "{sample}-profiled_metagenome.txt", sample = sample_list),
        expand(centrifuger + "{sample}_centrifuger.txt", sample = sample_list),
        expand(centrifuger + "{sample}_centrifugert_quant.tsv", sample = sample_list),
        expand(centrifuger + "{sample}_centrifguer_kreport.tsv", sample = sample_list),
        results + "merged_abundance_table.txt",
        results + "merged_genus_abundance_table.txt",
        results + "merged_species_abundance_table.txt",

        #----- Function outputs
        expand(function + "{sample}-humann.log", sample = sample_list),
        expand(function + "{sample}_genefamilies.tsv", sample = sample_list),
        expand(function + "{sample}_pathcoverage.tsv", sample = sample_list),
        expand(function + "{sample}_pathabundance.tsv", sample = sample_list),
        expand(function + "{sample}_bowtie2_aligned.tsv", sample = sample_list),
        function + "merged_genefamilies.txt",
        function + "merged_pathabundances.txt",
        function + "merged_pathcoverages.txt",
        function + "merged_genefamilies_cpm.txt",
        function + "merged_pathabundances_cpm.txt",
        function + "merged_pathcoverages_cpm.txt",

        #----- Functional results
        results + "merged_genefamilies_cpm_named.txt",
        results + "merged_pathabundances_cpm_named.txt",
        results + "merged_pathcoverages_cpm_named.txt",
        results + "merged_taxonomic_profiles.tsv"
    output:
        results + "Species_Shannon_Diversity_Per_Sample.png",
        results + "Genus_Shannon_Diversity_Per_Sample.png"
    params:
        sample_csv = config["sample_csv"]
    conda: "r-plotting"
    resources: cpus="10", maxtime="7:00:00", mem_mb="60gb"
    shell: """
    
        #----- Run Alpha diversity
        Rscript scripts/plot_alpha.R {params.sample_csv}

    """


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

#----- Rule to plot kneadData diagnostics
rule plot_filtering:
    input:
        expand(filtered + "{sample}.fastq.gz", sample = sample_list)    
    output:
        figures + "KneadData_Diagnostics.png"
    conda: "r-plotting"
    resources: cpus="40", maxtime="20:00:00", mem_mb="60gb"
    shell: """
    
        #----- Run plotting script
        Rscript scripts/plot_kneadData.R

    """

#----- Rule to run centrifuger
rule centrifuger:
    output:
        tax_txt = centrifuger + "{sample}_centrifuger.txt",
    params:
        sample = lambda wildcards: wildcards.sample,
        R1 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_1"],
        R2 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_2"] if config["layout"]=="paired" else "None",
        layout = config["layout"],
        cfrPath = config["centrifuger"],
        cfrIdx = config["cfr_idx"],
    conda: "centrifuger"
    threads: 8
    resources: cpus="40", maxtime="20:00:00", mem_mb="60gb"
    shell: """

        echo {params.sample}
        if [ {params.layout} == "single" ]
        then
            #== Run in single-end mode
            {params.cfrPath} \
                -1 {params.R1} \
                -x {params.cfrIdx} \
                -t {threads} > {output.tax_txt}

        else

            #== Run in paired-end mode
            {params.cfrPath} \
                -1 {params.R1} \
                -2 {params.R2} \
                -x {params.cfrIdx} \
                -t {threads} > {output.tax_txt}

        fi

    """

#----- Rule to run centrifuger-quant
rule centrifuger_quant:
    input:
        tax_txt = centrifuger + "{sample}_centrifuger.txt",
    output:
        quant = centrifuger + "{sample}_centrifugert_quant.tsv",
    params:
        cfrIdx = config["cfr_idx"],
        minScore = config["minScore"]
    conda: "centrifuger"
    resources: cpus="40", maxtime="20:00:00", mem_mb="60gb",
    shell: """
    
        #== Run quant
        centrifuger-quant \
            -c {input.tax_txt} \
            -x {params.cfrIdx} > {output.quant}
    
    """

#----- Rule to run centrifuger-kreport
rule centrifuger_kreport:
    input:
        tax_txt = centrifuger + "{sample}_centrifuger.txt",
    output:
        kreport = centrifuger + "{sample}_centrifguer_kreport.tsv"
    params:
        cfr_kreport = config["cfr_kreport"],
        cfrIdx = config["cfr_idx"],
        minScore = config["minScore"]
    conda: "centrifuger"
    resources: cpus="40", maxtime="20:00:00", mem_mb="60gb"
    shell: """
    
        #== Run Kreport
        {params.cfr_kreport} \
            -x {params.cfrIdx} \
            --min-score {params.minScore} \
            {input.tax_txt} > {output.kreport}
    
    """


#----- Rule to assign taxonomy
rule taxonomy:
    input:
        filtered + "{sample}.fastq.gz"
    output: 
        taxonomy + "{sample}-profiled_metagenome.txt"
    params:
        metaphlan_path = config["metaphlan_path"],
        metaphlan_db = config["metaphlan_db"],
        metaphlan_analysis_type = config["metaphlan_analysis_type"],
        mpIndex = config["mpIndex"]
    conda:
        "metaphlan",
    resources: cpus="40", maxtime="20:00:00", mem_mb="60gb",
    shell: """
        {params.metaphlan_path} \
            {input[0]} \
            --input_type fastq \
            -o {output[0]} \
            --bowtie2db {params.metaphlan_db} \
            -x {params.mpIndex} \
            -t {params.metaphlan_analysis_type} \
            --offline
      
"""

#----- Rule to merge abundance tables
rule merge_abundances:
    input: 
        expand(taxonomy + "{sample}-profiled_metagenome.txt", sample = sample_list)
    output:
        results + "merged_abundance_table.txt"
    conda: "metaphlan"
    resources: cpus="40", maxtime="20:00:00", mem_mb="60gb"
    shell: """
    
        #----- Merge metaphlan tables
        merge_metaphlan_tables.py {input} > {output}
    
    """

#----- Rule to split into genus and species abundances
rule extract_genus_species:
    input:
        results + "merged_abundance_table.txt"
    output:
        genus = results + "merged_genus_abundance_table.txt",
        species = results + "merged_species_abundance_table.txt"
    conda: "metaphlan"
    resources: cpus="40", maxtime="20:00:00", mem_mb="60gb"
    shell: """
    
        #----- Split based on taxonomic ranks to species
        grep -E "s__|SRS" {input[0]} \
        | grep -v "t__" \
        | sed "s/^.*|//g" \
        | sed "s/SRS[0-9]*-//g" \
        > {output.species}

        #----- Split based on taxonomic ranks to genus
        grep -E "g__|SRS" {input[0]} \
        | grep -v "t__" \
        | grep -v "s__" \
        | sed "s/^.*|//g" \
        | sed "s/SRS[0-9]*-//g" \
        > {output.genus}
    
    """

#----- Rule to assign function
rule function:
    input: 
        filtered + "{sample}.fastq.gz",
        taxonomy + "{sample}-profiled_metagenome.txt"  
    output: 
        function + "{sample}-humann.log",
        function + "{sample}_genefamilies.tsv",
        function + "{sample}_pathcoverage.tsv",
        function + "{sample}_pathabundance.tsv",
        function + "{sample}_bowtie2_aligned.tsv"
    params:
        humann_path = config["humann_path"],
        sample = lambda wildcards: wildcards.sample,
        mode = config["mode"],
        nt_db = config["nt_db"],
        aa_db = config["aa_db"]
    conda:
        "humann3",
    resources: cpus="40", maxtime="25:00:00", mem_mb="60gb",
    shell: """
		{params.humann_path} -i {input[0]} \
			-o {function} \
            --output-basename {params.sample} \
			--threads 16 --o-log {output[0]} \
			--taxonomic-profile {input[1]} \
			--search-mode {params.mode} \
			--nucleotide-database {params.nt_db} \
			--protein-database {params.aa_db} &&
        
        #----- Move temp files
        mv function/{params.sample}_humann_temp/{params.sample}_bowtie2_aligned.tsv function/{params.sample}_bowtie2_aligned.tsv &&

        #----- Remove temp folders
        rm -r function/{params.sample}_humann_temp
   """

#----- Rule to join functional tables
rule merge_function:
    input:
        expand(function + "{sample}_genefamilies.tsv", sample = sample_list),
        expand(function + "{sample}_pathabundance.tsv", sample = sample_list),
        expand(function + "{sample}_pathcoverage.tsv", sample = sample_list)
    output:
        gfs = function + "merged_genefamilies.txt",
        pas = function + "merged_pathabundances.txt",
        pcs = function + "merged_pathcoverages.txt"
    conda: "humann3"
    resources: cpus="40", maxtime="25:00:00", mem_mb="60gb",
    shell: """
    
        #----- Merge humann tables
        humann_join_tables --input function --output {output.gfs} --file_name genefamilies &&
        humann_join_tables --input function --output {output.pas} --file_name pathabundance &&
        humann_join_tables --input function --output {output.pcs} --file_name pathcoverage

    """

#----- Rule to normalize tables
rule normalize_function:
    input:
        merged_gfs = function + "merged_genefamilies.txt",
        merged_pas = function + "merged_pathabundances.txt",
        merged_pcs = function + "merged_pathcoverages.txt"
    output:
        cpm_gfs = function + "merged_genefamilies_cpm.txt",
        cpm_pas = function + "merged_pathabundances_cpm.txt",
        cpm_pcs = function + "merged_pathcoverages_cpm.txt"
    params:
        normMethod = config["normMethod"]
    conda: "humann3"
    resources: cpus="40", maxtime="25:00:00", mem_mb="60gb",
    shell: """
    
        #----- Run normalization
        humann_renorm_table --input {input.merged_gfs} --units {params.normMethod} --output {output.cpm_gfs} &&
        humann_renorm_table --input {input.merged_pas} --units {params.normMethod} --output {output.cpm_pas} &&
        humann_renorm_table --input {input.merged_pcs} --units {params.normMethod} --output {output.cpm_pcs}
    
    """
    
#----- Rule to rename tables
rule rename_function:
    input:
        cpm_gfs = function + "merged_genefamilies_cpm.txt",
        cpm_pas = function + "merged_pathabundances_cpm.txt",
        cpm_pcs = function + "merged_pathcoverages_cpm.txt"
    output:
        rn_gfs = results + "merged_genefamilies_cpm_named.txt",
        rn_pas = results + "merged_pathabundances_cpm_named.txt",
        rn_pcs = results + "merged_pathcoverages_cpm_named.txt"
    params:
        names = config["mode"]
    conda: "humann3"
    resources: cpus="40", maxtime="25:00:00", mem_mb="60gb",
    shell: """
    
        #----- Run renaming
        humann_rename_table --input {input.cpm_gfs} --names {params.names} --output {output.rn_gfs} &&
        humann_rename_table --input {input.cpm_pas} --names {params.names} --output {output.rn_pas} &&
        humann_rename_table --input {input.cpm_pcs} --names {params.names} --output {output.rn_pcs}
    
    """

#----- Rule to infer taxonomy from function
rule infer_taxonomy_from_function:
    input:
        rn_gfs = results + "merged_genefamilies_cpm_named.txt",
    output:
        taxProfiles = results + "merged_taxonomic_profiles.tsv"
    params:
        taxons = config["taxons"]
    conda: "humann3"
    resources: cpus="40", maxtime="25:00:00", mem_mb="60gb",
    shell: """
    
        #----- Infer taxonomy
        humann_infer_taxonomy -i {input.rn_gfs} -d {params.taxons} > {output.taxProfiles}
    
    """

