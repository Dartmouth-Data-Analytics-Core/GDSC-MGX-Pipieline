#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Snakemake module for linking taxonomy to function
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# SET GLOBAL SCOPE PYTHON VARIABLES (EXECUTED BEFORE SNAKEMAKE)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
import pandas as pd
import re

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
tax2func = projDir + "tax2Func/"
results = projDir + "results/"
figures = projDir + "Figures/"

#----- Make necessary directories
tmp = config["tmp_dir"]
if not os.path.isdir(tmp):
    os.mkdir(tmp)

#----- Expose terms as a wildcard
def term_to_key(term: str) -> str:
    key = re.sub(r'[^0-9a-zA-Z]+', '_', term.strip().lower())
    return re.sub(r'_+', '_', key).strip('_')

#----- Get terms of interest
terms = config.get("function_term")
if terms is None:
    terms = [config["function_term"]]

#----- Multi-term support
term_keys = [term_to_key(t) for t in terms]
term_key_map = dict(zip(term_keys, terms))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# SNAKEMAKE RULES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#    
rule final:
    input:
        expand(tax2func + "function_subset_{term}_complete.txt", term = term_keys),
        expand(tax2func + "uniref_ID_list_{term}.txt", term = term_keys),
        expand(tax2func + "{sample}_readIDs_{term}.txt", sample = sample_list, term = term_keys),
        expand(tax2func + "{sample}_centrifuger_{term}_matches.txt", sample = sample_list, term = term_keys)
        
#----- Rule to subset gene-families to specific functions
rule subset_function:
    input:
        results + "merged_genefamilies_cpm_named.txt"
    output:
        tax2func + "function_subset_{term}_complete.txt"
    params:
        function_term = lambda wildcards: term_key_map[wildcards.term]
    resources: cpus="10", maxtime="7:00:00", mem_mb="60gb"
    shell: """

    #----- Replace spaces
    term="{params.function_term}"
    term_safe="${{term// /_}}"
    
    #----- Grab function of interest from full results
    grep -i "{params.function_term}" {input[0]} > {output[0]}
    
    """

#----- Rule to grab uniref IDs
rule extract_uniref_IDs:
    input:
        tax2func + "function_subset_{term}_complete.txt"
    output:
        tax2func + "uniref_ID_list_{term}.txt"
    params:
        function_term = lambda wildcards: term_key_map[wildcards.term]
    resources: cpus="10", maxtime="7:00:00", mem_mb="60gb"
    shell: """
    
    #----- Get uniref IDs
    cut -f1 {input[0]} | sed 's/:.*//' | sort -u > {output[0]}
    
    """

#----- Rule to grep for uniref IDs in the bowtie2 temp files
#! This was tested with weird read IDs sequenced from another place? Might need to modify
rule get_readIDs:
    input:
        bt2 = function + "{sample}_bowtie2_aligned.tsv",
        unirefs = tax2func + "uniref_ID_list_{term}.txt"
    output:
        tax2func + "{sample}_readIDs_{term}.txt"
    params:
        sample = lambda wildcards: wildcards.sample
    resources: cpus="10", maxtime="7:00:00", mem_mb="60gb"
    shell: """

    grep -F -f {input.unirefs} {input.bt2} \
    | cut -d'.' -f1 \
    | cut -d':' -f1-7 \
    | sed 's/\.[0-9]\+$//' \
    > {output[0]}
    
    """

#----- Rule to search centrifuger output for read IDs
rule search_centrifuger:
    input:
        readIDs = tax2func + "{sample}_readIDs_{term}.txt",
        cfr = centrifuger + "{sample}_centrifuger.txt",
        quant = centrifuger + "{sample}_centrifuger_quant.tsv"
    output:
        matches = tax2func + "{sample}_centrifuger_{term}_matches.txt"
    params:
        sample = lambda wildcards: wildcards.sample
    resources: cpus="10", maxtime="7:00:00", mem_mb="60gb"
    shell: """
    
    #----- Search centrifuger file
    grep -F -f {input.readIDs} {input.cfr} > tax2Func/{params.sample}_centrifuger_temp.txt &&
    cut -f3 tax2Func/{params.sample}_centrifuger_temp.txt | sort | uniq -c > tax2Func/{params.sample}_centrifuger_matches_temp.txt &&
    rm tax2Func/{params.sample}_centrifuger_temp.txt &&

    #----- Align taxa
    awk 'NR==FNR {{quant[$2]=$1; next}} {{print $0, ($2 in quant ? quant[$2] : 0)}}' \
        {input.quant} \
        tax2Func/{params.sample}_centrifuger_matches_temp.txt > tax2Func/{params.sample}_results_no_header.txt &&
    rm tax2Func/{params.sample}_centrifuger_matches_temp.txt &&

    #----- Add header
    (echo -e "Count\\ttaxID\\ttaxa"; sort -k1,1nr tax2Func/{params.sample}_results_no_header.txt) > {output.matches} &&
    rm tax2Func/{params.sample}_results_no_header.txt
    
    """
    


