# ************************************
# Snakefile for Kraken2 and Bracken  pipeline *
# Compiled by Hena R. Ramay
# Adapted from metaplhan pipeline by Alana Schick
# Bioinformatician @ IMC
# March 2020
# ************************************


configfile: "config.yaml"


## Imports

import pandas as pd
SAMPLES = pd.read_csv(config["list_files"], header = None)
SAMPLES = SAMPLES[0].tolist()



rule all:
    input:
#        config["outputRscript"] + "/species_melted.csv",
        config["output_dir"] + "/plots/abundance_heatmap_species.png"


## Runs kraken2 with the specified database

rule kraken2:
    input:
        r1 = config["path"]+"{sample}"+config["for"],
        r2 = config["path"]+"{sample}"+config["rev"]
    output:
         kr = config["output_dir"]+"/kraken2/{sample}_kreport.tab",
#        pr = config["output_dir"]+"/metaphlan2/{sample}_profile.tab"
    params:
        db = config["kraken_db"],
	threads=1
    log: config["output_dir"]+"/logs/{sample}_kraken.log"
                
    conda: "utils/envs/kb.yaml"
    shell:
            "kraken2 --use-names --db {params.db} --threads {params.threads} --report {output.kr} --paired {input.r1} {input.r2} >> {log} "




##  Runs bracken on above generate kraken reports

rule bracken:
    input:
        r1 = config["output_dir"]+"/kraken2/{sample}_kreport.tab"
    output:
        br1 = config["output_dir"]+"/bracken/{sample}_breport",
        br2 = config["output_dir"]+"/kraken2/{sample}_kreport_bracken.tab"
    params:
        db = config["kraken_db"],
        level = config["level"],
	readlen=config["readlen"],
	threshold=config["threshold"]
    conda: "utils/envs/kb.yaml"
    log: config["output_dir"]+"/logs/{sample}_braken.log"
    shell:
          "bracken -d {params.db} -i {input.r1} -l {params.level} -t {params.threshold} -r {params.readlen} -o {output.br1} >> {log}"


## Convert braken reprot to Metaphlan style output
rule kreport2mpa:
    input:
        config["output_dir"]+"/kraken2/{sample}_kreport_bracken.tab"
    output:
        config["output_dir"]+"/kreport2mpa/{sample}_mpa.tab"
    conda:"utils/envs/kb.yaml"
    shell:
          "kreport2mpa.py -r {input} -o {output}"



## Normalize raw abundances generated by kraken

rule normalize_mpa:
    input:
        config["output_dir"]+"/kreport2mpa/{sample}_mpa.tab"
    output:
        config["output_dir"]+"/kreport2mpa_norm/{sample}_profile.tab"
    shell:
        """
        sum=$(grep -vP "\\|" {input} | cut -f 2 | awk '{{sum += $1}} END {{printf ("%.2f\\n", sum/100)}}')
        awk -v sum="$sum" 'BEGIN {{FS="\\t"}} {{OFS="\\t"}} {{print $1,$2/sum}}' {input} > {output}
        """


## Merge all sample normalize results into one file 
## and also generate a species files

rule merge_mpa:
    input:
        expand(config["output_dir"]+"/kreport2mpa_norm/{sample}_profile.tab", sample=SAMPLES)
    output:
        merge = config["output_dir"]+"/kreport2mpa_norm/merged_metakraken_abundance_table.txt",
        merge_species = config["output_dir"]+"/kreport2mpa_norm/merged_metakraken_abundance_species.txt"
    conda:"utils/envs/metaphlan2_env.yaml"
    shell:
        """
        merge_metaphlan_tables.py {input} >  {output.merge}
        sed -i 's/d_/k_/g' {output.merge}; sed -i 's/_/__/g' {output.merge}
        grep  -a -E "(s__)|(^ID)"  {output.merge} | grep -v "t__" | sed 's/^.*s__//g' >  {output.merge_species}

        """

## Rscript to generate basic plots

rule make_plots:
    input:
       tab=config["output_dir"]+"/kreport2mpa_norm/merged_metakraken_abundance_table.txt",
       sp=config["output_dir"]+"/kreport2mpa_norm/merged_metakraken_abundance_species.txt"
    output:
        res=config["output_dir"] + "/plots/species_melted.csv"
    conda:"utils/envs/rscript.yaml"
    log: config["output_dir"]+"/logs/make_plots.log"
    script:
        "utils/scripts/plot_profile.R"



## clustering heatmap

rule heatmap:
    input: config["output_dir"]+"/kreport2mpa_norm/merged_metakraken_abundance_species.txt"
    output: config["output_dir"] + "/plots/abundance_heatmap_species.png"
    params:
       topn = config["topN"],
    conda: "utils/envs/hclust_env.yaml"
    shell:
            """
            hclust2.py -i {input} -o {output} --ftop {params.topn}  --f_dist_f braycurtis --s_dist_f braycurtis --cell_aspect_ratio 0.5 -l --flabel_size 4 --slabel_size 4 --max_flabel_len 100 --max_slabel_len 100 --minv 0.1 --dpi 300
            """
