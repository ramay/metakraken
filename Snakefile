# ************************************
# * Snakefile for metaphlan pipeline *
# ************************************

# **** Variables ****

configfile: "config.yaml"

# **** Imports ****

import pandas as pd
SAMPLES = pd.read_csv(config["list_files"], header = None)
SAMPLES = SAMPLES[0].tolist()

# **** Rules ****

rule all:
    input:
        config["outputRscript"] + "/species_melted.csv"

rule kraken2:
    input:
        r1 = config["path"]+"{sample}"+config["for"],
        r2 = config["path"]+"{sample}"+config["rev"]
    output:
         kr = "output/kraken2/{sample}_kreport.tab",
#        pr = "output/metaphlan2/{sample}_profile.tab"
    params:
        db = config["kraken_db"],
	threads=28
    log: "output/logs/{sample}_kraken.log"
                
    conda: "utils/envs/kraken2_bracken_final.yaml"
    shell:
            "kraken2 --use-names --db {params.db} --threads {params.threads} --report {output.kr} --paired {input.r1} {input.r2} >> {log} "





rule bracken:
    input:
        r1 = "output/kraken2/{sample}_kreport.tab"
    output:
        br1 = "output/bracken/{sample}_breport",
        br2 = "output/kraken2/{sample}_kreport_bracken.tab"
    params:
        db = config["kraken_db"],
        level = config["level"],
	readlen=config["readlen"],
	threshold=config["threshold"]
    conda: "utils/envs/kraken2_bracken_final.yaml"
    shell:
          "bracken -d {params.db} -i {input.r1} -l {params.level} -t {params.threshold} -r {params.readlen} -o {output.br1}"


rule kreport2mpa:
    input:
        "output/kraken2/{sample}_kreport_bracken.tab"
    output:
        "output/kreport2mpa/{sample}_mpa.tab"
    conda:"utils/envs/kraken2_bracken_final.yaml"
    shell:
          "kreport2mpa.py -r {input} -o {output}"



rule norm_mpa:
    input:
        "output/kreport2mpa/{sample}_mpa.tab"
    output:
        "output/kreport2mpa_norm/{sample}_mpa_norm.tab"
    shell:
        """
        sum=$(grep -vP "\\|" {input} | cut -f 2 | awk '{{sum += $1}} END {{printf ("%.2f\\n", sum/100)}}')
        awk -v sum="$sum" 'BEGIN {{FS="\\t"}} {{OFS="\\t"}} {{print $1,$2/sum}}' {input} > {output}
        """


rule merge_mpa:
    input:
        expand("output/kreport2mpa_norm/{sample}_mpa_norm.tab", sample=SAMPLES)
    output:
        merge = "output/kreport2mpa_norm/merge_metaphlan.txt",
        merge_species = "output/kreport2mpa_norm/merge_metaphlan_species.txt"
    conda:"utils/envs/metaphlan2_env.yaml"
    shell:
        """
        merge_metaphlan_tables.py {input} >  {output.merge}
        sed -i 's/d_/k_/g' {output.merge}; sed -i 's/_/__/g' {output.merge}
        grep -E "(s__)|(^ID)"  {output.merge} | grep -v "t__" | sed 's/^.*s__//g' >  {output.merge_species}
        """

rule make_plots:
    input:
       tab="output/kreport2mpa_norm/merge_metaphlan.txt",
       sp="output/kreport2mpa_norm/merge_metaphlan_species.txt"
    output:
        res=config["outputRscript"] + "/species_melted.csv"
    conda:"utils/envs/rscript.yaml"
    script:
        "utils/scripts/plot_profile.R"


rule heatmap:
    input: "output/kreport2mpa_norm/merged_abundance_table.txt"
    output: "output/plots/abundance_heatmap_species.png"
    conda: "utils/envs/hclust_env.yaml"
    shell:
            """
            hclust2.py -i output/merged_abundance_table_species.txt -o output/abundance_heatmap_species.png --ftop 25 --f_dist_f braycurtis --s_dist_f braycurtis --cell_aspect_ratio 0.5 -l --flabel_size 6 --slabel_size 6 --max_flabel_len 100 --max_slabel_len 100 --minv 0.1 --dpi 300
            """
