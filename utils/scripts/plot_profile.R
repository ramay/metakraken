# ************************************
# A script to generate taxonomic profile plots 
# Developed by Hena R. Ramay
# Bioinformatician @ IMC
# March 2020
# ************************************



## ****IMPORTANT*** 
## Make sure you have modified the varibales  variableX & variableFacet
## associated with this Rscript in the config.yaml file.
## This script works with a snakemake file. 
## It receive a S4 object called snakemake from the snakemake
## file with all the needed information
## To see what information is passed on to this script from snakemake 'print(str(snakemake))'



#saveRDS(snakemake,"snakemake.rds")


# Required packages
suppressMessages(library(phyloseq))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(limma))


log <- file(snakemake@log[[1]], open="wt")
sink(log,type = "message")

message("*** Executing plot_profile.R script........")

message("*** Input from snakemake to R script :")

capture.output(str(snakemake),file=log)


# Read in metadata file
metadata<-read.csv(snakemake@config$metadata_file,header = T)

# An assumption that sample names are stored in column 1 of the file.
rownames(metadata)<- metadata[,1]

# Read in clade file generated by kraken2 snakefile or metaphlan2
tab<- read.table(snakemake@input$tab, header = TRUE,sep="\t")
colnames(tab) <- strsplit2(colnames(tab),split = "_")[,1]


# Define variables for plots

variablex=snakemake@config$variableX
variabley=snakemake@config$variableFacet
output_dir=snakemake@config$outputRscript
topn=snakemake@config$topN


## Change these variables if the plots don't look good
txtsize=9
height=8
width=13



## Modified plot_bar function
plot_bar<-function (physeq, x = "Sample", y = "Abundance", fill = NULL, 
                    title = NULL, facet_grid = NULL) 
{
  mdf = psmelt(physeq)
  fill=reorder(mdf[,fill],-mdf[,y])
  x=reorder(mdf[,x],-mdf[,y])
  p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
  p = p + geom_col()
  p = p + theme(axis.text.x = element_text(angle = 90, hjust = 0))
  if (!is.null(facet_grid)) {
    p <- p + facet_grid(facet_grid)
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}


plot_taxa<-function(ps,taxaRank,variablex,variabley,cols,topn=NULL)
{
  
 ps_p<-tax_glom(ps,taxrank = taxaRank) 

 message("*** Generating ",taxaRank," level plots ....") 

   if(!is.null(topn))
  {
    top20 <- names(sort(taxa_sums(ps_p), TRUE))[1:topn]
    ps_p <- prune_taxa(top20, ps_p)
  }
  plot_bar(ps_p, x = "Samples", fill = taxaRank) +
    scale_fill_manual(values = cols) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
    ylab("Relative abundance") +
    facet_wrap(as.formula(paste(variablex,"~", variabley)),
               labeller = label_both,scales="free") +
    theme(axis.text.x = element_text(angle = 90, hjust = 0)) +xlab("Samples")+labs(fill=taxaRank)+
    scale_y_continuous(position = "right")
}


# Taken from https://rdrr.io/github/mmclaren42/metaphlanr/src/R/metaphlan.R
parse_taxonomy <- function (clade, derep = TRUE) {
  # The clade string has the form
  # "k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae|g__Methanosphaera|s__Methanosphaera_stadtmanae|t__GCF_000012545"
  # truncated at whatever the smallest rank is. We want to be able to parse
  # the clade string for any possible smallest rank, from Kingdom to Strain.
  rank_letters <- c("k", "p", "c", "o", "f", "g", "s", "t")
  tax_pattern  <- paste0("(?:", rank_letters, "__(\\w+))?") %>%
    paste(collapse = "\\|?")
  if (derep) {
    tax <- clade %>%
      unique %>%
      stringr::str_match(tax_pattern)
  } else {
    tax <- clade %>%
      stringr::str_match(tax_pattern)
  }
  colnames(tax) <- c("Clade", "Kingdom", "Phylum", "Class", "Order", "Family",
                     "Genus", "Species", "Strain")
  tax %>% as_tibble()
}


# Taxa table
taxa <-parse_taxonomy(tab$ID) %>% select(-Clade)

# Seqtable
seqtab<-tab[,2:ncol(tab)]


# Make Phyloseq object
ps <- phyloseq(otu_table(seqtab, taxa_are_rows=TRUE), 
	sample_data(metadata), tax_table(as.matrix(taxa))) %>% 
  	subset_taxa(.,!is.na(Species)) %>%
  	subset_taxa(.,is.na(Strain))


# Colors for plots
cols=c("#080068","#f76e40","#169989","#db2335","#efcb4a", "#4363d8", "#42d4f4",
       "#911eb4", "#bfef45", "#fabebe", "#9A6324", "#000075", '#3cb44b', '#ffe119',
       '#4363d8', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe','#008080',
       '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000','#ffd8b1',
       '#000075', '#808080', '#000000', "#DC050C", "#FB8072", "#1965B0", "#7BAFDE",
        "#882E72")


# Make plots and save them as pdf

king <- tax_glom(ps, "Kingdom") %>% 
  transform_sample_counts(., function (x) x / sum(x)) %>% 
  plot_taxa(.,"Kingdom",variablex,variabley,cols=cols)

ggsave( king+theme(text = element_text(size=txtsize)),
	filename = paste0(output_dir,"/kingdom_barplots.pdf"),
	width = width,height = height)

Phylum <-tax_glom(ps, "Phylum") %>% 
  transform_sample_counts(., function (x) x / sum(x)) %>% 
  plot_taxa(.,"Phylum",variablex,variabley,cols=cols)

text = element_text(size=9)
ggsave(Phylum+theme(text = element_text(size=txtsize)),
	filename = paste0(output_dir,"/Phylum_barplots.pdf"),
	width = width,height = height)


Family <-tax_glom(ps, "Family") %>% 
  transform_sample_counts(., function (x) x / sum(x)) %>% 
  plot_taxa(.,"Family",variablex,variabley,cols=cols,topn)

ggsave(Family+theme(text = element_text(size=txtsize)),
	filename = paste0(output_dir,"/Family_barplots_top",topn,".pdf"),
	width = width,height = height)


Genus <-tax_glom(ps, "Genus") %>% 
  transform_sample_counts(., function (x) x / sum(x)) %>% 
  plot_taxa(.,"Genus",variablex,variabley,cols=cols,topn)

ggsave(Genus+theme(text = element_text(size=txtsize)),
	filename = paste0(output_dir,"/Genus_barplots_top",topn,".pdf"),
	width = width,height = height)

Species <-tax_glom(ps, "Species") %>% 
  transform_sample_counts(., function (x) x / sum(x)) %>% 
  plot_taxa(.,"Species",variablex,variabley,cols=cols,topn)

ggsave(Species+theme(text = element_text(size=txtsize)),
	filename = paste0(output_dir,"/Species_barplots_top",topn,".pdf"),
	width = width,height = height)


# Save the melted phyloseq object in a file.
top_sp <-tax_glom(ps, "Species")%>% 
  transform_sample_counts(., function (x) x / sum(x)) %>% psmelt()

write.csv(top_sp,file = paste0(output_dir,"/species_melted.csv"),row.names = F,quote=F)
