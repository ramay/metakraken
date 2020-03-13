print(str(snakemake))
#print(snakemake@input$tab)
#print(snakemake@output$fake)

###### Read in 
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(limma)
# Read in metadata file
print(snakemake@config$metadata_file)

metadata<-read.csv(snakemake@config$metadata_file,header = T)

rownames(metadata)<- metadata[,1]


# Read in clade file 

print(snakemake@input$tab)

tab<- read.table(snakemake@input$tab, header = TRUE,sep="\t")
colnames(tab) <- strsplit2(colnames(tab),split = "_")[,1]

print(colnames(tab))

# Read in species file
#sp <- read.table("control_merged_species.txt", header = TRUE, row.names = 1,sep="\t")
#colnames(sp)<-strsplit2(colnames(sp),split = "_")[,1]

# Define variables



variablex=snakemake@config$variablex
variabley=snakemake@config$variabley
output_dir=snakemake@config$outputRscript
topn=snakemake@config$topn




## R script to analyze the output of metaphlan
## MOdified plot_bar function
plot_bar<-function (physeq, x = "Sample", y = "Abundance", fill = NULL, 
                    title = NULL, facet_grid = NULL) 
{
  mdf = psmelt(physeq)
  fill=reorder(mdf[,fill],-mdf[,y])
  x=reorder(mdf[,x],-mdf[,y])
  p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
  #p = p + geom_bar(stat = "identity", position = "stack", color = fill)
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

plot_taxa<-function(ps,rank,variablex,variabley,cols,topn=NULL){
  
  ps_p<-tax_glom(ps,taxrank = rank) 
  print(rank)
    if(!is.null(topn))
  {
    top20 <- names(sort(taxa_sums(ps_p), TRUE))[1:topn]
    ps_p <- prune_taxa(top20, ps_p)
    print("in topn")
  }
  plot_bar(ps_p, x = "Samples", fill = rank) +
    scale_fill_manual(values = cols) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
    ylab("Relative abundance") +
    facet_wrap(as.formula(paste(variablex,"~", variabley)),
               labeller = label_both,scales="free") +
    theme(axis.text.x = element_text(angle = 90, hjust = 0)) +xlab("Samples")+labs(fill=rank)+
    scale_y_continuous(position = "right")
  #return(ps_p)
}

#Taken from https://rdrr.io/github/mmclaren42/metaphlanr/src/R/metaphlan.R
#Modified to save 
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


## Taxa table

taxa <-parse_taxonomy(tab$ID) %>% select(-Clade)

# Seqtable
seqtab<-tab[,2:ncol(tab)]

match(colnames(seqtab) , metadata$Samples)

### Make Phyloseq object

ps <- phyloseq(otu_table(seqtab, taxa_are_rows=TRUE), sample_data(metadata), tax_table(as.matrix(taxa))) %>% 
  subset_taxa(.,!is.na(Species)) %>%
  subset_taxa(.,is.na(Strain))

#Make a overall Kingdom level plot

king <- tax_glom(ps, "Kingdom") %>% 
  transform_sample_counts(., function (x) x / sum(x)) %>% 
  plot_taxa(.,"Kingdom",variablex,variabley,cols=cols)

ggsave( king,filename = paste0(output_dir,"/kingdom_barplots.pdf"))

Phylum <-tax_glom(ps, "Phylum") %>% 
  transform_sample_counts(., function (x) x / sum(x)) %>% 
  plot_taxa(.,"Phylum",variablex,variabley,cols=cols)

ggsave(Phylum,filename = paste0(output_dir,"/Phylum_barplots.pdf"))

Family <-tax_glom(ps, "Family") %>% 
  transform_sample_counts(., function (x) x / sum(x)) %>% 
  plot_taxa(.,"Family",variablex,variabley,cols=cols)

ggsave(Family,filename = paste0(output_dir,"/Family_barplots.pdf"))

Genus <-tax_glom(ps, "Genus") %>% 
  transform_sample_counts(., function (x) x / sum(x)) %>% 
  plot_taxa(.,"Genus",variablex,variabley,cols=cols,topn)

ggsave(Genus,filename = paste0(output_dir,"/Genus_barplots_top",topn,".pdf"))

Species <-tax_glom(ps, "Species") %>% 
  transform_sample_counts(., function (x) x / sum(x)) %>% 
  plot_taxa(.,"Species",variablex,variabley,cols=cols,topn)
ggsave(Species,filename = paste0(output_dir,"/Species_barplots_top",topn,".pdf"))

top_sp <-tax_glom(ps, "Species")%>% 
  transform_sample_counts(., function (x) x / sum(x))
write.csv(psmelt(top_sp),file = "species_melted.csv",row.names = F)

