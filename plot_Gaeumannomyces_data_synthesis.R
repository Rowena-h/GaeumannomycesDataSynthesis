library(tidyverse)      #2.0.0
library(ape)            #5.8-1
library(ggtree)         #3.99.2
library(aplot)          #0.2.9
library(caper)          #1.0.3
library(cowplot)        #1.1.3
library(ggh4x)          #0.2.8
library(ggplotify)      #0.1.2
library(ggpp)           #0.5.9
library(ggpubr)         #0.6.0
library(ggrepel)        #0.9.5
library(ggsankey)       #0.0.99999
library(rnaturalearth)  #1.0.1
library(sf)             #1.0-16
library(tgutil)         #0.1.15
library(treeio)         #1.29.1

#Set separate legends for composite plots
options("aplot_guides"="keep")

dir <- "W:/GaeumannomycesMetaanalysis/"

#Read in sequence, sample and SH metadata
samples <- read.csv(paste0(dir, "data/sample_list.txt"), sep="\t")
seqs <- read.csv(paste0(dir, "data/seqs.list"), sep="\t", header=FALSE)
sh.map <- read.csv(paste0(dir, "data/sh_list.txt"), sep="\t")

#Read in map data
world.map <- ne_countries(type="map_units", returnclass="sf")


################################################################################
############################### Taxonomy sankey ################################
################################################################################

taxonomy.df <- data.frame(
  orig=c("G. graminis var. tritici",
         "G. graminis var. avenae",
         "G. graminis var. graminis",
         "G. graminis var. graminis",
         "G. graminis var. graminis",
         "G. graminis var. graminis",
         "G. graminis var. graminis",
         "G. graminis var. graminis",
         "G. graminis var. graminis",
         "G. graminis var. graminis",
         "G. graminis var. graminis",
         "G. graminis var. graminis",
         "G. graminis var. graminis",
         "G. graminis var. graminis"),
  new=c("G. tritici",
        "G. avenae",
        "G. graminis",
        "G. graminicola",
        "G. ellisiorum",
        "G. setariicola",
        "G. californicus",
        "G. oryzinus",
        "G. oryzicola",
        "G. floridanus",
        "G. glycinicola",
        "G. fusiformis",
        "G. arxii",
        "G. australiensis")
) %>%
  make_long(orig, new)

taxonomy.df$node <- factor(taxonomy.df$node, levels=c(
  "G. tritici", "G. setariicola", "G. oryzinus", "G. oryzicola", "G. graminis var. tritici", 
  "G. graminis var. graminis", "G. graminis var. avenae", "G. graminis", "G. graminicola", 
  "G. glycinicola", "G. fusiformis", "G. floridanus", "G. ellisiorum", "G. californicus", 
  "G. australiensis", "G. arxii", "G. avenae"
))
taxonomy.df$next_node <- factor(taxonomy.df$next_node, levels=c(
  "G. tritici", "G. setariicola", "G. oryzinus", "G. oryzicola", "G. graminis var. tritici", 
  "G. graminis var. graminis", "G. graminis var. avenae", "G. graminis", "G. graminicola", 
  "G. glycinicola", "G. fusiformis", "G. floridanus", "G. ellisiorum", "G. californicus", 
  "G. australiensis", "G. arxii", "G. avenae"
))

#Plot sankey
gg.taxonomy.sankey <- ggplot(taxonomy.df,
                             aes(x=x, 
                                 next_x=next_x, 
                                 node=node, 
                                 next_node=next_node,
                                 label=node,
                                 fill=node)) +
  geom_sankey(flow.alpha=0.5, colour=NA) +
  geom_sankey_text(aes(
    # Shift labels conditional on position
    x = stage(x,
              after_stat = x + .1 *
                dplyr::case_when(
                  x == 1 ~ -1,
                  x == 2 ~ 1,
                  .default = 0
                )
    ),
    # Align labels conditional on position
    hjust = dplyr::case_when(
      x == "orig" ~ 1,
      x == "new" ~ 0,
      .default = .5
    )
  ),
  size=2.5, fontface="italic") +
  scale_x_discrete(labels=c("Previous\ntaxonomy", "Current\ntaxonomy"),
                   position="top") +
  scale_fill_manual(values=c("#777777", "grey", "grey90", "#777777",
                             "grey", "grey", "grey", "grey",
                             "grey", "grey", "grey", "grey", "grey",
                             "grey", "grey", "grey", "grey90",
                             "grey", "grey")) +
  coord_cartesian(clip="off") +
  theme_void() +
  theme(legend.position="none",
        plot.margin=margin(5.5, 20, 5.5, 50),
        axis.text.x.top=element_text(size=7, face="bold", margin=margin(b=0))) +
  ggpreview(width=3, height=2)

#Write taxonomy sankey to file
#pdf(paste0(dir, "taxonomy_sankey-", Sys.Date(), ".pdf"), width=3, height=2)
gg.taxonomy.sankey
#dev.off()


################################################################################
######################## EPA results on reference tree #########################
################################################################################

#Read in metadata
metadata <- read.csv(paste0(dir, "data/metadata.csv"))

#Format tip labels
metadata <- metadata %>%
  mutate(new.label=paste0('italic("', name, '")'),
         new.label=ifelse(grepl("sp\\.", new.label),
                          sub(' sp\\.")', '"), " sp\\."', new.label),
                          new.label),
         new.label=paste0(new.label, ', " ', strain, '"'),
         new.label=paste0('paste(', new.label, ')'))

#Read in EPA placement results
jplace <- read.jplace(paste0(dir, "epang/its/epa_result.jplace"))
#Extract tree
tree <- jplace@phylo

#Truncate excessively long branch
shortened.edge <- tree$edge[which.max(tree$edge.length), 2]
tree$edge.length[which.max(tree$edge.length)] <- 
  tree$edge.length[which.max(tree$edge.length)] / 2

#Plot base tree
gg.tree <- ggtree(tree, linetype=NA) %<+% metadata +
  xlim(0, 0.1)

#Make dataframe of clade nodes
clades.df <- data.frame(clade=gg.tree$data %>%
                          filter(clade != "outgroup" & !is.na(clade) & isTip == "TRUE") %>%
                          arrange(y) %>%
                          pull(clade) %>% 
                          unique(),
                        node=NA)

#Find the most recent common ancestor for each clade
for (i in 1:length(clades.df$clade)) {
  clades.df$node[i] <- 
    MRCA(gg.tree,
         metadata$tip[metadata$clade == clades.df$clade[i]])
}

#Make alternated coding for highlights on tree
clades.df$box <- rep(c(0,1), length.out=length(clades.df$clade))

#Add clade highlights
gg.tree2 <- gg.tree +
  geom_highlight(data=clades.df, 
                 aes(node=node, fill=as.factor(box)),
                 alpha=1, extend=0.1,
                 show.legend=FALSE) +
  geom_cladelab(data=clades.df,
                mapping=aes(node=node, label=clade),
                fontsize=1.8,
                fontface="italic",
                barsize=0.5,
                align=TRUE,
                offset=0.01,
                offset.text=0.001) +
  scale_fill_manual(values=c("#F5F5F5", "#ECECEC")) +
  geom_tree(aes(linetype=ifelse(as.numeric(label) < 70, "insig", NA)),
            linewidth=0.1,
            show.legend=FALSE) +
  scale_linetype_manual(values="dashed", 
                        na.value="solid")

#Get node coordinates for placement values
placements <- jplace@data %>%
  left_join(gg.tree$data, by="node") %>%
  mutate(nplace=as.numeric(sub(0, NA, nplace)))

#Add raw EPA placements on tree
gg.tree3 <- gg.tree2 +
  geom_point(data=placements,
             aes(x=x, y=y, size=nplace),
             alpha=0.3) +
  geom_label2(aes(x=branch, subset=node == shortened.edge),
              label="//",
              size=2,
              label.padding=unit(0, "pt"),
              label.size=0) +
  scale_size_continuous(breaks=c(1, 100, 1000, 3000)) +
  guides(size=guide_legend(title="Number of EPA\nplacements",
                           title.position="top")) +
  theme(legend.position=c(0.2, 0.7),
        legend.direction="vertical",
        legend.text=element_text(size=5),
        legend.title=element_text(size=6, face="bold"))

#Filter placements for best likelihood weight ratio
placements.filt <- jplace@placements %>%
  group_by(name) %>% 
  filter(like_weight_ratio == max(like_weight_ratio),
         !name %in% jplace@phylo$tip.label) %>%
  arrange(name, like_weight_ratio)

#If there is more than one best likelihood, get the MRCA of the placements
unique.seqs <- unique(placements.filt$name)
seq.mrca <- list()

for (i in 1:length(unique.seqs)){
  
  seq.placements <- placements.filt[grep(unique.seqs[i], placements.filt$name, fixed=TRUE),]
  
  if (length(unique(seq.placements$node)) > 1) {
    
    seq.mrca[unique(seq.placements$name)] <- getMRCA(tree, unique(seq.placements$node))
    
  }
  
}

#Summarise node assignment for sequence variants
seq.node <- placements.filt %>% 
  group_by(name) %>% 
  mutate(num=n_distinct(node)) %>% 
  filter(num == 1) %>%
  dplyr::select(name, node) %>%
  bind_rows(data.frame(name=names(seq.mrca),
                       node=unlist(seq.mrca))) %>%
  mutate(region=ifelse(grepl("ITS1", name), "nplace.ITS1", "nplace.ITS2"),
         clade=NA) %>%
  ungroup()

#Add clade
for (i in 1:length(clades.df$clade)){
  
  members <- clade.members(clades.df$node[i], tree, include.nodes=TRUE)
  
  seq.node$clade[seq.node$node %in% unlist(members)] <- clades.df$clade[i]
  
}

#Summarise number of placements per node
placements.filt2 <- seq.node %>%
  group_by(node, region) %>%
  summarise(nplace=n()) %>%
  pivot_wider(names_from=region, values_from=nplace) %>%
  mutate(nplace=sum(c_across(any_of(c("nplace.ITS1", "nplace.ITS2"))), na.rm=TRUE)) %>%
  left_join(gg.tree$data, by="node")


#Create pie charts to show split of ITS1/2 placements per node
range02 <- function(x){ (x - min(x))/(max(x)-min(x)) * (0.08 - 0.02) + 0.02 }

placements.filt2$plot.size <- range02(placements.filt2$nplace)

ggpie_edited <- function(data, y, fill, color, alpha=1) {
  p <- ggplot(data, aes(x=1, y=value, fill=type)) +
    geom_bar(stat='identity', alpha=alpha, colour="black", linewidth=0.1) +
    coord_polar(theta='y') + 
    theme_inset() +
    scale_fill_manual(values=c("black", "white"))
  
  return(p)
}

all.nodes <- gather(placements.filt2, type, value, !! 2:3) %>% split(., .$node)

pies <- lapply(all.nodes, function(df) {
  ggpie_edited(data=df, y=~value, fill=~type, color=c("black", "white"), alpha=0.5) 
})

#Add pies
gg.tree.4 <- gg.tree2 +
  geom_label2(aes(x=branch, subset=node == shortened.edge),
              label="//",
              size=2,
              label.padding=unit(0, "pt"),
              label.size=NA) +
  geom_point(data=placements.filt2,
             aes(x=x, y=y, size=nplace),
             colour="black",
             alpha=0.3) +
  scale_size_continuous(breaks=c(1, 100, 1000, 3000),
                        range=c(1, 6)) +
  guides(size=guide_legend(title="Number of EPA\nplacements",
                           title.position="top")) +
  theme(legend.position=c(0.2, 0.7),
        legend.direction="vertical",
        legend.text=element_text(size=5),
        legend.title=element_text(size=6, face="bold"),
        legend.background=element_blank())

for (node in placements.filt2$node) {
  
  #Add filtered placements to tree
  gg.tree.4 <- gg.tree.4 +
    geom_inset(pies[which(as.numeric(names(pies)) == node)],
               width=placements.filt2$plot.size[placements.filt2$node == node],
               height=placements.filt2$plot.size[placements.filt2$node == node])
  
}

#Create custom legend for pies
gg.pie.legend <- ggplot(all.nodes[[6]], aes(x=1, y=value, fill=type)) +
  geom_bar(stat='identity', alpha=0.5, colour="black", linewidth=0.1) +
  geom_text(aes(x=2, label=sub("nplace.", "", type)),
            position=position_stack(vjust=0.5),
            fontface="bold",
            size=1.7) +
  coord_polar(theta='y') + 
  theme_inset() +
  scale_fill_manual(values=c("black", "white"))

pie.legend <- as.grob(gg.pie.legend)

gg.tree.5 <- gg.tree.4 + 
  annotation_custom(pie.legend, xmin=0.025, xmax=0.035, ymin=65, ymax=85) +
  ggpreview(width=4, height=4.5)


#Proportion of sequences placed within genus
jplace@placements %>%
  left_join(gg.tree$data, by="node") %>%
  filter(grepl("Gaeumannomyces", label)) %>%
  pull(name.x) %>%
  unique() %>%
  length() / length(unique(jplace@placements$name))


################################################################################
############################### Summary of seqs ################################
################################################################################

#Add sample and UNITE SH data to sequence variants
seq.samples <- seq.node %>%
  mutate(seq.var=sub("\\|.*", "", name),
         sh=sub("sh_", "", seqs$V3[match(seq.var, seqs$V1)]),
         unite.spec=replace_na(sub("Gaeumannomyces", "G.", sh.map$Species[match(sh, sh.map$SH)]), "Unassigned"),
         id=as.numeric(sub("\\|.*", "", sub(".*SampleID_", "", name)))) %>%
  left_join(samples, by="id") %>%
  mutate(MAP=MAP/1000)

#Check that duplicated sequence variants have the same EPA classification
seq.samples %>% group_by(seq.var) %>% summarise(num=n_distinct(clade)) %>% filter(num > 1)

#Number of unique sequence variants
seqs %>% summarise(n_distinct(V1))
seqs %>% group_by(V4) %>% summarise(n_distinct(V1))

#Number of samples
seq.samples %>% group_by(region) %>% summarise(n_distinct(id))

#Studies for ITS1 sequences
seq.samples %>% filter(region == "nplace.ITS1") %>% pull(continent) %>% unique()

#Proportion of sequence variants assigned to species
seq.samples %>% filter(!is.na(clade)) %>% summarise(n_distinct(seq.var)) / 
  seq.samples %>% summarise(n_distinct(seq.var))

#Proportion of sequence variants in agreement across methods
seq.samples %>% filter(unite.spec == clade) %>% summarise(n_distinct(seq.var)) / 
  seq.samples %>% filter(paste(unite.spec, clade) != "G. graminis G. tritici",
                         unite.spec != "Unassigned") %>% summarise(n_distinct(seq.var))

#Proportion of G. graminis variants allocated to G. tritici
seq.samples %>% filter(unite.spec == "G. graminis", clade == "G. tritici") %>% summarise(n_distinct(seq.var)) / 
  seq.samples %>% filter(unite.spec == "G. graminis") %>% summarise(n_distinct(seq.var))

#Proportion of G. californicus variants allocated to G. australiensis
seq.samples %>% filter(unite.spec == "G. californicus", clade == "G. australiensis") %>% summarise(n_distinct(seq.var)) / 
  seq.samples %>% filter(unite.spec == "G. californicus") %>% summarise(n_distinct(seq.var))

#Proportion of UNITE unassigned variants allocated to species
seq.samples %>% filter(unite.spec == "Unassigned", !is.na(clade)) %>% summarise(n_distinct(seq.var)) / 
  seq.samples %>% filter(unite.spec == "Unassigned") %>% summarise(n_distinct(seq.var))

#Proportion of EPA internal node placements that had a species-level UNITE SH
seq.samples %>% filter(is.na(clade), unite.spec != "Unassigned", unite.spec != "G. sp.") %>% summarise(n_distinct(seq.var)) /
  seq.samples %>% filter(is.na(clade)) %>% summarise(n_distinct(seq.var))

#Most frequently recovered species
seq.samples %>% group_by(clade) %>% summarise(num=n(), sv=n_distinct(seq.var)) %>% arrange(desc(num))
seq.samples %>% group_by(clade, region) %>% summarise(num=paste(n(), ",", n_distinct(seq.var))) %>% pivot_wider(id_cols=clade, values_from=num, names_from=region)


#Arrange for file
seq.samples.file <- seq.samples %>%
  dplyr::select("unite.sh"=sh, unite.spec, "EPA.spec"=clade, name, seq.var, "sample_id"=id, paper, permanent_id, sample_type, latitude, longitude, continent, year_of_sampling_from, year_of_sampling_to, Biome, primers, MAT, MAP, pH, SOC, ITS_total, manipulated, abundances)

#Write to file
#write.table(seq.samples.file, paste0(dir, "globalfungi_gaeumannomyces_EPA_results_", Sys.Date(), ".tsv"), sep="\t", row.names=FALSE, quote=FALSE)


################################################################################
################################ Summary trees #################################
################################################################################

#Outgroups
outgroups <- 
  c("Magnaporthe_rhizophila_isolate_M23", "Magnaporthiopsis_maydis_strain_CBS_662.82A",
    "Magnaporthiopsis_sp._CPC_26038", "Magnaporthe_poae_isolate_M48",
    "Magnaporthiopsis_maydis_strain_CBS_133165", "Gaeumannomyces_incrustans_isolate_M35")

#Make summary tree for species found in samples
tree.ingroup <- drop.tip(tree, (metadata %>% group_by(clade) %>% slice(-1) %>% pull(tip)))
tree.ingroup <- drop.tip(tree.ingroup, outgroups)
tree.ingroup$tip.label <- metadata$clade[match(tree.ingroup$tip.label, metadata$tip)]
tip.col <- data.frame(tip=tree.ingroup$tip.label,
                      seqs=tree.ingroup$tip.label %in% placements.filt2$clade)

#Plot summary tree
gg.summary <- ggtree(tree.ingroup, linewidth=0.3, branch.length="none") %<+% tip.col +
  xlim(0, 18)

#Arrange branches to match species tree
gg.summary +
  geom_nodelab(aes(label=node)) +
  geom_tiplab()

gg.summary2 <- flip(gg.summary, 34, 22) %>% 
  ggtree::rotate(28) %>% 
  ggtree::rotate(38) %>% 
  ggtree::rotate(30)

#Make dataframe of clade nodes
summary.clades.df <- gg.summary2$data %>% 
  filter(isTip) %>% 
  arrange(y) %>% 
  dplyr::select(clade=label, node)

#Make alternated coding for highlights on tree
summary.clades.df$box <- rep(c(0,1), length.out=length(summary.clades.df$clade))

#Add clade highlights
gg.summary3 <- gg.summary2 +
  geom_highlight(data=summary.clades.df, 
                 aes(node=node, fill=as.factor(box)),
                 alpha=1, extend=100,
                 show.legend=FALSE) +
  scale_fill_manual(values=c("#F5F5F5", "#ECECEC"))


################################################################################
################################## ITS trees ###################################
################################################################################

## ENTIRE ITS ##

tree.its <- read.tree(paste0(dir, "raxmlng/gaeumannomyces_ITS.raxml.support"))
tree.its <- root(tree.its, outgroup=outgroups)

#Truncate excessively long branch
shortened.edge.its <- tree.its$edge[which.max(tree.its$edge.length), 2]
tree.its$edge.length[which.max(tree.its$edge.length)] <-
  tree.its$edge.length[which.max(tree.its$edge.length)] / 2

#Plot base tree
gg.tree.its <- ggtree(tree.its, linetype=NA) %<+% metadata +
  xlim(0, 0.2)

#Make dataframe of clade nodes
clades.df <- data.frame(clade=gg.tree.its$data %>%
                          filter(clade != "outgroup" & !is.na(clade) & isTip == "TRUE") %>%
                          arrange(y) %>%
                          pull(clade) %>% 
                          unique(),
                        node=NA)

#Find the most recent common ancestor for each clade
for (i in 1:length(clades.df$clade)) {
  clades.df$node[i] <- 
    MRCA(gg.tree.its,
         metadata$tip[metadata$clade == clades.df$clade[i]])
}

#Make alternated coding for highlights on tree
clades.df$box <- rep(c(0,1), length.out=length(clades.df$clade))

#Add clade highlights
gg.tree.its.2 <- gg.tree.its +
  geom_cladelab(data=clades.df,
                mapping=aes(node=node, label=clade),
                fontsize=1.8,
                fontface="italic",
                barsize=0.5,
                align=TRUE,
                offset=0.07,
                offset.text=0.001) +
  scale_fill_manual(values=c("#F5F5F5", "#ECECEC")) +
  geom_tree(aes(linetype=ifelse(as.numeric(label) < 70, "insig", NA)),
            linewidth=0.1,
            show.legend=FALSE) +
  geom_label2(aes(x=branch, subset=node == shortened.edge.its),
              label="//",
              size=2,
              label.padding=unit(0, "pt"),
              label.size=0) +
  geom_tiplab(aes(label=new.label),
              parse=TRUE,
              size=1.8) +
  annotate(geom="label", x=max(gg.tree.its$data$x) * 0.15, y=max(gg.tree.its$data$y) * 0.8, label="ITS") +
  scale_linetype_manual(values="dashed", 
                        na.value="solid") +
  ggpreview(width=3, height=6)


## ITS1 ##

tree.its1 <- read.tree(paste0(dir, "raxmlng/gaeumannomyces_ITS1.raxml.support"))
tree.its1 <- root(tree.its1, outgroup=outgroups)

#Plot base tree
gg.tree.its1 <- ggtree(tree.its1, linetype=NA) %<+% metadata +
  xlim(0, 0.4)

#Make dataframe of clade nodes
clades.df <- data.frame(clade=gg.tree.its1$data %>%
                          filter(!clade %in% c("outgroup", "G. tritici", "G. avenae", "G. hyphopodioides"),
                                 !is.na(clade),
                                 isTip == "TRUE") %>%
                          arrange(y) %>%
                          pull(clade) %>% 
                          unique(),
                        node=NA)

#Find the most recent common ancestor for each clade
for (i in 1:length(clades.df$clade)) {
  clades.df$node[i] <- 
    MRCA(gg.tree.its1,
         metadata$tip[metadata$clade == clades.df$clade[i]])
}

#Make alternated coding for highlights on tree
clades.df$box <- rep(c(0,1), length.out=length(clades.df$clade))

#Add clade highlights
gg.tree.its1.2 <- gg.tree.its1 +
  geom_cladelab(data=clades.df,
                mapping=aes(node=node, label=clade),
                fontsize=1.8,
                fontface="italic",
                barsize=0.5,
                align=TRUE,
                offset=0.1,
                offset.text=0.001) +
  geom_strip("Gaeumannomyces_tritici_strain_CBS_186.65", "Gaeumannomyces_tritici_strain_CBS_131293",
             label="G. tritici",
             fontsize=1.8,
             fontface="italic",
             offset=0.1,
             offset.text=0.001) +
  geom_strip("Gaeumannomyces_avenae_strain_CPC_26255", "Gaeumannomyces_tritici_strain_CPC_26268",
             label="G. avenae",
             fontsize=1.8,
             fontface="italic",
             offset=0.1,
             offset.text=0.001) +
  geom_strip("Gaeumannomyces_hyphopodioides_strain_CPC_26267", "Gaeumannomyces_hyphopodioides_strain_CBS_541.86",
             label="G. hyphopodioides",
             fontsize=1.8,
             fontface="italic",
             offset=0.1,
             offset.text=0.001) +
  scale_fill_manual(values=c("#F5F5F5", "#ECECEC")) +
  geom_tree(aes(linetype=ifelse(as.numeric(label) < 70, "insig", NA)),
            linewidth=0.1,
            show.legend=FALSE) +
  geom_tiplab(aes(label=new.label),
              parse=TRUE,
              size=1.6) +
  annotate(geom="label", x=max(gg.tree.its1$data$x) * 0.15, y=max(gg.tree.its1$data$y) * 0.8, label="ITS1") +
  scale_linetype_manual(values="dashed", 
                        na.value="solid") +
  ggpreview(width=3, height=6)


## ITS2 ##

tree.its2 <- read.tree(paste0(dir, "raxmlng/gaeumannomyces_ITS2.raxml.support"))
tree.its2 <- root(tree.its2, outgroup=outgroups)

#Truncate excessively long branch
shortened.edge.its2 <- tree.its2$edge[which.max(tree.its2$edge.length), 2]
tree.its2$edge.length[which.max(tree.its2$edge.length)] <-
  tree.its2$edge.length[which.max(tree.its2$edge.length)] / 2

#Plot base tree
gg.tree.its2 <- ggtree(tree.its2, linetype=NA) %<+% metadata +
  xlim(0, 0.3)

#Make dataframe of clade nodes
clades.df <- data.frame(clade=gg.tree.its2$data %>%
                          filter(!clade %in% c("outgroup", "G. oryzinus",
                                               "G. oryzicola", "G. arxii",
                                               "G. fusiformis", "G. graminicola"),
                                 !is.na(clade),
                                 isTip == "TRUE") %>%
                          arrange(y) %>%
                          pull(clade) %>% 
                          unique(),
                        node=NA)

#Find the most recent common ancestor for each clade
for (i in 1:length(clades.df$clade)) {
  clades.df$node[i] <- 
    MRCA(gg.tree.its2,
         metadata$tip[metadata$clade == clades.df$clade[i]])
}

#Make alternated coding for highlights on tree
clades.df$box <- rep(c(0,1), length.out=length(clades.df$clade))

#Add clade highlights
gg.tree.its2.2 <- gg.tree.its2 +
  geom_cladelab(data=clades.df,
                mapping=aes(node=node, label=clade),
                fontsize=1.8,
                fontface="italic",
                barsize=0.5,
                align=TRUE,
                offset=0.1,
                offset.text=0.001) +
  geom_strip("Gaeumannomyces_oryzinus_strain_CPC_26031", "Gaeumannomyces_oryzinus_strain_CPC_26065",
             label="G. oryzinus",
             fontsize=1.8,
             fontface="italic",
             offset=0.1,
             offset.text=0.001) +
  scale_fill_manual(values=c("#F5F5F5", "#ECECEC")) +
  geom_tree(aes(linetype=ifelse(as.numeric(label) < 70, "insig", NA)),
            linewidth=0.1,
            show.legend=FALSE) +
  geom_label2(aes(x=branch, subset=node == shortened.edge.its2),
              label="//",
              size=2,
              label.padding=unit(0, "pt"),
              label.size=0) +
  geom_tiplab(aes(label=new.label),
              parse=TRUE,
              size=1.8) +
  annotate(geom="label", x=max(gg.tree.its2$data$x) * 0.15, y=max(gg.tree.its2$data$y) * 0.8, label="ITS2") +
  scale_linetype_manual(values="dashed", 
                        na.value="solid") +
  ggpreview(width=3, height=6)

#Combine all ITS trees
gg.its.trees <- ggarrange(gg.tree.its.2, gg.tree.its1.2, gg.tree.its2.2, nrow=1)

#Write ITS trees to file
#pdf(paste0(dir, "its_trees-", Sys.Date(), ".pdf"), width=9, height=6)
gg.its.trees
#dev.off()


################################################################################
####################### Published host/lifestyle reports #######################
################################################################################

#Read in literature reports
hosts <- read.csv(paste0(dir, "data/hosts.csv"))

#Format for grid
hosts2 <- hosts %>%
  mutate(value="Y",
         label=ifelse(group == "Crop", common.name, host),
         label=ifelse(group == "Crop",
                      paste0('paste("', label, '")'),
                      paste0('paste(italic("', label, '"))')),
         label=ifelse(grepl("sp\\.", label),
                      sub(' sp."))', '"), " sp.")', label),
                      label),
         label=ifelse(grepl("spp\\.", label),
                      sub(' spp."))', '"), " spp.")', label),
                      label)) %>%
  complete(species, label, group, taxonomy) %>%
  mutate(lifestyle=replace_na(lifestyle, "unspecified"),
         group=paste0('paste("', group, '")'),
         taxonomy=paste0('paste(bolditalic("', taxonomy, '"))')) %>%
  group_by(label, group, taxonomy) %>% 
  filter("Y" %in% value)

# #For common names 
# parse_names <- function(x) {
#   parse(text=ifelse(x %in% c("Ammopiptanthus mongolicus", "Ctenanthe"),
#                     paste0('paste(italic("', x, '"))'),
#                     paste0('paste("', x, '")')))
# }

#Plot grid
gg.host.grid <- ggplot(hosts2, aes(y=species, x=label, fill=common.name, alpha=value)) +
  facet_nested(~taxonomy+group, scales="free", space="free",
               nest_line=element_line(), solo_line=TRUE,
               strip=strip_nested(clip="off"),
               labeller=label_parsed) +
  geom_tile(colour="white", linewidth=0.5) +
  geom_point(data=hosts2,
             aes(shape=lifestyle),
             size=2,
             colour="white") +
  scale_shape_manual(values=c(16, 13, 30),
                     labels=c("non-pathogen", "pathogen", "unspecified")) +
  scale_alpha_manual(values=1, na.value=0, guide="none") +
  scale_fill_manual(breaks=c("Soybean",
                             "Barley", "Maize", "Oat",
                             "Rice", "Rye", "Triticale", "Wheat",
                             "Annual meadow grass", "Blackgrass", "Buffalo grass", "Brome grass",
                             "Canary grass", "Centipede grass", "Cockspur", "Cocksfoot grass", "Common bent grass", "Couch grass",
                             "Foxtail millet", "Yellow foxtail", "Goosegrass", "Hybrid Bermuda grass", "Kikuyu grass", "Red fescue", "Ryegrass",
                             "Timothy", "Tufted hairgrass", "Vetiver grass", "Wheatgrass", "Yorkshire fog",
                             "Malacca ginger", "Ctenanthe", "Wild ginger",
                             "Ammopiptanthus mongolicus", "Prairie milkvetch"),
                    values=c("#44BB99",
                             "#EE8866", "#DDAA33", "#EEDD88",
                             "#77AADD", "#F1932D", "#AAAA00", "#117733",
                             "#BBCC33", "#BBCC33", "#BBCC33", "#BBCC33",
                             "#BBCC33", "#BBCC33", "#BBCC33", "#BBCC33",
                             "#BBCC33", "#BBCC33", "#BBCC33", "#BBCC33",
                             "#BBCC33", "#BBCC33", "#BBCC33", "#BBCC33",
                             "#BBCC33", "#BBCC33", "#BBCC33", "#BBCC33",
                             "#BBCC33", "#BBCC33",
                             "#FFAABB", "#FFAABB", "#FFAABB",
                             "#BA8DB4", "#BA8DB4"),
                    guide="none") +
  scale_x_discrete(labels=function(x) parse(text=x)) +
  labs(shape="Reported lifestyle") +
  theme_void() +
  theme(legend.position=c(-0.14, -0.18),
        legend.direction="vertical",
        legend.title=element_text(size=7, face="bold"),
        legend.text=element_text(size=6),
        legend.key=element_rect(fill="dimgrey", colour=NA),
        legend.key.spacing.y=unit(2, "pt"),
        legend.key.size=unit(10, "pt"),
        strip.placement="outside",
        strip.text=element_text(size=6, margin=margin(2, 0, 2, 0), face="bold"),
        panel.border=element_rect(fill=NA, colour="black", linewidth=0.5),
        axis.text.x=element_text(size=5, angle=90, hjust=1, margin=margin(2, 0, 0, 0)))

#Format summary tree
gg.summary.hosts <- gg.summary3 +
  xlim(0, 20) +
  geom_tiplab(fontface="italic", size=2.5) +
  coord_cartesian(clip="off")

#Combine
gg.hosts <- gg.host.grid %>% insert_left(gg.summary.hosts, width=0.4)
ggpreview(gg.hosts, height=4.3, width=7)

#Write host/lifestyle grid to file
pdf(paste0(dir, "hosts-", Sys.Date(), ".pdf"), width=7, height=4.3)
gg.hosts
dev.off()


################################################################################
############################ Classification sankey #############################
################################################################################

#Format classification changes for sankey
sankey.df <- seq.samples %>%
  mutate(unite.spec=sub("G. sp.", "Unassigned", unite.spec),
         clade=replace_na(clade, "Unassigned")) %>%
  dplyr::select(name, unite.spec, clade) %>%
  make_long(unite.spec, clade) %>%
  mutate(node=factor(node,
                     levels=rev(c("G. tritici", "G. graminis", "G. avenae", "G. arxii",
                                  "G. glycinicola", "G. amomi", "G. caricis",
                                  "G. ellisiorum", "G. graminicola", "G. fusiformis",
                                  "G. floridanus", "G. nanograminis", "G. oryzicola",
                                  "G. oryzinus", "G. australiensis", "G. californicus",
                                  "G. hyphopodioides", "G. radicicola", "G. setariicola",
                                  "G. wongoonoo", "Unassigned"))),
         next_node=factor(next_node, levels=rev(c(rev(gg.summary$data %>%
                                                        filter(isTip) %>% 
                                                        arrange(y) %>% 
                                                        pull(label)),
                                                  "G. caricis", "Unassigned"))))

#Plot sankey
gg.sankey <- ggplot(sankey.df, aes(x=x, 
                                   next_x=next_x, 
                                   node=node, 
                                   next_node=next_node,
                                   fill=node,
                                   label=node)) +
  geom_sankey(flow.alpha=0.5, colour=NA, space=500) +
  geom_sankey_text(aes(
    # Shift labels conditional on position
    x=stage(x,
            after_stat=x + .1 *
              dplyr::case_when(
                x == 1 ~ -1,
                x == 2 ~ 1,
                .default=0
              )
    ),
    # Align labels conditional on position
    hjust=dplyr::case_when(
      x == "unite.spec" ~ 1,
      x == "clade" ~ 0,
      .default=.5
    )
  ),
  space=500, size=2, fontface="italic") +
  scale_x_discrete(labels=c("UNITE SH", "EPA placement"),
                   position="top") +
  scale_fill_manual(values=c("#D1BBD7", "#BA8DB4", "#AA6F9E", "#882E72",
                             "#1965B0", "#5289C7", "#7BAFDE", "grey", "grey", 
                             "grey", "grey", "#CAE0AB", "#F7F056", "#F6C141",
                             "#F1932D", "#777777", "#E8601C", "#DC050C",
                             "grey", "#4EB265", "grey"),
                    breaks=c("G. graminis", "G. avenae", "G. arxii", "G. glycinicola",
                             "G. caricis", "G. ellisorium", "G. graminicola", "G. fusiformis", "G. floridanus"
                             , "G. nanograminis", "G. oryzicola", "G. californicus", "G. hyphopodioides", "G. radicicola",
                             "G. setariicola", "Unassigned", "G. tritici", "G. amomi",
                             "G. oryzinus", "G. australiensis", "G. wongoonoo")) +
  coord_cartesian(clip="off") +
  theme_void() +
  theme(legend.position="none",
        plot.margin=margin(5.5, 10, 5.5, 10),
        axis.text.x.top=element_text(size=6, face="bold", margin=margin(b=-7))) +
  ggpreview(width=3, height=4.5)

#Write EPA classification results to file
#pdf(paste0(dir, "epa_results-", Sys.Date(), ".pdf"), width=7, height=4.5)
plot_grid(gg.tree.5, gg.sankey, rel_widths=c(4, 3), labels="auto")
#dev.off()


################################################################################
############################### Sample metadata ################################
################################################################################

#Format summary tree
gg.summary.samples <- gg.summary3 +
  geom_tiplab(aes(colour=seqs), fontface="italic", size=3, show.legend=FALSE) +
  scale_colour_manual(values=c("grey", "black")) +
  coord_cartesian(clip="off")

#Summarise number of sequences/samples for each species
seq.samples.num <- seq.samples %>%
  filter(!is.na(clade)) %>%
  group_by(clade) %>%
  summarise(num.sequences=n(), num.samples=n_distinct(id)) %>%
  bind_rows(data.frame(clade=tree.ingroup$tip.label[!tree.ingroup$tip.label %in% unique(seq.samples$clade)],
                       num.sequences=0, num.samples=0))

#Filter for unique samples classified to species-level
samples.filt <- seq.samples %>%
  filter(!is.na(clade)) %>%
  distinct(id, clade, .keep_all=TRUE) %>%
  add_row(clade=seq.samples.num$clade[seq.samples.num$num.samples == 0])

#Filter for non-manipulated samples for species with more than 10 sequences
samples.filt2 <- samples.filt %>%
  filter(manipulated == 0)

#Summarise means for each abiotic measurement
samples.sum <- samples.filt %>%
  filter(manipulated == 0) %>%
  dplyr::select(clade, pH, MAT, MAP, SOC) %>%
  group_by(clade) %>%
  summarise(pH=mean(na.omit(pH)), MAT=mean(na.omit(MAT)),
            MAP=mean(na.omit(MAP)), SOC=mean(na.omit(SOC))) %>%
  filter(!is.na(pH))

#Generate theme for plots
theme_aplot <- theme(axis.text.y=element_blank(),
                     axis.title.y=element_blank(),
                     axis.ticks.y=element_blank(),
                     axis.title.x=element_text(size=7),
                     axis.text.x=element_text(size=6),
                     panel.grid=element_blank(),
                     plot.background=element_blank(),
                     panel.background=element_blank(),
                     legend.position="top",
                     legend.direction="horizontal",
                     legend.text=element_text(size=6, margin=margin(0, 0, 0, 2)),
                     legend.title=element_text(size=7, face="bold"),
                     legend.key.spacing.y=unit(0, "pt"),
                     legend.key.size=unit(0.5, "lines"))

#Plot bar graph of sample types
gg.sampletype <- ggplot(samples.filt, aes(y=clade, fill=sample_type)) +
  geom_bar(width=0.8) +
  scale_fill_manual(breaks=c("air", "dust", "water", "sediment", "litter", "shoot", "root",
                             "root + rhizosphere soil", "rhizosphere soil", "soil", "topsoil"),
                    values=c('#99DDFF', '#DDDDDD', '#77AADD', "#663333", "dimgrey" , '#44BB99',
                             '#BBCC33', '#AAAA00','#EEDD88', '#EE8866', '#FFAABB'),
                    na.value=NA) +
  labs(x="Sample types", fill="Sample types") +
  guides(fill=guide_legend(title.position="top",
                           nrow=3)) +
  theme_aplot +
  theme(legend.position=c(3.6, 1.17),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

#Plot bar graph of sample biomes
gg.biome <- ggplot(samples.filt, aes(y=clade, fill=Biome)) +
  geom_bar(width=0.8) +
  scale_fill_manual(breaks=c("anthropogenic", "cropland", "desert", "grassland", "shrubland",
                             "woodland", "forest", "wetland", "aquatic"),
                    values=c('#AA4499', '#CC6677','#DDCC77', '#999933', '#882255',  
                             '#117733', '#44AA99',  '#332288', '#88CCEE'),
                    na.value=NA) +
  labs(x="Biomes", fill="Biomes") +
  guides(fill=guide_legend(title.position="top",
                           nrow=3)) +
  theme_aplot +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

#Plot boxplot of sample pH
gg.pH <- ggplot(samples.filt2, aes(x=clade, y=pH)) +
  geom_violin(linewidth=0.2,
              colour="dimgrey",
              fill="white") +
  geom_boxplot(width=0.1,
               linewidth=0.2,
               outlier.size=0.2,
               outlier.colour="dimgrey",
               colour="dimgrey") +
  geom_point(data=samples.sum,
             shape=4,
             stroke=0.8,
             size=0.4,
             colour="black") +
  geom_line(data=samples.sum,
            aes(group=1),
            linetype="dashed",
            linewidth=0.2) +
  scale_y_continuous(limits=c(4, 9)) +
  coord_flip() +
  theme_aplot +
  theme(panel.grid.major.x=element_line(colour="white"),
        panel.grid.minor.x=element_line(colour="white"))

#Plot boxplot of sample site temperatures
gg.MAT <- ggplot(samples.filt2, aes(x=clade, y=MAT)) +
  geom_violin(linewidth=0.2,
              colour="dimgrey",
              fill="white") +
  geom_boxplot(width=0.1,
               linewidth=0.2,
               outlier.size=0.2,
               outlier.colour="dimgrey",
               colour="dimgrey") +
  geom_point(data=samples.sum,
             shape=4,
             stroke=0.8,
             size=0.4,
             colour="black") +
  geom_line(data=samples.sum,
            aes(group=1),
            linetype="dashed",
            linewidth=0.2) +
  scale_y_continuous(limits=c(-10, 30)) +
  labs(y="Mean annual\ntemperature (°C)") +
  coord_flip() +
  theme_aplot +
  theme(panel.grid.major.x=element_line(colour="white"),
        panel.grid.minor.x=element_line(colour="white"))

#Plot boxplot of sample site precipitation
gg.MAP <- ggplot(samples.filt2, aes(x=clade, y=MAP)) +
  geom_violin(linewidth=0.2,
              colour="dimgrey",
              fill="white") +
  geom_boxplot(width=0.1,
               linewidth=0.2,
               outlier.size=0.2,
               outlier.colour="dimgrey",
               colour="dimgrey") +
  geom_point(data=samples.sum,
             shape=4,
             stroke=0.8,
             size=0.4,
             colour="black") +
  geom_line(data=samples.sum,
            aes(group=1),
            linetype="dashed",
            linewidth=0.2) +
  scale_y_continuous(limits=c(0, 5)) +
  labs(y="Mean annual\nprecipitation (m)") +
  coord_flip() +
  theme_aplot +
  theme(panel.grid.major.x=element_line(colour="white"),
        panel.grid.minor.x=element_line(colour="white"))

#Generate plot spacer
gg.spacer <- ggplot() +
  theme_void()

#Combine all plots
gg.samples <- gg.biome %>%
  insert_left(gg.summary.samples, width=3.2) %>%
  insert_right(gg.sampletype, width=1) %>%
  insert_right(gg.pH, width=1.4) %>%
  insert_right(gg.spacer, width=0.2) %>%
  insert_right(gg.MAT, width=1.4) %>%
  insert_right(gg.spacer, width=0.2) %>%
  insert_right(gg.MAP, width=1.4)

ggpreview(gg.samples, width=6, height=4)

#pdf(paste0(dir, "sample_comp-", Sys.Date(), ".pdf"), width=6, height=4)
gg.samples
#dev.off()

## Interesting samples to check original papers for ##

#Samples from air
samples.filt %>%
  filter(sample_type == "air") %>%
  dplyr::select(clade, paper, continent, Biome, year_of_sampling_from, longitude, latitude) %>%
  group_by(clade, Biome, paper, continent, year_of_sampling_from) %>%
  summarise(n()) %>%
  arrange(paper) %>%
  print(n=30)

#Samples from shoots
samples.filt %>%
  filter(sample_type == "shoot") %>%
  dplyr::select(clade, paper) %>%
  arrange(clade, paper) %>%
  unique() %>%
  print(n=30)

#Samples from water
samples.filt %>%
  filter(sample_type == "water") %>%
  dplyr::select(clade, paper, continent, Biome, year_of_sampling_from, longitude, latitude) %>%
  unique()

#Samples from forest or woodland
samples.filt %>%
  filter(Biome == "forest" | Biome == "woodland") %>%
  dplyr::select(clade, sample_type, paper) %>%
  arrange(clade) %>%
  unique() %>%
  print(n=150)

#Abiotic measurement ranges
samples.filt %>%
  filter(manipulated == 0) %>%
  dplyr::select(clade, pH, MAT, MAP, SOC) %>%
  summarise(mean.pH=mean(na.omit(pH)), min.pH=min(na.omit(pH)), max.pH=max(na.omit(pH)),
            mean.MAT=mean(na.omit(MAT)), min.MAT=min(na.omit(MAT)), max.MAT=max(na.omit(MAT)),
            mean.MAP=mean(na.omit(MAP)), min.MAP=min(na.omit(MAP)), max.MAP=max(na.omit(MAP)))

samples.filt %>%
  filter(manipulated == 0, sample_type != "water|air") %>%
  dplyr::select(clade, pH, MAT, MAP, SOC) %>%
  summarise(mean.SOC=mean(na.omit(SOC)), min.SOC=min(na.omit(SOC)), max.SOC=max(na.omit(SOC)))


################################################################################
########################### All GlobalFungi samples ############################
################################################################################

#Read in all samples
all.samples <- read.csv(paste0(dir, "data/GF5_allsamples.csv"))

#Fix ITS1/ITS2 assignment
all.samples.split <- all.samples %>%
  filter(target_gene == "ITSboth") %>%
  group_by(database_ID) %>%
  slice(rep(1, 2)) %>%
  mutate(target_gene=c("ITS1", "ITS2")) %>%
  dplyr::select(target_gene, everything())

all.samples.comb <- all.samples %>%
  filter(target_gene != "ITSboth") %>%
  bind_rows(all.samples.split)

#Format for barplot
all.samples.bar <- all.samples.comb %>%
  mutate(sample_type=factor(sample_type, 
                            levels=c("air", "dust", "water", "glacial ice debris", "deadwood", "litter",
                                     "coral", "fungal sporocarp", "lichen", "mosses", "shoot",
                                     "rhizosphere soil", "soil", "topsoil", "sediment", "root + rhizosphere soil", "root")),
         Biome=factor(Biome, 
                      levels=c("anthropogenic", "cropland", "grassland", "shrubland",
                               "woodland", "forest", "mangrove", "wetland", "aquatic", "tundra", "desert")))

#Format for line plot
all.sample.values <- all.samples.comb %>%
  filter(manipulated == 0) %>%
  dplyr::select(Biome, pH, MAT, MAP, target_gene) %>%
  mutate(MAP=as.numeric(MAP)/1000) %>%
  gather(measurement, value, -c(Biome, target_gene)) %>%
  mutate(value=as.numeric(sub("NA_", NA, value)))

#Summarise means for each abiotic measurement
all.sample.sum <- all.samples.comb %>%
  filter(manipulated == 0) %>%
  dplyr::select(Biome, pH, MAT, MAP, target_gene) %>%
  mutate(MAP=as.numeric(MAP)/1000,) %>%
  group_by(Biome, target_gene) %>%
  summarise(pH=mean(na.omit(as.numeric(pH))), MAT=mean(na.omit(as.numeric(MAT))),
            MAP=mean(na.omit(as.numeric(MAP)))) %>%
  ungroup() %>%
  gather(measurement, value, -c(Biome, target_gene))

#Plot bar graph of sample biomes 
gg.sample.biome.its1 <- ggplot(all.samples.bar %>% filter(target_gene == "ITS1"), aes(y=Biome, fill=sample_type)) +
  geom_bar(position="fill", width=0.5) +
  scale_fill_manual(breaks=c("air", "dust", "water", "glacial ice debris", "deadwood", "litter",
                             "coral", "fungal sporocarp", "lichen", "mosses", "shoot",
                             "rhizosphere soil", "soil", "topsoil", "sediment", "root + rhizosphere soil", "root"),
                    values=c('#99DDFF', '#99DDFF', '#77AADD', "dimgrey", '#EEDD88', '#EEDD88',
                             '#FFAABB',  '#DDDDDD', '#DDDDDD', '#BBCC33', '#BBCC33',
                             '#EE8866', '#EE8866', '#EE8866', '#EE8866', '#EE8866', '#EE8866')) +
  labs(fill=NULL, x="Sample types") +
  scale_x_continuous(position="top") +
  scale_y_discrete(limits=rev) +
  guides(fill=guide_legend(nrow=3)) +
  theme(axis.title.x.top=element_text(size=9, vjust=-6),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        plot.background=element_blank(),
        panel.background=element_blank(),
        legend.position="bottom",
        legend.margin=margin(0, 0, 0, 140),
        legend.background=element_blank(),
        legend.direction="horizontal",
        legend.text=element_text(size=6, margin=margin(0, 5, 0, 3)),
        legend.key.spacing.y=unit(0, "pt"),
        legend.key.size=unit(0.5, "lines"))

gg.sample.biome.its2 <- ggplot(all.samples.bar %>% filter(target_gene == "ITS2"), aes(y=Biome, fill=sample_type)) +
  geom_bar(position="fill", width=0.5) +
  scale_fill_manual(breaks=c("air", "dust", "water", "glacial ice debris", "deadwood", "litter",
                             "coral", "fungal sporocarp", "lichen", "mosses", "shoot",
                             "rhizosphere soil", "soil", "topsoil", "sediment", "root + rhizosphere soil", "root"),
                    values=c('#99DDFF', '#99DDFF', '#77AADD', "dimgrey", '#EEDD88', '#EEDD88',
                             '#FFAABB',  '#DDDDDD', '#DDDDDD', '#BBCC33', '#BBCC33',
                             '#EE8866', '#EE8866', '#EE8866', '#EE8866', '#EE8866', '#EE8866')) +
  labs(fill=NULL, x="Sample types") +
  scale_x_continuous(position="top") +
  scale_y_discrete(limits=rev) +
  guides(fill=guide_legend(nrow=3)) +
  theme(axis.title.x.top=element_text(size=9, vjust=-6),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        plot.background=element_blank(),
        panel.background=element_blank(),
        legend.position="none")

#Plot line plots of sample biomes
gg.sample.values.its1 <- ggplot(all.sample.values %>% filter(target_gene == "ITS1"), aes(x=Biome, y=value)) +
  facet_grid(~measurement, scales="free") +
  geom_violin(linewidth=0.2,
              colour="dimgrey",
              fill="white") +
  geom_boxplot(width=0.1,
               linewidth=0.2,
               outlier.size=0.2,
               outlier.colour="dimgrey",
               colour="dimgrey") +
  geom_point(data=all.sample.sum,
             shape=4,
             stroke=0.8,
             size=0.4,
             colour="black") +
  geom_line(data=all.sample.sum  %>% filter(target_gene == "ITS1"),
            aes(group=1),
            linetype="dashed",
            linewidth=0.2) +
  coord_flip() +
  theme_aplot +
  theme(axis.title.x=element_blank(),
        strip.text=element_text(size=9),
        panel.grid.major.x=element_line(colour="lightgrey"),
        panel.grid.minor.x=element_line(colour="lightgrey"))

gg.sample.values.its2 <- ggplot(all.sample.values %>% filter(target_gene == "ITS2"), aes(x=Biome, y=value)) +
  facet_grid(~measurement, scales="free") +
  geom_violin(linewidth=0.2,
              colour="dimgrey",
              fill="white") +
  geom_boxplot(width=0.1,
               linewidth=0.2,
               outlier.size=0.2,
               outlier.colour="dimgrey",
               colour="dimgrey") +
  geom_point(data=all.sample.sum,
             shape=4,
             stroke=0.8,
             size=0.4,
             colour="black") +
  geom_line(data=all.sample.sum  %>% filter(target_gene == "ITS2"),
            aes(group=1),
            linetype="dashed",
            linewidth=0.2) +
  coord_flip() +
  theme_aplot +
  theme(axis.title.x=element_blank(),
        strip.text=element_text(size=9),
        panel.grid.major.x=element_line(colour="lightgrey"),
        panel.grid.minor.x=element_line(colour="lightgrey"))

#Combine
gg.sample.data.its1 <- gg.sample.biome.its1 %>% insert_right(gg.sample.values.its1, width=2)
gg.sample.data.its2 <- gg.sample.biome.its2 %>% insert_right(gg.sample.values.its2, width=2)

#Summarise proportion of samples per biome
all.samples.sum.its1 <- all.samples.bar %>%
  filter(target_gene == "ITS1") %>%
  mutate(total=nrow(.)) %>%
  group_by(Biome, total, target_gene) %>%
  summarise(n=n()) %>%
  mutate(perc=as.character(round(n/total*100)),
         perc=ifelse(perc == 0, "<1", perc))

all.samples.sum.its2 <- all.samples.bar %>%
  filter(target_gene == "ITS2") %>%
  mutate(total=nrow(.)) %>%
  group_by(Biome, total, target_gene) %>%
  summarise(n=n()) %>%
  mutate(perc=as.character(round(n/total*100)),
         perc=ifelse(perc == 0, "<1", perc))

#Plot biome proportions
gg.biome.prop.its1 <- ggplot(all.samples.sum.its1, aes(x=1, y=n, fill=Biome)) +
  geom_bar(stat="identity",
           position="fill",
           width=0.3) +
  geom_label_repel(aes(label=paste0(Biome, " (", perc, "%)"), x=1.15, colour=Biome),
                   position=position_fillnudge(vjust=0.5, x=0.3),
                   hjust=0,
                   force=3,
                   direction="y",
                   fill="white",
                   size=3,
                   label.size=1,
                   label.padding=unit(3, "pt")) +
  geom_label_repel(aes(label=paste0(Biome, " (", perc, "%)"), x=1.15, group=Biome),
                   position=position_fillnudge(vjust=0.5, x=0.3),
                   hjust=0,
                   force=3,
                   direction="y",
                   fill="white",
                   size=3,
                   label.size=NA,
                   segment.size=NA,
                   label.padding=unit(3, "pt")) +
  scale_fill_manual(breaks=c("anthropogenic", "cropland", "grassland", "shrubland",
                             "woodland", "forest", "mangrove", "wetland", "aquatic", "tundra", "desert"),
                    values=c('#AA4499', '#CC6677', '#999933', '#882255',  
                             '#117733', '#44AA99', "dimgrey",  '#332288', '#88CCEE', '#DDDDDD', '#DDCC77'),
                    na.value=NA) +
  scale_colour_manual(breaks=c("anthropogenic", "cropland", "grassland", "shrubland",
                               "woodland", "forest", "mangrove", "wetland", "aquatic", "tundra", "desert"),
                      values=c('#AA4499', '#CC6677', '#999933', '#882255',  
                               '#117733', '#44AA99', "dimgrey",  '#332288', '#88CCEE', '#DDDDDD', '#DDCC77'),
                      na.value=NA) +
  scale_x_continuous(limits=c(0.8, 3)) +
  scale_y_continuous(expand=c(0, 0)) +
  labs(title="Proportion of GlobalFungi\nsamples in each biome", y="ITS1") +
  theme_minimal() +
  theme(legend.position="none",
        plot.title=element_text(size=10, hjust=0.5),
        plot.margin=margin(5.5, 5.5, 15, 5.5),
        axis.text=element_blank(),
        axis.title.x=element_blank(),
        panel.grid=element_blank())

gg.biome.prop.its2 <- ggplot(all.samples.sum.its2, aes(x=1, y=n, fill=Biome)) +
  geom_bar(stat="identity",
           position="fill",
           width=0.3) +
  geom_label_repel(aes(label=paste0(Biome, " (", perc, "%)"), x=1.15, colour=Biome),
                   position=position_fillnudge(vjust=0.5, x=0.3),
                   hjust=0,
                   force=3,
                   direction="y",
                   fill="white",
                   size=3,
                   label.size=1,
                   label.padding=unit(3, "pt")) +
  geom_label_repel(aes(label=paste0(Biome, " (", perc, "%)"), x=1.15, group=Biome),
                   position=position_fillnudge(vjust=0.5, x=0.3),
                   hjust=0,
                   force=3,
                   direction="y",
                   fill="white",
                   size=3,
                   label.size=NA,
                   segment.size=NA,
                   label.padding=unit(3, "pt")) +
  scale_fill_manual(breaks=c("anthropogenic", "cropland", "grassland", "shrubland",
                             "woodland", "forest", "mangrove", "wetland", "aquatic", "tundra", "desert"),
                    values=c('#AA4499', '#CC6677', '#999933', '#882255',  
                             '#117733', '#44AA99', "dimgrey",  '#332288', '#88CCEE', '#DDDDDD', '#DDCC77'),
                    na.value=NA) +
  scale_colour_manual(breaks=c("anthropogenic", "cropland", "grassland", "shrubland",
                               "woodland", "forest", "mangrove", "wetland", "aquatic", "tundra", "desert"),
                      values=c('#AA4499', '#CC6677', '#999933', '#882255',  
                               '#117733', '#44AA99', "dimgrey",  '#332288', '#88CCEE', '#DDDDDD', '#DDCC77'),
                      na.value=NA) +
  scale_x_continuous(limits=c(0.8, 3)) +
  scale_y_continuous(expand=c(0, 0)) +
  labs(y="ITS2") +
  theme_minimal() +
  theme(legend.position="none",
        plot.margin=margin(15, 5.5, 5.5, 5.5),
        axis.text=element_blank(),
        axis.title.x=element_blank(),
        panel.grid=element_blank())

#Combine
gg.biome.all.its1 <- plot_grid(gg.biome.prop.its1, as.grob(gg.sample.data.its1),
                               rel_widths=c(2, 5), labels="auto")
gg.biome.all.its2 <- plot_grid(gg.biome.prop.its2, as.grob(gg.sample.data.its2),
                               rel_widths=c(2, 5))

gg.biome.all <- plot_grid(gg.biome.all.its1, gg.biome.all.its2, nrow=2, rel_heights=c(1, 0.92))

#Check ITSboth samples that have ITS1 and ITS2 results
samples.filt2 %>%
  filter(permanent_id %in% (all.samples %>%
                              filter(target_gene == "ITSboth") %>%
                              pull(PermanentID))) %>%
  group_by(permanent_id, clade, region) %>%
  summarise(n()) %>%
  group_by(permanent_id) %>%
  filter(n()>1) %>%
  print(n=100)


################################################################################
################################# World maps ###################################
################################################################################

## Base GlobalFungi sampling ##

#Plot GlobalFungi samples on map
gg.map.all.blank <- ggplot(world.map) +
  geom_hex(bins=60,
           data=all.samples,
           aes(x=longitude, y=latitude),
           colour=NA,
           inherit.aes=FALSE) +
  scale_fill_gradient(trans="pseudo_log",
                      low="snow2", high="#DDAA33",
                      breaks=c(1, 10, 100, 2500),
                      limits=c(0, 2700),
                      labels=c(1, 10, 100, 2500)) +
  scale_alpha(range=c(0.2, 1), guide="none") +
  geom_sf(linewidth=0.15,
          fill=NA) +
  labs(title=expression(paste("All ", italic("Gaeumannomyces"))), fill="Number of samples") +
  coord_sf() +
  theme_void() +
  theme(legend.position=c(0.15, 0.4),
        legend.direction="horizontal",
        legend.title.position="top",
        legend.key.size=unit(8, "pt"),
        legend.title=element_text(size=6, face="bold"),
        legend.text=element_text(size=5, margin=margin(2, 0, 0, 0), angle=90, hjust=1, vjust=0.5),
        plot.title=element_text(size=8, hjust=0.5)) +
  ggpreview(width=5, height=3)

#tiff(paste0(dir, "globalfungi_allsamples_blank-", Sys.Date(), ".tiff"), width=4, height=2.3, units="in", res=600, compression="lzw")
gg.map.all.blank
#dev.off()

## ITS1 vs ITS2 ##

#Facet by ITS region
gg.map.all <- ggplot(world.map) +
  facet_wrap(~target_gene, ncol=2) +
  geom_hex(bins=60,
           data=all.samples.comb,
           aes(x=longitude, y=latitude),
           colour=NA,
           inherit.aes=FALSE) +
  scale_fill_gradient(trans="pseudo_log",
                      low="snow2", high="#DDAA33",
                      breaks=c(1, 10, 100, 2500),
                      limits=c(0, 2700),
                      labels=scales::comma(c(1, 10, 100, 2500))) +
  scale_alpha(range=c(0.2, 1), guide="none") +
  geom_sf(linewidth=0.15,
          fill=NA) +
  labs(fill="Number of samples") +
  coord_sf() +
  theme_void() +
  theme(legend.position="bottom",
        legend.key.size=unit(7, "pt"),
        legend.title=element_text(size=6, face="bold", vjust=0.5),
        legend.text=element_text(size=5, margin=margin(2, 0, 0, 0), vjust=0.5),
        strip.text=element_text(size=8, hjust=0.5),
        plot.title=element_text(size=8, hjust=0.5)) +
  ggpreview(width=5, height=2)

gg.allsamples <- plot_grid(gg.biome.all, gg.map.all, ncol=1, rel_heights=c(1, 0.3), labels=c("", "c"))

ggpreview(gg.allsamples, width=7, height=9)

#pdf(paste0(dir, "allsample_comp-", Sys.Date(), ".pdf"), width=7, height=9)
gg.allsamples
#dev.off()

## All Gaeumannomyces on one map ##

#Overlay Gaeumannomyces samples
gg.map.all <- ggplot(world.map) +
  geom_hex(bins=60,
           data=all.samples,
           aes(x=longitude, y=latitude),
           colour=NA,
           inherit.aes=FALSE) +
  scale_fill_gradient(trans="pseudo_log",
                      low="snow2", high="#DDAA33",
                      breaks=c(1, 10, 100, 2500),
                      limits=c(0, 2700),
                      labels=c(1, 10, 100, 2500)) +
 scale_alpha(range=c(0.2, 1), guide=FALSE) +
  geom_sf(linewidth=0.15,
          fill=NA) +
  geom_point(data=samples, 
             aes(y=latitude, x=longitude),
             fill="red",
             colour="black",
             stroke=0.2,
             size=0.6,
             alpha=0.3,
             shape=21) +
  labs(title=expression(paste("All ", italic("Gaeumannomyces"))), fill="Number of samples") +
  coord_sf() +
  theme_void() +
  theme(legend.position=c(0.15, 0.4),
        legend.direction="horizontal",
        legend.title.position="top",
        legend.key.size=unit(8, "pt"),
        legend.title=element_text(size=6, face="bold"),
        legend.text=element_text(size=5, margin=margin(2, 0, 0, 0), angle=90, hjust=1, vjust=0.5),
        plot.title=element_text(size=8, hjust=0.5)) +
  ggpreview(width=5, height=3)

#tiff(paste0(dir, "globalfungi_allsamples-", Sys.Date(), ".tiff"), width=4, height=2.3, units="in", res=600, compression="lzw")
gg.map.all
#dev.off()

## G. tritici ##

#Overlay G. tritici samples
gg.map.gt <- ggplot(world.map) +
  geom_hex(bins=60,
           data=all.samples,
           aes(x=longitude, y=latitude),
           colour=NA,
           inherit.aes=FALSE) +
  scale_fill_gradient(trans="pseudo_log",
                      low="snow2", high="#DDAA33",
                      breaks=c(1, 10, 100, 2500),
                      limits=c(0, 2700),
                      labels=c(1, 10, 100, 2500)) +
  scale_alpha(range=c(0.2, 1), guide=FALSE) +
  geom_sf(linewidth=0.15,
          fill=NA) +
  geom_point(data=seq.samples %>% 
               filter(clade %in% c("G. tritici")), 
             aes(y=latitude, x=longitude),
             fill="red",
             colour="black",
             stroke=0.2,
             size=0.6,
             alpha=0.3,
             shape=21) +
  labs(title=expression(paste(italic("G. tritici"))), fill="Number of samples") +
  coord_sf() +
  theme_void() +
  theme(legend.position=c(0.15, 0.4),
        legend.direction="horizontal",
        legend.title.position="top",
        legend.key.size=unit(8, "pt"),
        legend.title=element_text(size=6, face="bold"),
        legend.text=element_text(size=5, margin=margin(2, 0, 0, 0), angle=90, hjust=1, vjust=0.5),
        plot.title=element_text(size=8, hjust=0.5)) +
  ggpreview(width=5, height=3)

#tiff(paste0(dir, "globalfungi_gt-", Sys.Date(), ".tiff"), width=4, height=2.3, units="in", res=600, compression="lzw")
gg.map.gt
#dev.off()

## G. tritici and G. hyphopodioides ##

#Filter for Gt/Gh
seq.samples.GtGh <- seq.samples %>% 
  filter(clade %in% c("G. tritici", "G. hyphopodioides"))

#Overlay Gt/Gh samples
gg.map.subset <- ggplot(world.map) +
  facet_wrap(~clade, ncol=2) +
  geom_hex(bins=60,
           data=all.samples,
           aes(x=longitude, y=latitude),
           colour=NA,
           inherit.aes=FALSE) +
  scale_fill_gradient(trans="pseudo_log",
                      low="snow2", high="#DDAA33",
                      breaks=c(1, 10, 100, 2500),
                      limits=c(0, 2600),
                      labels=c(1, 10, 100, 2500)) +
  scale_alpha(range=c(0.2, 0.9), guide=FALSE) +
  geom_sf(linewidth=0.15,
          fill=NA) +
  geom_point(data=seq.samples.GtGh, 
             aes(y=latitude, x=longitude),
             fill="red",
             colour="black",
             stroke=0.2,
             size=0.6,
             alpha=0.3,
             shape=21) +
  coord_sf() +
  theme_void() +
  theme(legend.position="none",
        strip.text=element_text(size=8, face="italic", hjust=0.5),
        strip.clip="off") +
  ggpreview(width=5, height=1.5)

#Write maps to file
#pdf(paste0(dir, "globalfungi_maps-", Sys.Date(), ".pdf"), width=4, height=3.5)
plot_grid(gg.map.all, gg.map.subset, nrow=2, rel_heights=c(2, 1))
#dev.off()

#Check South African Gt/Gh samples
gg.map.subset +
  coord_sf(xlim=c(-20, 45), ylim=c(-39, -20), expand=FALSE)

## All Gaeumannomyces on individuals maps ##

#Plot facets for each species individually
gg.map.facet <- ggplot(world.map) +
  facet_wrap(~clade, ncol=5) +
  geom_hex(bins=60,
           data=all.samples,
           aes(x=longitude, y=latitude),
           colour=NA,
           inherit.aes=FALSE) +
  scale_fill_gradient(trans="pseudo_log",
                      low="snow2", high="#DDAA33",
                      breaks=c(1, 10, 100, 2500),
                      limits=c(0, 2700),
                      labels=c(1, 10, 100, 2500)) +
  scale_alpha(range=c(0.2, 0.9), guide=FALSE) +
  geom_sf(linewidth=0.15,
          fill=NA) +
  geom_point(data=seq.samples %>% filter(!is.na(clade)), 
             aes(y=latitude, x=longitude),
             fill="red",
             colour="black",
             stroke=0.1,
             size=0.6,
             alpha=0.3,
             shape=21) +
  labs(fill="Number of samples") +
  theme_void() +
  theme(strip.text=element_text(size=7, face="italic"),
        strip.clip="off",
        legend.position="top",
        legend.direction="horizontal",
        legend.title.position="top",
        legend.key.size=unit(10, "pt"),
        legend.title=element_text(size=6, face="bold"),
        legend.text=element_text(size=5, margin=margin(2, 0, 0, 0), angle=90, hjust=1, vjust=0.5))

#Function to put legend in empty facet cell https://stackoverflow.com/a/54443955
library(gtable)
library(lemon)
shift_legend2 <- function(p) {
  # ...
  # to grob
  gp <- ggplotGrob(p)
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]
  
  # establish name of empty panels
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  names <- empty.facet.panels$name
  
  # now we just need a simple call to reposition the legend
  reposition_legend(p, 'center', panel=names)
}

gg.map.facet.legend <- shift_legend2(gg.map.facet)

#pdf(paste0(dir, "globalfungi_maps_supp-", Sys.Date(), ".pdf"), width=10, height=5)
plot(gg.map.facet.legend)
#dev.off()


################################################################################
################################# Continents ###################################
################################################################################

#Summarise number of sequences per species per continent
seq.samples.continents <- seq.samples %>%
  mutate(continent=sub("Pacific Ocean", "Asia", sub("Atlantic Ocean", "Asia", continent))) %>%
  group_by(clade, continent) %>% 
  summarise(num=n()) %>%
  left_join(all.samples %>% group_by(continent) %>% summarise(samples.num=n()), by="continent") %>%
  mutate(perc=num/samples.num*100) %>%
  bind_rows(gg.summary$data %>% 
              filter(isTip, !label %in% seq.samples$clade) %>% 
              dplyr::select(clade=label) %>% 
              mutate(num=0, continent="Europe")) %>%
  filter(!is.na(clade)) %>%
  mutate(continent=factor(continent,
                          levels=c("Africa", "Antarctica", "Asia", "Australia", "Europe",
                                   "North America", "South America")))

#Plot dotplot of sequences across continents
gg.cont <- ggplot(seq.samples.continents, aes(y=clade, x=continent)) +
  geom_segment(aes(x=continent, xend=continent),
               y=0.5, yend=20.5,
               colour="grey",
               linewidth=0.25,
               linetype="dashed") +
  geom_point(aes(size=perc), shape=21, fill="white") +
  scale_size_continuous(breaks=c(0.07, 0.7, 7),
                        range=c(0.07, 7),
                        name="Percentage of total samples") +
  scale_x_discrete(position="top") +
  coord_cartesian(clip="off") +
  theme_aplot +
  theme(legend.position=c(-0.55, 1.06),
        legend.title.position="top",
        legend.title=element_text(margin=margin(0, 0, 3, 0)),
        legend.box.margin=margin(-15, 0, 0, 0),
        legend.margin=margin(0, 0, 0, 0),
        axis.title.x=element_blank(),
        axis.ticks=element_blank(),
        axis.text.x.top=element_text(angle=-90, hjust=1, vjust=0.5)) +
  ggpreview(height=3.5, width=1.5)

#Combine
gg.cont.tree <- gg.cont %>% insert_left(gg.summary.samples + xlim(0, 30))
ggpreview(gg.cont.tree, width=3, height=3.5)

#pdf(paste0(dir, "continents-", Sys.Date(), ".pdf"), width=3, height=3.5)
gg.cont.tree
#dev.off()

#Combine hosts and continents
gg.host.cont <- gg.host.grid %>%
  insert_left(gg.summary.hosts, width=0.35) %>%
  insert_right(gg.cont + 
                 scale_size_continuous(breaks=c(0.07, 0.7, 7),
                                       range=c(0.07, 7),
                                       name="Percentage of GlobalFungi\nsamples") +
                 theme(legend.position=c(0.5, -0.2)), width=0.3)

ggpreview(gg.host.cont, height=4, width=7.5)

#png(paste0(dir, "hosts-cont-", Sys.Date(), ".png"), width=7.5, height=4, units="in", res=600)
gg.host.cont
#dev.off()

#Check sample types for surprising countries with minimal Gt
all.samples %>%
  group_by(country, sample_type) %>% 
  summarise(num=n()) %>%
  filter(country %in% c("Canada", "Russia", "India")) %>%
  arrange(country, desc(num)) %>%
  print(n=22)

#Check biome per continent
all.samples %>%
  group_by(continent, Biome) %>% 
  summarise(num=n()) %>%
  arrange(continent, desc(num)) %>%
  print(n=67)
