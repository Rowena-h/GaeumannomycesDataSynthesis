library(ape)

tree <- read.tree("W:/GaeumannomycesMetaanalysis/raxmlng/gaeumannomyces_concat.raxml.support")

tree.rooted <- root(tree, c("Magnaporthe_rhizophila_isolate_M23", "Magnaporthiopsis_maydis_strain_CBS_662.82A",
                            "Magnaporthiopsis_sp._CPC_26038", "Magnaporthe_poae_isolate_M48",
                            "Magnaporthiopsis_maydis_strain_CBS_133165", "Gaeumannomyces_incrustans_isolate_M35"), edgelabel=TRUE, resolve.root=TRUE)

write.tree(tree.rooted, "W:/GaeumannomycesMetaanalysis/raxmlng/gaeumannomyces_concat.raxml.support_rooted")