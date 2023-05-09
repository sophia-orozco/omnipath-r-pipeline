setwd("/Users/sophia/Documents/MET_project/omnipath-r-pipeline")
library(OmnipathR)
source("./omnipath_myfunctions.R")

interactions <- import_all_interactions()
interactions_mirna <- import_mirnatarget_interactions()


#filter by references
interactions <- with_references(interactions)

#filter by resource
signor <- filter_by_resource(interactions, resources = 'SIGNOR')

#get upstream and downstream interactions of a gene
p63=search_gene("TP63", interactions)
p63$targets_positive$transcriptional

#define nodes for path
start_nodes <-c("TP63","TCF4","CEBPA")
end_nodes <-c("TP63","TCF4","CEBPA")

start_nodes <-c("ZEB1","ZEB2","hsa-miR-200b-3p")
end_nodes <-c("ZEB1","hsa-miR-200b-3p")

gr_graph <- interaction_graph(interactions)
paths <-find_all_paths(graph = gr_graph, 
                       start = start_nodes, 
                       end = end_nodes, attr = 'name')

#return a list of the signed interactions
find_signed_paths(start_nodes, end_nodes, interactions)

#return a list of the signed interactions and their interaction type
find_signed_paths(start_nodes, end_nodes, interactions, type_flag = 1)

#return a list of the signed interactions filtered by interaction type
find_signed_paths(start_nodes, end_nodes, interactions, type_flag = 1,type_option = "transcriptional")

#plot a graph with interactions 
plot_paths(start_nodes, interactions)

#plot a graph with interactions and their type
plot_paths(start_nodes, interactions, 1)

#plot_neighbours("TFAP2A",interactions,type=0)

search_gene("RUNX1", interactions)

##################
#get upstream and downstream interactions of a gene
res=search_gene("GSK3B", interactions)
res$targets_positive$transcriptional

mirnas=c("hsa-miR-22-3p", "hsa-miR-30a", "hsa-miR-203-3p", "hsa-miR-222-3p")

gene=mirnas[2]
gene="ZEB1"
gene_list=c("ZEB1","GRHL2","TP63")

plot_neighbours(gene,interactions,type=0)
plot_neighbours(gene,interactions,type=1)
plot_neighbours(gene,interactions,type=0, filter = "sources")
plot_neighbours(gene,interactions,type=0, filter = "positive_sources")
plot_neighbours(gene,interactions,type=0, filter = "negative_sources")
plot_neighbours(gene,interactions,type=0, filter = "targets")
plot_neighbours(gene,interactions,type=0, filter = "positive_targets")
plot_neighbours(gene,interactions,type=0, filter = "negative_targets")
#search_gene(gene, interactions_mirna)

search_gene(gene, interactions, "sources")
search_gene(gene, interactions, "positive_sources")
search_gene(gene, interactions, "negative_sources")
search_gene(gene, interactions, "targets")
search_gene(gene, interactions, "positive_targets")
search_gene(gene, interactions, "negative_targets")


set1=c("GRHL2","ZEB1","hsa-miR-200c","TP63")
a=get_gene_interaction(set1,interactions)#[1,]$references
get_gene_interaction(set1,interactions, TRUE)

#build a network from gene list
nodes=read.csv("mymodel.bnet", skip = 2)$targets
a=get_gene_interaction(nodes,interactions)
plot_network(a,type=0)
  

plot_paths(nodes, interactions)

gr_graph <- interaction_graph(interactions)
find_all_paths(graph = gr_graph, 
               start = nodes, 
               end = nodes, attr = 'name')

save_bnet(set1, interactions, include_indirect=FALSE)
  
list_c(search_gene(gene, interactions)$positive_sources)
list_c(search_gene(gene, interactions)$positive_negative)



#test paths: plots and bnet file should have the same output
genes=c("GRHL2","ZEB1")
genes=c("GRHL2","ZEB1","hsa-miR-200c")
genes=c("GRHL2","ZEB1","hsa-miR-200c","TP63","CDH1","TWIST1")


#Only direct interactions
#plot
plot_direct(genes,interactions)
#bnet
save_bnet(genes, interactions,"test1.bnet",include_indirect=FALSE)

#Indluce Indirect interactions
#plot
plot_paths(genes, interactions)
#bnet
save_bnet(genes, interactions,"test2.bnet",include_indirect=TRUE)


