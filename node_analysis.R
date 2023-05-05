setwd("/Users/sophia/Documents/MET_project/omnipath-r-pipeline")
library(OmnipathR)
source("./omnipath_myfunctions.R")

interactions <- import_all_interactions()

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
plot_paths(start_nodes, end_nodes, interactions)

#plot a graph with interactions and their type
plot_paths(start_nodes, end_nodes, interactions, 1)

#plot_neighbours("TFAP2A",interactions,type=0)

search_gene("RUNX1", interactions)

##################
#get upstream and downstream interactions of a gene
res=search_gene("GSK3B", interactions)
res$targets_positive$transcriptional

mirnas=c("hsa-miR-22-3p", "hsa-miR-30a", "hsa-miR-203-3p", "hsa-miR-222-3p")

gene=mirnas[2]
gene="ZEB1"

plot_neighbours(gene,interactions,type=0)
plot_neighbours(gene,interactions,type=1)
plot_neighbours(gene,interactions,type=0, filter = "sources")
plot_neighbours(gene,interactions,type=0, filter = "positive_sources")
plot_neighbours(gene,interactions,type=0, filter = "negative_sources")
plot_neighbours(gene,interactions,type=0, filter = "targets")
plot_neighbours(gene,interactions,type=0, filter = "negative_targets")
plot_neighbours(gene,interactions,type=0, filter = "positive_targets")

search_gene(gene, interactions)$sources_negative
