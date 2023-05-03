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


##################
#get upstream and downstream interactions of a gene
res=search_gene("GSK3B", interactions)
res$targets_positive$transcriptional
