library(dplyr)
library(OmnipathR)
library(igraph)
library(igraphdata)

##################################################################################
# find paths nodes interactions in any db.
find_paths_db <- function(paths,db){
  interac_db = NULL
  for (j in 1:length(paths)) {
    for (i in 1:(length(paths[[j]])-1)) {
      interac_db <- bind_rows(interac_db, db[which(db$source_genesymbol==paths[[j]][i]
                                                   &db$target_genesymbol==paths[[j]][i+1]),]) 
    }
  }
  return(interac_db)
}
##################################################################################
# get interactions between genes. Return interaction sign, type and reference
#get_gene_interaction <- function(gene1, gene2,db){
  
#    upstream=grep(gene1, db$target_genesymbol)
#    gene_names=c(gene_names, db[upstream,]$target_genesymbol)
#    upstream_names=db[upstream,]$source_genesymbol
#    downstream=grep(gene[i], db$source_genesymbol)
#    gene_names=c(gene_names, db[downstream,]$source_genesymbol)
#    downstream_names=db[downstream,]$target_genesymbol
#    match=c(match,upstream,downstream)
  
#  gene_names<-unique(gene_names)
#  db= db[match,]
#  return(interac_db)
#}
##################################################################################
#plot network
plot_network <- function(db,type=0, filter="none",gene="none"){
  
  #check database type (mirna vs all)
  
  db_final=db[-which(db$is_stimulation==0&db$is_inhibition==0),]
  
  #check for mirna in db
  mirna=grep("hsa",db$source_genesymbol)
  #if no mirna in database keep only stimulation or inhibition
  if(length(mirna)){
    db_mirna=db[mirna,]
    db_final=bind_rows(db_final,db_mirna)
  }
  
  #filter option only for plot_neighbours
  if(filter!="none"){
    if(filter=="sources"){
      db_final=filter(db_final, (target_genesymbol==gene))
    }else if(filter=="targets"){
      db_final=filter(db_final, (source_genesymbol==gene))
    }else if(filter=="positive_targets"){
      db_final=filter(db_final, (source_genesymbol==gene& is_stimulation==1))
    }else if(filter=="negative_targets"){
      #if gene is mirna
      if(length(grep("hsa", gene))){ 
        db_final=filter(db_final, (source_genesymbol==gene& is_directed==1))
      }else{
        db_final=filter(db_final, (source_genesymbol==gene& is_inhibition==1))
      }
    }else if(filter=="positive_sources"){
      db_final=filter(db_final, (target_genesymbol==gene& is_stimulation==1))
    }else if(filter=="negative_sources"){
      db_final=filter(db_final, (target_genesymbol==gene& (is_inhibition==1)))
      #if source is an mirna
      if(length(mirna)){
        db_mirna=filter(db_mirna, (target_genesymbol==gene& (is_directed==1)))
        db_final=bind_rows(db_final,db_mirna)
      }
    }
  }
  
  gr_graph <- interaction_graph(db_final)
  
  ecol=(rep("gray", ecount(gr_graph)))
  ecol[E(gr_graph)$is_directed==1]="red1"
  ecol[E(gr_graph)$is_stimulation==1]="chartreuse4"
  ecol[E(gr_graph)$is_inhibition==1]="brown1"
  #if it has both positive and negative interactions
  ecol[E(gr_graph)$is_inhibition==1&E(gr_graph)$is_stimulation==1]="blue"
  
  etype=(rep("solid", ecount(gr_graph)))
  if(type){
    etype[(E(gr_graph)$type=="transcriptional"|
               E(gr_graph)$type=="mirna_transcriptional")]="solid"
    etype[(E(gr_graph)$type=="post_transcriptional"|
            E(gr_graph)$type=="post_translational"|
              E(gr_graph)$type=="lncrna_post_transcriptional")]="dashed"
  } 
  
  plot(gr_graph, edge.color=ecol, vertex.label.color="black",
       vertex.label.cex=.8, vertex.size=20, vertex.shape="circle",vertex.color="white",
       vertex.label.family="Helvetica",edge.curved=0.8,
       edge.arrow.size=0.5,layout=layout_with_fr,edge.lty=etype)
  if(type){
    legend("bottomleft",legend = c("transcriptional","post_transcriptional/translational"),
           lty=c("solid","dashed"), cex = 0.8, bty = "n", ncol = 1)
  }
  
  plot_res <- recordPlot()
  
  return(plot_res)
}

##################################################################################
# plot gene neighbours
plot_neighbours <- function(gene,db,type=0, filter="none"){
  db=search_gene(gene, db, reg=0)$db
  res=plot_network(db, type, filter, gene)
  return(res)
}

##################################################################################
# plot graph from interactions
plot_paths <- function(start_nodes,end_nodes,db,type=0){
  gr_graph <- interaction_graph(db)
  paths <-find_all_paths(graph = gr_graph, 
                         start = start_nodes, 
                         end = end_nodes, attr = 'name')
  paths_db=find_paths_db(paths,db)
  res=plot_network(paths_db, type)
  return(res)
}
##################################################################################
# Takes any interaction db and returns those with inhibition or stimulation, 
# consensus direction and removes interactions from "u_source".
conf_interactions_filter <-function(interact_db, u_source){
  interact_conf=unique(interact_db[which((interact_db$is_stimulation
                                          |interact_db$is_inhibition)
                                         &interact_db$consensus_direction
                                         &interact_db$sources!=u_source),])
  return(interact_conf)
}

##################################################################################
# From paths list, returns a list with the interactions (pos, neg or undetermined)
# Type=1 returns interaction type
# make sure the db used to find the paths is the same one used here 
find_signed_paths <- function(start_nodes,end_nodes,db,type_flag=0,type_option="none"){
  gr_graph <- interaction_graph(db)
  paths <-find_all_paths(graph = gr_graph, 
                 start = start_nodes, 
                 end = end_nodes, attr = 'name')
  signed_paths <- paths
  for (j in 1:length(paths)) {
    k=1
    for (i in 1:(length(paths[[j]])-1)) { 
      # select db for each interaction
      int_db=db[which(db$source_genesymbol==paths[[j]][i]&db$target_genesymbol==paths[[j]][i+1]),]
      
      # if no interaction is found (when list is only one dimensional or not in db)
      if(!dim(int_db)[1]){
        break
      }
      
      # if more types of interaction are found
      if(length(int_db$is_stimulation)>1){
      # if the interactions do not have the same sign make them have the same sign
          if(!dim(unique(int_db[,1:10]))[1]==1){
              if(1 %in% int_db$is_stimulation){ #if at least one is a stimulation
                int_db$is_stimulation==rep(1, dim(int_db)[1])
              }
              if(1 %in% int_db$is_inhibition){ #if at least one is an inhibition
                int_db$is_inhibition==rep(1, dim(int_db)[1])
              }
          }
        # keep one row and collapse interaction type
        int_db$type=paste(int_db$type,collapse='/')
        int_db$type=unique(int_db$type)
        int_db=int_db %>% distinct(type, .keep_all = TRUE)
      }
      
      if(int_db$is_stimulation&!int_db$is_inhibition){
        signed_interaction= "--(+)-->"
      } else if (!int_db$is_stimulation&int_db$is_inhibition){
        signed_interaction= "--(-)-->"
      } else if (int_db$is_stimulation&int_db$is_inhibition){
        signed_interaction= "--(b)-->"
      } else{
        signed_interaction= "--(u)-->"
      }
      if(type_flag){
        type=int_db$type
        signed_interaction=paste(signed_interaction, type, sep = " ")
      }
      signed_paths[[j]][k]<- paths[[j]][i]
      signed_paths[[j]][k+1]<- signed_interaction
      k=k+2
      
      if(i==(length(paths[[j]])-1)){
        signed_paths[[j]][k]<- paths[[j]][i+1] 
      }
    }
    if(type_flag&type_option!="none"){ #evaluate for every path 
      if(!grepl(type_option, signed_paths[j])){ #if not interaction type, empty (with lenght 1)
        signed_paths[j]<-0
      }
    }
  }
  #returns only paths longer than 1 element 
  signed_paths=signed_paths[lengths(signed_paths)>1]
  return(signed_paths)
}

##################################################################################
#three levels for regular expressions: 
    #reg=0 -> pattern=^gene$ (default, exact word)
    #reg=1 -> pattern=^gene//character or digit
    #reg=2 -> pattern=gene 
#db must have the following columns: "source_genesymbol"     "target_genesymbol"    
# "is_directed"           "is_stimulation"        "is_inhibition"
search_gene <-function(gene, db, reg=0){
  match=NULL
  gene_names=NULL
  
  #regular expressions
  if (reg==1){
    gene=c(paste("^",gene, "\\d", sep = ""),paste("^",gene, "\\w", sep = ""))
  }else if (!reg){
    gene=paste("^",gene, "$", sep = "")
  }
  
  #identify targets and sources
  for (i in 1:length(gene)) {
    upstream=grep(gene[i], db$target_genesymbol)
    gene_names=c(gene_names, db[upstream,]$target_genesymbol)
    upstream_names=db[upstream,]$source_genesymbol
    downstream=grep(gene[i], db$source_genesymbol)
    gene_names=c(gene_names, db[downstream,]$source_genesymbol)
    downstream_names=db[downstream,]$target_genesymbol
    match=c(match,upstream,downstream)
  }
  
  gene_names<-unique(gene_names)
  db= db[match,]#-c(1,2,5,9,10)
  targets=db %>% filter(target_genesymbol %in% downstream_names)
  sources=db %>% filter(source_genesymbol %in% upstream_names)
  names(targets)[names(targets)=="target_genesymbol"] <- "name"
  names(sources)[names(sources)=="source_genesymbol"] <- "name"
  
  targets_positive<- list(
    "transcriptional"=unique(filter(targets, (type=="transcriptional"|type=="mirna_transcriptional") & is_stimulation==1))$name,
    "post_translational"=unique(filter(targets, type=="post_translational" & is_stimulation==1)$name))
  targets_negative<- list(
    "transcriptional"=unique(filter(targets, (type=="transcriptional"|type=="mirna_transcriptional") & is_inhibition==1))$name,
    "post_translational"=unique(filter(targets, type=="post_translational" & is_inhibition==1)$name))
  sources_positive<- list(
    "transcriptional"=unique(filter(sources, (type=="transcriptional"|type=="mirna_transcriptional") & is_stimulation==1)$name,
    "post_translational"=unique(filter(sources, type=="post_translational" & is_stimulation==1))$name))
  sources_negative<- list(
    "transcriptional"=unique(filter(sources, (type=="transcriptional"|type=="mirna_transcriptional") & is_inhibition==1))$name,
    "post_translational"=unique(filter(sources, type=="post_translational" & is_inhibition==1))$name)
  
  gene_db <-list(gene_names,db, targets_positive,targets_negative, sources_positive, sources_negative)
  names(gene_db) <- c("gene_name","db","targets_positive", "targets_negative","sources_positive", "sources_negative")
  
  return(gene_db) 
}

#another function with more info about two genes


