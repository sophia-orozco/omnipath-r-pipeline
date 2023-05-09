library(dplyr)
library(purrr)
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
  #...
  
  #remove genes with no stimuation and no inhibition
  db_final=db[-which(db$is_stimulation==0&db$is_inhibition==0),]
  
  #check for mirna in db
  mirna=grep("hsa",db$source_genesymbol)
  #if no mirna in database keep only stimulation or inhibition
  if(length(mirna)){
    db_mirna=db[mirna,]
    db_final=bind_rows(db_final,db_mirna)
  }
  
  #filter option only for plot_neighbours
  if(filter!="none"& gene!="none"){
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
  #ecol[E(gr_graph)$is_inhibition==1]="brown1"
  #if it has both positive and negative interactions
  ecol[E(gr_graph)$is_inhibition==1&E(gr_graph)$is_stimulation==1]="blue"
  
  etype=(rep("solid", ecount(gr_graph)))
  if(type){
    etype[(E(gr_graph)$type=="transcriptional"|
             E(gr_graph)$type=="mirna_transcriptional")]="solid"
    etype[(E(gr_graph)$type=="transcriptional"|
               E(gr_graph)$type=="mirna_transcriptional")]="solid"
    etype[(E(gr_graph)$type=="post_transcriptional"|
            E(gr_graph)$type=="post_translational"|
              E(gr_graph)$type=="lncrna_post_transcriptional")]="dashed"
  } 
  #plot(gr_graph, edge.color=ecol)
  g=plot(gr_graph, edge.color=ecol,edge.lty=etype, vertex.label.color="black",
       vertex.label.cex=.8, vertex.size=20, vertex.shape="circle",vertex.color="white",
       vertex.label.family="Helvetica",edge.arrow.size=0.5,layout=layout_with_fr) #edge.curved=0.8
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
# plot direct interactions between a set of genes
plot_direct <- function(genes,db,type=0){
  db=get_gene_interaction(genes,db)
  res=plot_network(db,type)
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
search_gene <-function(gene, db, filter="none", reg=0){
  match=NULL
  gene_names=NULL
  #list_temp=NULL
  
  #regular expressions
  if (reg==1){
    gene=c(paste("^",gene, "\\d", sep = ""),paste("^",gene, "\\w", sep = ""))
    for (i in 1:length(gene)) {
      #identify targets and sources
      upstream=grep(gene[i], db$target_genesymbol)
      gene_names=c(gene_names, db[upstream,]$target_genesymbol)
      upstream_names=db[upstream,]$source_genesymbol
      downstream=grep(gene[i], db$source_genesymbol)
      gene_names=c(gene_names, db[downstream,]$source_genesymbol)
      downstream_names=db[downstream,]$target_genesymbol
      match=c(match,upstream,downstream)
    }
  }else{
    #identify targets and sources
    upstream=which(db$target_genesymbol==gene)
    gene_names=c(gene_names, db[upstream,]$target_genesymbol)
    upstream_names=db[upstream,]$source_genesymbol
    downstream=which(db$source_genesymbol==gene)
    gene_names=c(gene_names, db[downstream,]$source_genesymbol)
    downstream_names=db[downstream,]$target_genesymbol
    match=c(match,upstream,downstream)
  }

  gene_names<-unique(gene_names)
  db= db[match,]#-c(1,2,5,9,10)
  
  #check for mirna in db
  mirna=grep("hsa",db$source_genesymbol)
  if(length(mirna)){
    db_mirna=db[mirna,]
  }
  
  targets=db %>% filter(target_genesymbol %in% downstream_names)
  sources=db %>% filter(source_genesymbol %in% upstream_names)
  names(targets)[names(targets)=="target_genesymbol"] <- "name"
  names(sources)[names(sources)=="source_genesymbol"] <- "name"
  
  targets_positive=list()
  targets_negative=list()
  sources_positive=list()
  sources_negative=list()
  
  if(length(downstream_names)){
    #positive
    transcriptional=unique(filter(targets, (type=="transcriptional"|type=="mirna_transcriptional") & is_stimulation==1))$name
    post_translational=unique(filter(targets, type=="post_translational" & is_stimulation==1)$name)
    
    targets_positive<- list("transcriptional"=transcriptional,"post_translational"=post_translational)
    targets_positive=targets_positive[lapply(targets_positive, length)>0]
    
    #negative
    transcriptional=unique(filter(targets, (type=="transcriptional"|type=="mirna_transcriptional") & is_inhibition==1))$name
    post_translational=unique(filter(targets, type=="post_translational" & is_inhibition==1)$name)
    
    if(length(mirna)){
     post_transcriptional=unique(filter(targets, type=="post_transcriptional" & is_directed==1))$name
    }else{
      post_transcriptional=list()
    }
    targets_negative<- list("transcriptional"=transcriptional,"post_translational"=post_translational,"post_transcriptional"=post_transcriptional)
    targets_negative=targets_negative[lapply(targets_negative, length)>0]
  }
  if(length(upstream_names)){
    #positive
    transcriptional=unique(filter(sources, (type=="transcriptional"|type=="mirna_transcriptional") & is_stimulation==1))$name
    post_translational=unique(filter(sources, type=="post_translational" & is_stimulation==1))$name
    
    sources_positive<- list("transcriptional"=transcriptional,"post_translational"=post_translational)
    sources_positive=sources_positive[lapply(sources_positive, length)>0]
    
    #negative
    transcriptional=unique(filter(sources, (type=="transcriptional"|type=="mirna_transcriptional") & is_inhibition==1))$name
    post_translational=unique(filter(sources, type=="post_translational" & is_inhibition==1))$name
    
    if(length(mirna)){
      post_transcriptional=unique(filter(sources, type=="post_transcriptional" & is_directed==1))$name
    }else{
      post_transcriptional=list()
    }
    sources_negative<- list("transcriptional"=transcriptional,"post_translational"=post_translational,"post_transcriptional"=post_transcriptional)
    sources_negative=sources_negative[lapply(sources_negative, length)>0]
  }

  
  gene_db<-list(gene_names,db, sources_positive, sources_negative, targets_positive,targets_negative)
  names(gene_db) <- c("gene_name","db","positive_sources","negative_sources","positive_targets","negative_targets")
  gene_db=gene_db[lapply(gene_db, length)>0]
  
  if(filter!="none"){
    if(filter=="targets"){
      gene_db<- list(gene_db$positive_targets,gene_db$negative_targets)
      names(gene_db) <- c("positive_targets","negative_targets")
    }else if(filter=="sources"){
      gene_db<- list(gene_db$positive_sources,gene_db$negative_sources)
      names(gene_db) <- c("positive_sources","negative_sources")
    }else{
      gene_db=gene_db[filter]
    }
  }
  
  return(gene_db) 
}

#another function with more info about two genes
get_gene_interaction<-function(gene1,db, short_db=FALSE){
  genes<-NULL
  gene2<-gene1
  #genes<-list(c(gene1,gene2),c(gene2,gene1))
  k=1
  for (i in 1:length(gene1)) {
    for (j in 1:length(gene2)) {
      genes_temp1=c(gene1[i],gene2[j])
      genes_temp2=c(gene2[j],gene1[i])
      genes [[k]]<-(genes_temp1)
      genes [[k+1]]<-(genes_temp2)
      k=k+2
    }
  }
  genes=unique(genes)
  gene_int=find_paths_db(genes,db)
  if(short_db){
    gene_int=gene_int[,c(3,4,5,6,7,14)] 
  }
  return(gene_int)
}

#generate model. using find all paths -> find paths_db -> search gene for each node
#option for just direct interactions vs find all paths

#if gene is both activator and repressor

save_bnet <-function(genes, db, filename, include_indirect=FALSE){
  if(include_indirect){
    #get full db
    gr_graph <- interaction_graph(db)
    paths <-find_all_paths(graph = gr_graph, 
                           start = start_nodes, 
                           end = end_nodes, attr = 'name')
    db=find_paths_db(paths,db)
    genes=unique(c(paths_db$target_genesymbol,paths_db$source_genesymbol))
  }else{
    db=get_gene_interaction(genes,db)
    
    #remove nodes without sources from genes list and db 
    genes=genes[genes %in% db$target_genesymbol]
    db=subset(db, source_genesymbol %in% genes)
    db=distinct(db, source_genesymbol, target_genesymbol, is_directed,
                is_stimulation,is_inhibition, .keep_all = TRUE)
  }
  
    positive=rep(list(c(0)),length(genes))
    negative=rep(list(c(0)),length(genes))
    factors=rep(0,length(genes))
    #search genes in targets and get its sources
    for (i in 1:length(genes)) {
      activators=0
      repressors=0
      pos=search_gene(genes[i], db)$positive_sources
      neg=search_gene(genes[i], db)$negative_sources
      
      #add positive and negative relations in a list
      if(!is.null(pos)){
        positive[[i]]=list_c(pos)
        positive[[i]]=unique(positive[[i]])
        if(length(positive[[i]])==1){
          activators=paste("(",positive[[i]][1],")", sep = "")
        }else{
          activators=paste("(",positive[[i]][1], sep = "")
          for (j in 2:length(positive[[i]])) {
            activators=paste(activators,"|",positive[[i]][j], sep = "")
          }
          activators=paste(activators,")", sep = "")
        }
      }
      if(!is.null(neg)){
        negative[[i]]=list_c(neg)
        negative[[i]]=unique(negative[[i]])
        if(length(negative[[i]])==1){
          repressors=paste("(",negative[[i]][1],")", sep = "")
        }else{
          repressors=paste("(",negative[[i]][1], sep = "")
          for (j in 2:length(negative[[i]])) {
            repressors=paste(repressors,"|",negative[[i]][j], sep = "")
          }
          repressors=paste(repressors,")", sep = "")
        }
      }
      
      if(activators!=0&repressors!=0){
        #if the same gene is both in activators and repressors, remove it
        
        factors[i]=paste(activators," & !",repressors, sep = "")
      }else if(activators!=0&repressors==0){
        factors[i]=activators
      }else{
        factors[i]=paste("!",repressors, sep = "")
      }
      
      }
    
  #save file
  bnet=data.frame(targets=genes, factors)
  return(bnet)
  #write(bnet, filename)
  }


