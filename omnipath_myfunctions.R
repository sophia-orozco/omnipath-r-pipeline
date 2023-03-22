library(dplyr)

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
find_signed_paths <- function(paths,db,type_flag=0,type_option="none"){
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
#targets and sources do not include direct feedback 
#only works with db with following columns: "source_genesymbol"     "target_genesymbol"    
# "is_directed"           "is_stimulation"        "is_inhibition"
search_gene <-function(gene, db,reg=0){
  match=NULL
  gene_names=NULL
  
  if (reg==1){
    gene=c(paste("^",gene, "\\d", sep = ""),paste("^",gene, "\\w", sep = ""))
  }else if (!reg){
    gene=paste("^",gene, "$", sep = "")
  }
    
  for (i in 1:length(gene)) {
    new_match=grep(gene[i], db$target_genesymbol)
    gene_names=c(gene_names,db[new_match,]$target_genesymbol)
    new_match2=grep(gene[i], db$source_genesymbol)
    gene_names=c(gene_names,db[new_match2,]$source_genesymbol)
    match=c(match,new_match,new_match2)
  }
  
  gene_names<-unique(gene_names)
  db= db[unique(match),-c(1,2,5,9,10)]
  
  targets=subset(db, !(target_genesymbol %in% gene_names))
  sources=subset(db, !(source_genesymbol %in% gene_names))
  
  targets_stimulation=unique(filter(targets, is_stimulation&!is_inhibition))
  targets_stimulation=select(targets_stimulation,target_genesymbol, type)
  
  targets_inhibition=unique(filter(targets, !is_stimulation&is_inhibition))
  targets_inhibition=select(targets_inhibition,target_genesymbol, type)
  
  sources_stimulation=unique(filter(sources, is_stimulation&!is_inhibition))
  sources_stimulation=select(sources_stimulation,source_genesymbol, type)
  
  sources_inhibition=unique(filter(sources, !is_stimulation&is_inhibition))
  sources_inhibition=select(sources_inhibition,source_genesymbol, type)
  
  others= subset(db, !(source_genesymbol %in%c(targets_stimulation,targets_inhibition,sources_stimulation,sources_inhibition,gene_names)))
  others= unique(c(others$target_genesymbol,others$source_genesymbol))
  
  dfs <- c("targets_stimulation", "targets_inhibition", "sources_stimulation", "sources_inhibition")
  
  for(df in dfs)
    assign(df, setNames(get(df),  c("name", "type")))
  
  genes_db <- list("db" = db,
                   "names" = gene_names, 
                   "activates"=targets_stimulation,
                   "inhibits"= targets_inhibition,
                   "activated_by" =sources_stimulation,
                   "inhibited_by" = sources_inhibition,
                   "others"=others)
  
  return(genes_db) 
}


if(dim(unique(int_db[,1:10]))[1]==1){
  print(a)
}
