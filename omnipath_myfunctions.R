library(dplyr)

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

# Takes any interaction db and returns those with inhibition or stimulation, 
# consensus direction and removes interactions from "u_source".
conf_interactions_filter <-function(interact_db, u_source){
  interact_conf=unique(interact_db[which((interact_db$is_stimulation
                                          |interact_db$is_inhibition)
                                         &interact_db$consensus_direction
                                         &interact_db$sources!=u_source),])
  return(interact_conf)
}

# From paths list, returns a list with the interactions (pos, neg or undetermined)
find_signed_paths <- function(paths,db){
  test=0
  signed_paths <- paths
  for (j in 1:length(paths)) {
    k=1
    for (i in 1:(length(paths[[j]])-1)) { #for each interaction
      int_db=db[which(db$source_genesymbol==paths[[j]][i]&db$target_genesymbol==paths[[j]][i+1]),]
      
      if(length(int_db$is_stimulation)>1){
        if(length(unique(int_db[1:9]))){
          unique()
        }
      }
      
      if(int_db$is_stimulation&!int_db$is_inhibition){
        signed_interaction= "--(+)-->"
      } else if (!int_db$is_stimulation&int_db$is_inhibition){
        signed_interaction= "--(-)-->"
      } else{
        signed_interaction= "--(u)-->"
      }
      
      signed_paths[[j]][k]<- paths[[j]][i]
      signed_paths[[j]][k+1]<- signed_interaction
      k=k+2
      
      if(i==(length(paths[[j]])-1)){
        signed_paths[[j]][k]<- paths[[j]][i+1] 
      }
    }
  }
  #return(signed_paths)
  return(test)
}

#three levels for regular expressions: 
    #reg=0 -> pattern=^gene$ (default, exact word)
    #reg=1 -> pattern=^gene//character or digit
    #reg=2 -> pattern=gene 
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
  
  genes_db <- list("db" = db, "names" = gene_names)

  
  return(genes_db) 
}
