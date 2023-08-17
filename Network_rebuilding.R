library(MetaboSignal)
library(xlsx)
library(readxl)


#####Load metabolic pathways#####
pae_paths<-data.frame(MS_getPathIds(organism_code = "pae"))
metabo_paths<-pae_paths[pae_paths$Path_type == "metabolic",]$Path_id


#####Metabolic network rebuilding#####
pae_network <- MS_keggNetwork(metabo_paths = metabo_paths, expand_genes = T)
pae_network<-gsub("cpd:","",pae_network)
pae_network<-gsub("pae:","",pae_network)
rows_to_remove <- duplicated(t(apply(pae_network, 1, sort))) & pae_network[,3] == "k_compound:reversible"
filtered_df <- pae_network[!rows_to_remove, ]
write.xlsx(pae_network, "pae_network.xlsx", 
           col.names = TRUE, row.names = F, append = FALSE)

#####Check for missing genes#####
missing_elements <- c()
for ( pathway in metabo_paths) {
  file = paste("https://rest.kegg.jp/link/pae/", pathway, sep = "")
  pathway_url = try(getURL(file), silent = TRUE)
  if (!inherits(pathway_url, "try-error")) {
    cat(paste("Got data for", pathway, "\n"))
    matches <- gregexpr("\\bPA\\d{4}\\b", pathway_url)
    result <- regmatches(pathway_url, matches)
    result_str <- paste(result[[1]], collapse = " ")
    result_split <- unlist(strsplit(result_str, " "))
    result_with_prefix <- paste("pae:", result_split, sep = "")
    cat(paste(pathway,"has",length(result_with_prefix),"genes", "\n"))
    cat(paste("The network miss",length(result_with_prefix[!result_with_prefix %in% network_vec]),"genes", "\n"))
    
    #check if the elements are in the network. if not, collect elements that are not in the network in a seperate list
    missing_elements <- c(missing_elements,result_with_prefix[!result_with_prefix %in% network_vec])
  } else {
    cat(paste("Failed to get data for pathway", pathway, "\n"))
  }
  
  
} 

cat(paste("The network miss total", length(unique(missing_elements)), "genes", "\n"))
all_pathway_missing<-unique(missing_elements)
length(all_pathway_missing) #98


missing_pathways<-c("pae01100","pae00190","pae00543")
missing_elements <- c()
for (pathway in missing_pathways) {
  file = paste("https://rest.kegg.jp/link/pae/", pathway, sep = "")
  pathway_url = try(getURL(file), silent = TRUE)
  if (!inherits(pathway_url, "try-error")) {
    cat(paste("Got data for", pathway, "\n"))
    matches <- gregexpr("\\bPA\\d{4}\\b", pathway_url)
    result <- regmatches(pathway_url, matches)
    result_str <- paste(result[[1]], collapse = " ")
    result_split <- unlist(strsplit(result_str, " "))
    result_with_prefix <- paste("pae:", result_split, sep = "")
    #result_with_prefix <- gsub("\"", "", result_with_prefix)
    #length(result_with_prefix)
    cat(paste(pathway,"has",length(result_with_prefix),"genes", "\n"))
    cat(paste("The network miss",length(result_with_prefix[!result_with_prefix %in% network_vec]),"genes", "\n"))
    
    #check if the elements are in the network. if not, collect elements that are not in the network in a seperate list
    missing_elements <- c(missing_elements,result_with_prefix[!result_with_prefix %in% network_vec])
  } else {
    cat(paste("Failed to get data for pathway", pathway, "\n"))
  }
  
  
} 

cat(paste("The network miss total", length(unique(missing_elements)), "genes", "\n"))
missing_pathway_missing<-unique(missing_elements)
length(missing_pathway_missing) #137


#Elements only in all_pathway_missing
length(all_pathway_missing[!all_pathway_missing %in% missing_pathway_missing])#42
#add these genes manually in the network to complete the network