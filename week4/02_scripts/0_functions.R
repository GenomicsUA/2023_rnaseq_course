dataframe_to_enrichr <- function(dataframe){
  tmp = dataframe %>% tidyr::separate("Overlap", sep = "/", into = c("Gene_number", "max"))
  tmp$Gene_number = as.numeric(tmp$Gene_number)
  tmp = tmp %>% arrange(desc(Gene_number))
}
