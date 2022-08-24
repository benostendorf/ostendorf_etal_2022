get_GO_genes <- function(GO_term, w_children = TRUE){
  
  require(GO.db)
  require(org.Mm.eg.db)
  require(biomaRt)
  
  ## Get gene names associated with GO terms & children
  ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  ## Get GO IDs
  go_ids <- GOID( GOTERM[ Term(GOTERM) == GO_term])
  
  ## Get all GO IDs that are children of top term
  if (w_children) {
    go_ids_exp <-  c(GOBPOFFSPRING[[go_ids]], go_ids)
  } else {
    go_ids_exp <- go_ids
  }
  
  ## Get genes associated with GO term
  genes_GO <- 
    getBM(attributes=c('mgi_symbol', 'go_id'),
          filters = 'go', values = go_ids_exp, mart = ensembl) %>%
    pull(mgi_symbol) %>%
    unique()
  
  return(genes_GO)
}