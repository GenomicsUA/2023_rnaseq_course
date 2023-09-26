# https://bioconductor.org/packages/release/bioc/vignettes/recount3/inst/doc/recount3-quickstart.html

library(recount)
library(recount3)

human_projects <- available_projects()

proj_info <- subset(
  human_projects,
  project == "SRP065865" & project_type == "data_sources"
)

rse_gene_MAP <- create_rse(proj_info)

rse_gene_MAP

assay(rse_gene_MAP, "counts") <- transform_counts(rse_gene_MAP)
metadata(rse_gene_MAP)
dim(rse_gene_MAP)
recount3_cols <- colnames(colData(rse_gene_MAP))

assays(rse_gene_MAP)$TPM <- recount::getTPM(rse_gene_MAP)
colSums(assay(rse_gene_MAP, "TPM")) / 1e6

tpms <- assays(rse_gene_MAP)$TPM %>% as.data.frame()
raw_counts <- assays(rse_gene_MAP)$raw_counts %>% as.data.frame()
