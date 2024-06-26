init <- function(){
    library(biomaRt)
    library(tidyverse)
}

init_mart_human <- function(reference = "hg38"){
     if (reference == "grch37"){
         host = "https://grch37.ensembl.org"
     }else{
         host = "https://useast.ensembl.org"
     }
  
    mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = host)
    mart <- useDataset(mart, dataset = "hsapiens_gene_ensembl")
    return(mart)
}


build_master_table <- function(){
     targets_csv <- "PANELTargets.genes.csv"
     targets <- read_csv(targets_csv) %>% dplyr::select(ensembl_gene_id, entrezgene_ids, 
                                                        gene_name_hg19, gene_name_hg38)
     targets$panel <- "PANELTargets"
     
     mito_csv <- "FullMITO.genes.csv"
     mito <- read_csv(mito_csv) %>% dplyr::select(ensembl_gene_id, entrezgene_ids, 
                                                        gene_name_hg19, gene_name_hg38)
     mito$panel <- "FullMITO"
     
     pid_csv <- "PID.genes.csv"
     pid <- read_csv(pid_csv) %>% dplyr::select(ensembl_gene_id, entrezgene_ids, 
                                                  gene_name_hg19, gene_name_hg38)
     pid$panel <- "PID"
     
     master_panel <- bind_rows(targets, mito, pid)
     
     write_excel_csv(master_panel, "master_panel.csv")
     
     genes_sanger <- master_panel %>% group_by(ensembl_gene_id, entrezgene_ids, gene_name_hg19,
                                               gene_name_hg38) %>%
                       summarise(panels = paste(panel, collapse = ";"))
     
     sanger_no <- read_csv("no_sanger.csv")
     genes_sanger <- genes_sanger %>% left_join(sanger_no, 
                                                by = c("gene_name_hg19" = "gene_name")) %>%
                  replace_na(list(sanger_backfill = "SANGER_YES"))
     
     write_excel_csv(genes_sanger, "genes_sanger.csv")
}

# add entrez ids
get_entregene_id <- function(){
    genes_csv <- "TARGETS.genes.csv"
    genes <- read_csv(genes_csv)
    entrez_ann <- as_tibble(getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
                                filters = c("ensembl_gene_id"),
                                values = list(genes$ensembl_gene_id),
                                mart = mart)) 
    entrez_ann$entrezgene_id <- as.numeric(entrez_ann$entrezgene_id)
    entrez_ann <- entrez_ann %>% group_by(ensembl_gene_id) %>% 
      arrange(entrezgene_id) %>% 
        summarise(entrezgene_ids = paste(entrezgene_id, collapse = ":"))
    
    genes <- genes %>% left_join(entrez_ann, by = c("ensembl_gene_id" = "ensembl_gene_id"))
    # genes$ensembl_gene_id38 <- genes$ensembl_gene_id
    write_excel_csv(genes, genes_csv)
}

# coordinates of the exon starts and ends
# some exons of the canonical isoform are non coding - they have NA in genomic_coding_start,
# they are excluded, for example in ACTA1
# -001 is NOT always the canonical isoform
# takes the longest cds from gencode_basic transcripts
# PLEC gene has two ensembl identifiers:
# ENSG00000178209 - we need this one, which is protein_coding
# ENSG00000261109
# this is a slow function, use it for individual genes/gene panels, when every exon is needed
# for bulk coverage analysis use the function above
get_all_exon_coordinates <- function(ensembl_gene_id, mart, genes){
    # ensembl_gene_id <- "ENSG00000073734"
    # ensembl_gene_id <- "ENSG00000117054" 
    # ensembl_gene_id <- "ENSG00000184056"
    # ensembl_gene_id <- "ENSG00000171759" #  two entrez ids
    # ensembl_gene_id <- "ENSG00000211899" # no refseq transcript no entrez
    # ensembl_gene_id <- "ENSG00000229164" 
    ensembl_gene_id <- "ENSG00000270141"
  
    # print(ensembl_gene_id)
    genes_info <- as_tibble(getBM(attributes = c("chromosome_name",
                                               "genomic_coding_start", "genomic_coding_end",
                                               "ensembl_exon_id", "external_gene_name", "ensembl_gene_id",
                                               "start_position", "end_position",
                                               "exon_chrom_start", "exon_chrom_end",
                                               "ensembl_transcript_id"),
                                filters = c("ensembl_gene_id"),
                                values = c(ensembl_gene_id), 
                                mart = mart)) %>%
                dplyr::filter(!grepl("PATCH", chromosome_name)) %>%
                dplyr::filter(!grepl("HSCHR", chromosome_name)) %>% 
                dplyr::mutate(chromosome_name = as.character(chromosome_name)) %>%
                dplyr::filter(!is.na(genomic_coding_start) & !is.na(genomic_coding_end)) %>%
                dplyr::select(chromosome_name, genomic_coding_start, genomic_coding_end, 
                  external_gene_name, ensembl_gene_id, ensembl_exon_id, ensembl_transcript_id)
    
    # this is non-coding gene
    if (nrow(genes_info) == 0){
        genes_info <- as_tibble(getBM(attributes = c("chromosome_name",
                                                   "ensembl_exon_id", "external_gene_name", "ensembl_gene_id",
                                                   "start_position", "end_position",
                                                   "exon_chrom_start", "exon_chrom_end",
                                                   "ensembl_transcript_id"),
                                    filters = c("ensembl_gene_id"),
                                    values = c(ensembl_gene_id), 
                                    mart = mart)) %>%
        dplyr::filter(!grepl("PATCH", chromosome_name)) %>%
        dplyr::filter(!grepl("HSCHR", chromosome_name)) %>% 
        dplyr::mutate(chromosome_name = as.character(chromosome_name)) %>%
        dplyr::mutate(genomic_coding_start = exon_chrom_start) %>%
        dplyr::mutate(genomic_coding_end = exon_chrom_end) %>%
        dplyr::filter(!is.na(genomic_coding_start) & !is.na(genomic_coding_end)) %>%
        dplyr::select(chromosome_name, genomic_coding_start, genomic_coding_end, 
                      external_gene_name, ensembl_gene_id, ensembl_exon_id, ensembl_transcript_id)
    }
    
    query_transcripts <- genes_info$ensembl_transcript_id %>% unique()
    
    # some ENS  transcripts have refseq_mrna, some don't
    # there could be >1 refseq ids for ENSEMBL_TRANSCRIPT, we collapse them with :
    refseq_mrna <- as_tibble(getBM(attributes = c("ensembl_transcript_id", "refseq_mrna"),
                         filters = c("ensembl_transcript_id"),
                         values = list(query_transcripts),
                         mart = mart))  %>%
                 dplyr::filter(!is.na(refseq_mrna)) %>%
                 dplyr::filter(refseq_mrna != "") %>%
                 group_by(ensembl_transcript_id) %>% 
                 summarise(refseq_ids = paste(refseq_mrna, collapse = ":")) %>%
                 dplyr::rename(refseq_mrna = refseq_ids)
    
    
    # NC genes don't have refseq_mrna they have refseq_ncrna
    if (nrow(refseq_mrna) == 0 ){
         refseq_mrna <- as_tibble(getBM(attributes = c("ensembl_transcript_id", "refseq_ncrna"),
                                     filters = c("ensembl_transcript_id"),
                                     values = list(query_transcripts),
                                     mart = mart))  %>%
              dplyr::filter(!is.na(refseq_ncrna)) %>%
              dplyr::filter(refseq_ncrna != "") %>%
              dplyr::mutate(refseq_mrna = refseq_ncrna) %>%
              group_by(ensembl_transcript_id) %>% 
              summarise(refseq_ids = paste(refseq_mrna, collapse = ":")) %>%
              dplyr::rename(refseq_mrna = refseq_ids)
    }
    
    # genes without REFseq in grch37 but with it in grch38
    if (nrow(refseq_mrna) == 0 && genes[genes$ensembl_gene_id == ensembl_gene_id,]$refseq_status == "grch38"){
        refseq_mrna[1,"ensembl_transcript_id"] <- genes[genes$ensembl_gene_id == ensembl_gene_id,]$ensembl_transcript_hardcoded
        refseq_mrna[1,"refseq_mrna"] <- genes[genes$ensembl_gene_id == ensembl_gene_id,]$refseq_hardcoded
    }
    
    # genes without REFSeq ids
    if (nrow(refseq_mrna) == 0 && genes[genes$ensembl_gene_id == ensembl_gene_id,]$refseq_status == "NOREFSEQ"){
        refseq_mrna <- genes_info$ensembl_transcript_id %>% unique() %>% as_tibble() %>% dplyr::rename(ensembl_transcript_id = value)
        refseq_mrna$refseq_mrna <- rep("NOREFSEQ", nrow(refseq_mrna))
    }
    
        
        # MAP3K14 transcript does not have mapping to refseq in grch37
        # we hardcode it
        # in grch38 it is a protein coding gene with mapping:
        # https://useast.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000006062;r=17:45263119-45317029
        # https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000006062;r=17:43340488-43394414
        #if (nrow(refseq_mrna) == 0  && ensembl_gene_id == "ENSG00000006062"){
        #      refseq_mrna[1,"ensembl_transcript_id"] <- "ENST00000344686"
        #      refseq_mrna[1,"refseq_mrna"] <- "NM_003954"
        #      
        #      genes_info <- genes_info %>% dplyr::filter(ensembl_transcript_id == "ENST00000344686")
        #}
         
        # IGHM does not have refseq id 
        #if (nrow(refseq_mrna) == 0 && ensembl_gene_id == "ENSG00000211899") {
        #  refseq_mrna[1,"ensembl_transcript_id"] <- "ENST00000390559"
        #  refseq_mrna[1,"refseq_mrna"] <- "NA"
          
        #  genes_info <- genes_info %>% dplyr::filter(ensembl_transcript_id == "ENST00000390559")
        #}
         
        # TRAC in PID does not have refseq_id
        # if (nrow(refseq_mrna) == 0 && ensembl_gene_id == "ENSG00000229164") {
        #   refseq_mrna[1,"ensembl_transcript_id"] <- "ENST00000478163"
        #   refseq_mrna[1,"refseq_mrna"] <- "NA"
           
        #   genes_info <- genes_info %>% dplyr::filter(ensembl_transcript_id == "ENST00000478163")
        # } 
        # }
    
    # some turn boolean
    refseq_mrna$refseq_mrna <- as.character(refseq_mrna$refseq_mrna)
    
    # also annotate with entrezgene_id for compatibility with Emedgene
    # ENSG00000171759 has two entrez ids, we pick up the small one
    entrez_ann <- as_tibble(getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
                                   filters = c("ensembl_gene_id"),
                                   values = ensembl_gene_id,
                                   mart = mart))
    
    entrez_ann$entrezgene_id <- as.numeric(entrez_ann$entrezgene_id)
    entrez_ann <- entrez_ann %>% arrange(entrezgene_id) %>% head(1)
    # some genes don't have entrez id, ENST00000390559 = IGHM
    if (nrow(entrez_ann) == 0 ){
        entrez_ann[1, "ensembl_gene_id"] <- ensembl_gene_id
        entrez_ann[1, "entrezgene_id"] <- "NA"
    }
    
    # Don't annotate with gencode_basic - returns only transcripts which are gencode basic
    # all our transcripts with Refseq ID and without, are gencode basic
    # query_transcripts <- refseq_mrna$ensembl_transcript_id %>% unique()
    # gencode_basic <- as_tibble(getBM(attributes = c("ensembl_transcript_id", "transcript_gencode_basic"),
    #                              filters = c("ensembl_transcript_id"),
    #                              values = list(query_transcripts),
    #                              mart = mart)) 
    #gencode_basic$transcript_gencode_basic <- as.character(gencode_basic$transcript_gencode_basic)
    #gencode_basic$transcript_gencode_basic <- ifelse(gencode_basic$transcript_gencode_basic == "1", "GENCODE_BASIC", "GENCODE_NA")
    #
    #refseq_mrna <- refseq_mrna %>% left_join(gencode_basic, 
    #                                         by = c("ensembl_transcript_id" = "ensembl_transcript_id"))
    
    # remove exons without refseq_mrna
    genes_info <- genes_info %>% left_join(refseq_mrna, 
                                           by = c("ensembl_transcript_id" = "ensembl_transcript_id")) %>%
                  dplyr::filter(refseq_mrna != "") %>%
                  left_join(entrez_ann, by = c("ensembl_gene_id" = "ensembl_gene_id")) %>%
                  dplyr::select(chromosome_name, genomic_coding_start, genomic_coding_end,
                                external_gene_name, refseq_mrna, 
                                entrezgene_id, ensembl_gene_id, ensembl_transcript_id, 
                                ensembl_exon_id)
    
    return(genes_info)
}

# default padding is 15 bp
get_exon_coordinates <- function(genes_csv, exons_bed, mart, padding = 15){
    genes <- read_csv(genes_csv)
    exons <- get_all_exon_coordinates(genes$ensembl_gene_id[1], mart, genes)
  
    for(gene in tail(genes$ensembl_gene_id, -1)){
        exons_buf <- get_all_exon_coordinates(gene, mart, genes)
        exons <- bind_rows(exons, exons_buf)
    }
  
   padding <- as.numeric(padding)
   # need one more base at 5'
   left_padding <- padding + 1
   exons$genomic_coding_start <- exons$genomic_coding_start - left_padding
   exons$genomic_coding_end <- exons$genomic_coding_end + padding
   
   write_tsv(exons, exons_bed, col_names = F)
}

# start and end of the gene, all exons
# input = vector of genes, either ENSEMBL_IDS or external names = disease_panel.list.txt, no header
# output = bed file with coordinates = disease_panel.list.bed
# output is not sorted please sort with bedtools or bash sort
get_gene_coordinates <- function(v_ensembl_gene_ids, panel_genes_bed, mart){
    genes <- getBM(
        attributes = c("ensembl_gene_id", "chromosome_name",
                       "start_position", "end_position",
                       "external_gene_name"),
        filters = c("ensembl_gene_id"),
        values = v_ensembl_gene_ids,
        mart = mart)

    genes_bed <- as_tibble(genes)

    genes_bed <- genes_bed %>%
                    dplyr::select(chromosome_name, start_position, end_position, external_gene_name, ensembl_gene_id) %>%
                    dplyr::filter(str_detect(chromosome_name, "PATCH", negate = T)) %>%
                    dplyr::filter(str_detect(chromosome_name, "HSCHR", negate = T))

    write_tsv(genes_bed, panel_genes_bed, col_names = F)
}

args <- commandArgs(trailingOnly = T)
if (length(args) != 5 || args[1] == "--help"){
    cat("Usage: Rscript grch37|hg38 genes.csv panel.raw.bed padding_bp panel.genes.raw.bed\n")
}else{
    reference <- args[1]
    genes_csv <- args[2]
    panel_bed <- args[3]
    padding <- args[4]
    panel_genes_bed <- args[5]
    init()
    mart <- init_mart_human(reference)
    get_exon_coordinates(genes_csv, panel_bed, mart, padding)
    genes <- read_csv(genes_csv)
    get_gene_coordinates(genes$ensembl_gene_id, panel_genes_bed, mart)
}

