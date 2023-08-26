#!/usr/bin/Rscript

# adapted from https://github.com/dpryan79/Answers/tree/master/SEQanswers_42420

library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)

GTFfile = "../data/Homo_sapiens.GRCh38.110.gtf"
FASTAfile = "../data/Homo_sapiens.GRCh38.dna.toplevel.fa"

# Load the annotation and reduce it
# asRangedData = F
GTF <- import.gff(GTFfile, format="gtf", genome="GRCh38", feature.type="exon")
grl <- reduce(split(GTF, elementMetadata(GTF)$gene_id))
reducedGTF <- unlist(grl, use.names=T)
elementMetadata(reducedGTF)$gene_id <- rep(names(grl), elementNROWS(grl))

#Open the fasta file
FASTA <- FaFile(FASTAfile)
open(FASTA)

#Add the GC numbers
elementMetadata(reducedGTF)$nGCs <- letterFrequency(getSeq(FASTA, reducedGTF), "GC")[,1]
elementMetadata(reducedGTF)$widths <- width(reducedGTF)

#Create a list of the ensembl_id/GC/length
calc_GC_length <- function(x) {
    nGCs = sum(elementMetadata(x)$nGCs)
    width = sum(elementMetadata(x)$widths)
    c(width, nGCs/width)
}
output <- t(sapply(split(reducedGTF, elementMetadata(reducedGTF)$gene_id), calc_GC_length))
colnames(output) <- c("Length", "GC")

write.table(output, file = "../tables/GC_lengths.tsv", sep="\t")
