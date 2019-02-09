.libPaths(.libPaths()[2])
library(tidyverse)

gene_df <- read_csv(snakemake@input[["mmseqs"]])
cluster_df <- read_tsv(snakemake@input[["subclusters"]])
colnames(cluster_df) <- c("fasta_id", "subcluster", "main_cluster_id", "subcluster_id", "clonal_complex")

cluster_df$subcluster <- as.double(cluster_df$subcluster)

orf_df <- read_csv(snakemake@input[["orf"]])
annotation_df <- read_csv(snakemake@input[["annotations"]])

cluster_df <- semi_join(cluster_df, gene_df, by="fasta_id")


for(sub_cluster_id in unique(cluster_df$subcluster)) {
  my_fasta_ids <- (cluster_df %>% filter(subcluster==sub_cluster_id))$fasta_id
  gene_df %>% group_by(cluster) %>% filter(setequal(fasta_id, my_fasta_ids))

  out_df <- gene_df %>% 
    group_by(cluster) %>%
    filter(setequal(fasta_id, my_fasta_ids)) %>%
    left_join(orf_df, suffix=c(".mmseqs", ".prodigal"), by=c("member"="protein_accession")) %>%
    left_join(annotation_df, by=c("member"="protein_accession"), suffix=c(".mmseqs",".interproscan")) %>%
    select(-X1, -contig, -fasta_id.prodigal, -orf.prodigal) %>%
    select(fasta_id.mmseqs, member, cluster, analysis, signature_description, interpro_description, everything()) %>%
    arrange(cluster, fasta_id.mmseqs)

  if (sub_cluster_id == 0){
    sub_cluster_id = '0.0'
  }
  out_df %>%
    write_csv(paste0(snakemake@params[["all_orfs"]], ".", sub_cluster_id))

  out_df %>%
    select(cluster, analysis, signature_accession, signature_description, interpro_accession, interpro_description) %>% ungroup() %>% distinct() %>%
    write_csv(paste0(snakemake@params[["orfs_summary"]], ".", sub_cluster_id))
}