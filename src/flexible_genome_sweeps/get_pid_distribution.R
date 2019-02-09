.libPaths(.libPaths()[2])
library(Biostrings)
library(multidplyr)
library(tidyverse)


gene_df <- read_csv(snakemake@input[["master_presence_absence"]])
orf_df <- read_csv(snakemake@input[["orfs"]])
cluster_df <- read_tsv(snakemake@input[["subclusters"]], col_names=c("fasta_id", "clonal_group", "subcluster"))
cluster_df$subcluster <- as.double(cluster_df$subcluster)
cluster_df <- semi_join(cluster_df, gene_df, by="fasta_id")

avg_pairwise_identity <- function(orf_seqs) {
    pairs <- combn(orf_seqs, 2)
    scores <- c()
    for (i in seq_len(ncol(pairs))) {
        scores <- append(scores, pid(pairwiseAlignment(DNAString(pairs[,i][1]), DNAString(pairs[,i][2])), type="PID2"))
    }
    mean(scores)
}


for(sub_cluster_id in unique(cluster_df$subcluster)) {
    special_df_fn <- paste0("output/", snakemake@wildcards[["organism"]], "/", snakemake@wildcards[["organism"]], ".unique_orfs_summary.csv.", sub_cluster_id)
    special_df <- read_csv(special_df_fn)
    my_fasta_ids <- (cluster_df %>% filter(subcluster==sub_cluster_id))$fasta_id

    party_df1 <- gene_df %>% group_by(cluster) %>% filter(all(my_fasta_ids %in% fasta_id)) %>% ungroup() %>% filter(fasta_id %in% my_fasta_ids) %>% group_by(cluster) %>% filter(n()==length(my_fasta_ids)) %>% ungroup() %>% left_join(orf_df, suffix=c(".mmseqs", ".prodigal"), by=c("member"="protein_accession"))
    if (nrow(party_df1) == 0){
      next
    }
    party_df1 <- party_df1 %>% partition(cluster)
    cluster_eval(get_default_cluster(), source("util.R"))
    party_df2 <- party_df1 %>% select(fasta_id.prodigal, cluster, member, orf_seq) %>% group_by(cluster) %>% summarise(avg_pid=avg_pairwise_identity(orf_seq))
    party_df3 <- collect(party_df2)

    special_df_fn <- paste0("output/", snakemake@wildcards[["organism"]], "/", snakemake@wildcards[["organism"]], ".unique_orfs_summary.csv.", sub_cluster_id)
    special_df <- read_csv(special_df_fn)

    g <- party_df3 %>% 
      mutate(unique_to_cluster=ifelse(cluster %in% special_df$cluster, "unique_to_cluster", "not_unique")) %>% 
      ggplot(aes(avg_pid)) +
      geom_histogram(aes(fill=unique_to_cluster), binwidth=1, boundary=0.5, position="dodge") +
      scale_fill_brewer(palette="Set1", direction=-1) +
      theme_linedraw()

    ggsave(plot=g, filename=paste0(snakemake@params[["plot"]], ".", sub_cluster_id, ".png"))

    g <- party_df3 %>% 
      mutate(unique_to_cluster=ifelse(cluster %in% special_df$cluster, "unique_to_cluster", "not_unique")) %>% 
      ggplot(aes(avg_pid)) +
      geom_density(aes(color=unique_to_cluster)) +
      scale_color_brewer(palette="Set1", direction=-1) +
      theme_linedraw()

    ggsave(plot=g, filename=paste0(snakemake@params[["plot_density"]], ".", sub_cluster_id, ".png"))

    write_csv(party_df3, paste0(snakemake@params[["csv"]], ".", sub_cluster_id))
}