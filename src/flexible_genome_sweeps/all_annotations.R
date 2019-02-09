.libPaths(.libPaths()[2])

library(tidyverse)

combined <- map_df(snakemake@input, ~read_tsv(., col_names=c("protein_accession", "md5", "sequence_length", "analysis", "signature_accession", "signature_description", "start", "stop", "score", "status", "date", "interpro_accession", "interpro_description", "go_annotation", "pathway_annotation")))
write_csv(combined, snakemake@output[["strain_dataframe"]])
