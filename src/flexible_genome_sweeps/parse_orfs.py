from Bio import SeqIO
import pandas as pd
import re


def parse_fasta(fasta):
    with open(fasta, "r") as f:
        for rec in SeqIO.parse(f, "fasta"):
            yield parse_record(rec)


def parse_record(record):
    strain_contig_orf, orf_start, orf_end, orf_strand, prod_info = record.description.split(" # ")
    strain, contig, orf = re.match(r"(.*)_([^_]+)_(\d+)$", strain_contig_orf).groups()

    (prod_id, prod_partial_left, prod_partial_right,
        prod_start_type, prod_rbs_motif, prod_rbs_spacer,
        prod_gc_cont) = parse_prod_info(prod_info)

    return {
        "fasta_id": strain,
        "contig": contig,
        "orf": orf,
        "protein_accession": strain_contig_orf,
        "orf_start": orf_start,
        "orf_end": orf_end,
        "orf_strand": orf_strand,
        "prod_id": prod_id,
        "prod_partial_left": prod_partial_left,
        "prod_partial_right": prod_partial_right,
        "prod_start_type": prod_start_type,
        "prod_rbs_motif": prod_rbs_motif,
        "prod_rbs_spacer": prod_rbs_spacer,
        "prod_gc_cont": prod_gc_cont,
        "orf_seq": str(record.seq)
    }


def parse_prod_info(prod_info):
    (prod_id, prod_partial,
        prod_start_type, prod_rbs_motif, prod_rbs_spacer,
        prod_gc_cont) = [field.split("=")[1] for field in prod_info.split(";")]
    prod_partial_left, prod_partial_right = prod_partial
    return (prod_id, prod_partial_left, prod_partial_right,
            prod_start_type, prod_rbs_motif, prod_rbs_spacer,
            prod_gc_cont)

df = pd.DataFrame(parse_fasta(snakemake.input.fasta))
df[['orf_start', 'orf_end', 'prod_gc_cont']] = df[['orf_start', 'orf_end', 'prod_gc_cont']].apply(pd.to_numeric)
df.to_csv(snakemake.output.csv)
