# Example: `snakemake -s PopCOGenT.Snakefile --use-conda -rp -c $SLURM_NTASKS`


import glob
import pandas as pd
import itertools
import os


MODULE_PATH = config.get("PopCOGenT-path", "./")
infomap_param = config.get("infomap_param", "")
single_cell = config.get("single_cell", False)


if config.get("test", False):
    shell("realpath ../../test/*.fasta > test_M16x_genomes.fasta.ls")

    clonal_cutoff = config.get("clonal_cutoff", 0.000355362)

    rule target:
        input:
            graphml=f"test_M16x_genomes-length_bias-{clonal_cutoff}.unclust.graphml",
            cluster=f"test_M16x_genomes-length_bias-{clonal_cutoff}.cluster.tsv",


rule cluster:
    input:
        edge='{prefix}-length_bias-{clonal_cutoff}.edge.tsv',
    output:
        cluster='{prefix}-length_bias-{clonal_cutoff}.cluster.tsv',
    params:
        module=MODULE_PATH,
        tmpdir="{prefix}-length_bias_proc/infomap",
        infomap_args=infomap_param,
    conda:
        "PopCOGenT.yml"
    shell:
        """
        mkdir -p {params.tmpdir}

        PYTHONPATH="{params.module}" python <<!EOF!
        \nif True:
        import cluster

        cluster.make_clusterfile(
            "{input.edge}",
            "{output.cluster}",
            "{params.tmpdir}",
            "{params.infomap_args}",
        )

        \n!EOF!
        """


rule cluster_make_edgefile:
    input:
        all_bias="{prefix}-length_bias.tsv",
    output:
        edge='{prefix}-length_bias-{clonal_cutoff}.edge.tsv',
    params:
        module=MODULE_PATH,
        clonal_cutoff="{clonal_cutoff}",
        single_cell=single_cell,
    conda:
        "PopCOGenT.yml"
    shell:
        """
        PYTHONPATH="{params.module}" python <<!EOF!
        \nif True:
        import cluster

        cluster.make_edgefile(
            "{input.all_bias}",
            "{output.edge}",
            clonal_cutoff={params.clonal_cutoff},
            single_cell={params.single_cell},
            linear_model=cluster.negative_selection_linear_fit()
        )

        \n!EOF!
        """


rule cluster_make_graphml:
    input:
        edge='{prefix}-length_bias-{clonal_cutoff}.edge.tsv',
    output:
        graphml="{prefix}-length_bias-{clonal_cutoff}.unclust.graphml",
    params:
        clonal_cutoff="{clonal_cutoff}",
        single_cell=False,
    conda:
        "PopCOGenT.yml"
    shell:
        """
        python <<!EOF!
        \nif True:
        import networkx as nx

        nx.write_graphml(
            nx.read_edgelist("{input.edge}", data=(('weight', float),)),
            "{output.graphml}",
        )

        \n!EOF!
        """


checkpoint ls_genome_files:
    input:
        ls="{any}.fasta.ls"
    output:
        ls="{any}-length_bias_proc/genomes.tsv"
    run:
        with open(input.ls) as fi:
            genomes = {
                os.path.splitext(os.path.basename(v))[0]: v
                for v
                in (i.strip() for i in fi)
            }
        with open(output.ls, "w") as fo:
            for g, p in genomes.items():
                if os.path.splitext(p)[1] in (".fa", ".fna", ".fasta"):
                    fo.write(f"{g}\t{p}\n")


def expand_genome_bias(bias="{any}-length_bias_proc/alignment/{g1}_@_{g2}.length_bias.txt"):
    def _f(wildcards):
        mags_file = checkpoints.ls_genome_files.get(**wildcards).output.ls
        with open(mags_file) as fi:
            genomes = {
                k: v
                for k, v
                in (i.strip().split() for i in fi)
            }
        return [
            bias.format(g1=g1, g2=g2, **wildcards)
            for g1, g2
            in itertools.combinations(genomes, 2)
        ]
    return _f


rule concat_length_bias_files:
    input:
        expand_genome_bias(
            bias="{any}-length_bias_proc/alignment/{g1}_@_{g2}.length_bias.txt"
        ),
    output:
        bias="{any}-length_bias.tsv",
    params:
        module=MODULE_PATH
    conda:
        "PopCOGenT.yml"
    shell:
        """
        PYTHONPATH="{params.module}" python <<!EOF!
        \nif True:
        import length_bias_functions

        length_bias_functions.concat_length_bias_files(
            "{output.bias}", "{input}".split()
        )
        \n!EOF!
        """


rule make_length_bias_file:
    input:
        g1 = "{any}-length_bias_proc/renamed_mugsy/{g1}.fa",
        g2 = "{any}-length_bias_proc/renamed_mugsy/{g2}.fa",
        maf = "{any}-length_bias_proc/mugsy/{g1}_@_{g2}.maf",
    output:
        bias = "{any}-length_bias_proc/alignment/{g1}_@_{g2}.length_bias.txt",
    params:
        module=MODULE_PATH
    conda:
        "PopCOGenT.yml"
    shell:
        """
        PYTHONPATH="{params.module}" python <<!EOF!
        \nif True:
        import length_bias_functions

        length_bias_functions.calculate_length_bias(
            "{input.maf}", "{input.g1}", "{input.g2}", "{output.bias}"
        )
        \n!EOF!
        """


def musgy_extract_genome(f_basename="{g1}"):
    def _f(wildcards):
        mags_file = checkpoints.ls_genome_files.get(**wildcards).output.ls
        with open(mags_file) as fi:
            for line in fi:
                genome_name, genome_path = line.strip().split()
                if genome_name == f_basename.format(**wildcards):
                    return genome_path
    return _f


rule musgy_format_input:
    input:
        g1 = musgy_extract_genome(f_basename="{g1}"),
    output:
        g1 = "{any}-length_bias_proc/renamed_mugsy/{g1}.fa",
    params:
        g1="{g1}",
        module=MODULE_PATH
    wildcard_constraints:
        suffix="fa|fna|fasta"
    conda:
        "PopCOGenT.yml"
    shell:
        """
        PYTHONPATH="{params.module}" python <<!EOF!
        \nif True:
        import length_bias_functions

        length_bias_functions.rename_for_mugsy(
            "{input.g1}", "{input.g2}", "{params.g1}"
        )
        \n!EOF!
        """


rule musgy_align:
    input:
        g1 = "{any}-length_bias_proc/renamed_mugsy/{g1}.fa",
        g2 = "{any}-length_bias_proc/renamed_mugsy/{g2}.fa",
    output:
        maf = "{any}-length_bias_proc/mugsy/{g1}_@_{g2}.maf",
    log:
        maf = "{any}-length_bias_proc/mugsy/{g1}_@_{g2}.mugsy.log",
    conda:
        "PopCOGenT.yml"
    shadow:
        "shallow"
    shell:
        """
        /bin/rm -rf smk-musgy
        mkdir smk-musgy

        cp {input.g1} smk-musgy/1.tempcontig
        cp {input.g2} smk-musgy/2.tempcontig

        (
            cd smk-musgy/
            mugsy \
                --directory `pwd` \
                --prefix out \
                1.tempcontig \
                2.tempcontig
        )

        mv smk-musgy/out.maf {output.maf}
        mv smk-musgy/out.mugsy.log {log.maf}
        """
