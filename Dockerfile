FROM conda/miniconda3
COPY ./Infomap /Infomap
COPY ./mugsy /mugsy
COPY ./test /test
COPY ./src/PopCOGenT /popcogent
RUN apt-get update && apt-get -y install build-essential perl  && cd Infomap && make
RUN conda install -c bioconda biopython=1.68 joblib=0.9.4 networkx=1.11 snakemake=3.11.2 statsmodels=0.8.0 numpy pandas scipy pyyaml