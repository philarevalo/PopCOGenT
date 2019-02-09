import glob
from Bio import SeqIO
import pandas as pd
import argparse


def main():
    parser = argparse.ArgumentParser(description=('Takes genbank files and creates a fasta file with contigs renamed with strain name.'),
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--table_name',
                        default=None,
                        type=str,
                        help='Mapping table filename')
    parser.add_argument('--genome_ext',
                        default='gbff',
                        type=str,
                        help='Extension for genome files (e.g., .gbff')
    parser.add_argument('--genome_dir',
                        default=None,
                        type=str,
                        help='Directory containing genome files.')
    parser.add_argument('--output_dir',
                        default='./',
                        type=str,
                        help='Directory for new genomes.')

    args = parser.parse_args()

    outfile_name = args.table_name
    genome_ext = args.genome_ext
    genome_dir = args.genome_dir
    output_dir = args.output_dir

    genome_files = glob.glob('{input_dir}/*.{genome_ext}'.format(input_dir=genome_dir,
                                                                 genome_ext=genome_ext))
    output = pd.DataFrame(index=genome_files, columns=['gbff_file', 'fasta_file', 'strain_name', 'iso_source'])

    if outfile_name and genome_dir and genome_ext:
        for f in genome_files:
            # Append genome file name
            seqs = [s for s in SeqIO.parse(f, 'genbank')]

            # Deletes forbidden character entirely
            try:
                strain = seqs[0].features[0].qualifiers['strain'][0]
                strain = strain.translate(({ord(c): '' for c in """ !@#$%^&*()[]{};:,./<>?\|`"'~-=+"""}))
            except KeyError:
                strain = seqs[0].id

            try:
                iso_source = seqs[0].features[0].qualifiers['isolation_source'][0]
            except KeyError:
                iso_source = 'Unknown'

            accession_0 = seqs[0].id
            # Not the best way to do this, will break if filename is contained
            # in directory path for some reason
            fasta_outfile = '{output_dir}/{strain}--{acc}.fasta'.format(strain=strain,
                                                                        output_dir=output_dir,
                                                                        acc=accession_0)

            # Makes new fasta file
            new_seqs = []
            for i, s in enumerate(seqs):
                s.id = '{strain}_contig_{num}_{acc}'.format(strain=strain,
                                                            num=str(i),
                                                            acc=s.id)
                s.description = ''
                new_seqs.append(s)

            SeqIO.write(new_seqs, fasta_outfile, 'fasta')

            output.loc[f, 'gbff_file'] = f
            output.loc[f, 'fasta_file'] = fasta_outfile
            output.loc[f, 'strain_name'] = strain
            output.loc[f, 'iso_source'] = iso_source

        output.to_csv(outfile_name, sep='\t', index=None)


if __name__ == '__main__':
    main()