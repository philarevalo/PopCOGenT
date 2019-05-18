import argparse
import networkx as nx
import pandas as pd
import itertools
import numpy as np
from collections import defaultdict
import os
import glob
import pandas as pd
import statsmodels.api as sm
import numpy as np


def main():
    parser = argparse.ArgumentParser(
        description=('Clusters genomes based on length bias measurement'),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--base_name',
                        type=str,
                        help='Base name for outfiles')
    parser.add_argument('--length_bias_file',
                        type=str,
                        default=None,
                        help='Table of length bias measurements.')
    parser.add_argument('--clonal_cutoff',
                        type=float,
                        default=0.000355362)
    parser.add_argument('--output_directory',
                        type=str,
                        default='./infomap_temp/',
                        help='Output directory.')
    parser.add_argument('--infomap_args',
                        type=str,
                        default='',
                        help='Arguments to pass to infomap command')
    parser.add_argument('--infomap_path',
                        type=str,
                        default=None,
                        help='Path to infomap binary')
    parser.add_argument('--single_cell', default=False, action='store_true')
    
    # Defining variables from input
    args = parser.parse_args()

    # Checks inputs to see if they're valid
    check_inputs(args)
    infile = args.length_bias_file
    infomap_path = args.infomap_path
    infomap_args = args.infomap_args
    outfile_base = args.base_name
    clonal_cutoff = args.clonal_cutoff
    output_dir = args.output_directory
    single_cell = args.single_cell

    # Initializing lists and dictionaries
    final_clusters = defaultdict(list)
    initial_edgefile = 'infomap_out/{base}_{mindiv}.txt'.format(base=outfile_base,
                                                                 mindiv=str(clonal_cutoff))
    cluster_file = '{output_dir}/{initial_edgefile}.cluster.tab.txt'.format(initial_edgefile=os.path.basename(initial_edgefile),
                                                                            output_dir=output_dir)
    graphml_unclust_name = '{output_dir}/{initial_edgefile}.unclust.graphml'.format(initial_edgefile=os.path.basename(initial_edgefile),
                                                                                    output_dir=output_dir)

    # Make the initial network edgefile, graph object, and graphml file
    make_edgefile(infile,
                  initial_edgefile,
                  clonal_cutoff=clonal_cutoff,
                  single_cell=single_cell,
                  linear_model=negative_selection_linear_fit())
    G_unclust = nx.read_edgelist(initial_edgefile, data=(('weight', float),))
    nx.write_graphml(G_unclust, graphml_unclust_name)

    # Loops over each connected component of the initial network
    for i, component in enumerate(nx.connected_component_subgraphs(G_unclust)):

        # First checks if the component contains more than 1 node
        if(len(component.nodes())) > 1:
            outname = '{initial_edgefile}_{clustnum}.net'.format(initial_edgefile=initial_edgefile,
                                                                 clustnum=str(i))
            infomap_outfile = '{initial_edgefile}_{clustnum}.tree'.format(initial_edgefile=initial_edgefile,
                                                                          clustnum=str(i))

            nx.write_pajek(component, outname)

            # Modifes pajek file to a format the infomap understands
            new_lines = []
            at_edges = False
            with open(outname, 'r') as old_pajek:
                for line in old_pajek:
                    if 'edges' in line:
                        at_edges = True
                    if '*' not in line and not at_edges:
                        splits = line.split()
                        splits[1] = '"' + splits[1] + '"'
                        new_lines.append(' '.join(splits + ['\n']))
                    elif 'vertices' in line:
                        new_lines.append(line.replace('vertices', 'Vertices'))
                    elif 'edges' in line:
                        new_lines.append(line.replace('edges', 'Edges'))
                    else:
                        new_lines.append(line)
            with open(outname, 'w') as outfile:
                outfile.writelines(new_lines)

            # Runs infomap
            cmd = '{infomap} -i pajek {infomap_args} {pajek_file} infomap_out/'.format(infomap=infomap_path,
                                                                                       pajek_file=outname,
                                                                                       infomap_args=infomap_args)
            os.system(cmd)

            # Parse the final cluster file

            # Get subclusters and number them
            subcluster_strings = []
            for line in open(infomap_outfile, 'r'):
                if line[0] != '#' and line != '':
                    cluster_string, flow, name, node = line.split()
                    subcluster_strings.append(':'.join(cluster_string.split(':')[0:-1]))
            subcluster_strings = sorted(list(set(subcluster_strings)))

            # Fill up the final cluster dictionary
            for line in open(infomap_outfile, 'r'):
                if line[0] != '#' and line != '':
                    cluster_string, flow, name, node = line.split()
                    subcluster_string = ':'.join(cluster_string.split(':')[0:-1])
                    cluster_id = '{cluster}.{subcluster}'.format(cluster=str(i),
                                                                 subcluster=str(subcluster_strings.index(subcluster_string)))
                    final_clusters[cluster_id].append(name.replace('"', ''))

        else:
            final_clusters[str(i) + '.0'].append(component.nodes()[0])

    # Enforces that final populations must be cliques if there is an obvious reason to cut out a node (i.e., a single max clique)
    cluster_dict = {}
    for clust, strains in final_clusters.items():
        for clonal_complex in strains:
            cluster_dict[clonal_complex] = clust
    remove = []
    for u, v in G_unclust.edges():
        if cluster_dict[u] != cluster_dict[v]:
            remove.append((u, v))
    for u, v in remove:
        G_unclust.remove_edge(u, v)
    clusters = G_unclust
    max_cluster = int(max([float(clust) for clust in final_clusters.keys()]))
    for c in nx.connected_components(clusters):
        G = clusters.subgraph(c)
        cliques = [(len(cli), cli) for cli in nx.find_cliques(G)]
        if len(cliques) != 1:
            max_clique_size = max(cliques)[0]
            exclude = []
            for size, clique in cliques:
                if len(clique) == max_clique_size:
                    exclude.append(set(G.nodes()) - set(clique))
            if len(exclude) == 1:
                for node in exclude[0]:
                    clust = cluster_dict[node]
                    final_clusters[clust].remove(node)
                    final_clusters[str(float(max_cluster + 1))] = [node]
                    max_cluster += 1
                    print(node)

    with open(cluster_file, 'w') as final:
        final.write('\t'.join(['Strain',
                               'Cluster_ID',
                               'Main_cluster',
                               'Sub_cluster',
                               'Clonal_complex\n']))
        for clust, strains in final_clusters.items():
            for clonal_complex in strains:
                for clone in clonal_complex.split(','):
                    final.write('\t'.join([clone,
                                           clust,
                                           clust.split('.')[0],
                                           clust.split('.')[1],
                                           clonal_complex + '\n']))
def check_inputs(args):

    # Check for output directory. Makes it if it isn't there.
    if not os.path.exists(args.output_directory):
        print('Ouput directory does not exist. Creating new directory.')
        os.makedirs(args.output_directory)

    if not os.path.exists('infomap_out/'):
        print('Ouput directory does not exist. Creating new directory.')
        os.makedirs('infomap_out')

    # Checks infomap
    if not os.path.exists(args.infomap_path):
        raise FileNotFoundError('Invalid infomap path.')


def negative_selection_linear_fit():
    # Hard-coded negative selection cutoff
    m = 1.8361818809179542
    b = 0.96292873154210135
    se = 0.3852360007448502

    X = np.array((641905.0, 1894187.75, 2329260.25, 4673513.5))
    Y = np.array((1.78799146525, 4.85288564669, 7.218796311329999, 9.45229717906))

    x = sm.add_constant(np.array(X)/1e6)
    model = sm.OLS(Y, x)
    results = model.fit()
    return results


def make_edgefile(infile,
                  outfile_name,
                  clonal_cutoff=0,
                  single_cell=False,
                  linear_model=None):

    final_edges = []

    trn_table = pd.read_table(infile)  # Table containing transfer measurement results from "calculate_ssd.py"
    all_strains = set(list(trn_table['Strain 1']) + list(trn_table['Strain 2']))

    for i in trn_table.index:
        if trn_table.loc[i, 'Genome 1 size'] > trn_table.loc[i, 'Genome 2 size']:
            trn_table.loc[i, 'Larger genome'] = trn_table.loc[i, 'Genome 1 size']
        else:
            trn_table.loc[i, 'Larger genome'] = trn_table.loc[i, 'Genome 2 size']

    predict_df = pd.DataFrame(columns=['constant', 'Genome_size'])
    predict_df['Genome_size'] = trn_table['Larger genome'] / 1e6
    predict_df['constant'] = 1

    trn_table['Negative selection cutoff'] = linear_model.get_prediction(predict_df).summary_frame(alpha=0.1)['obs_ci_upper']

    # Filter negative selection cutoff
    if single_cell:  # Special filtering just for single cell genomes
        neg_cutoff = max(trn_table['Negative selection cutoff'])
    else:  # Otherwise just use the negative selection index
         neg_cutoff = trn_table['Negative selection cutoff']

    trn_table = trn_table[trn_table['SSD 95 CI low'] > neg_cutoff]

    # Find clonal clusters
    clonal_df = trn_table[trn_table['Initial divergence'] < clonal_cutoff][['Strain 1', 'Strain 2']]
    print(clonal_df)
    clones = nx.from_pandas_dataframe(clonal_df,
                                      'Strain 1',
                                      'Strain 2')
    clonal_components = tuple(nx.connected_component_subgraphs(clones))
    flat_clones = [node for cluster in clonal_components for node in cluster.nodes()]

    # Get edges between clonal clusters
    for c1, c2 in itertools.combinations(clonal_components, 2):
        n1 = c1.nodes()
        n2 = c2.nodes()
        clonal_edges_1 = trn_table[(trn_table['Strain 1'].isin(n1) & trn_table['Strain 2'].isin(n2))]
        clonal_edges_2 = trn_table[(trn_table['Strain 1'].isin(n2) & trn_table['Strain 2'].isin(n1))]
        if len(clonal_edges_1) > 0 or len(clonal_edges_2) > 0:
            ave_edge = np.average(list(clonal_edges_1['Observed SSD']) + list(clonal_edges_2['Observed SSD']))
            final_edge = (clonal_components.index(c1), clonal_components.index(c2), ave_edge)
            final_edges.append(final_edge)

    # Get average edge weight between each clonal cluster and all non-clonal nodes
    for c in clonal_components:
        n = c.nodes()
        non_clonal = trn_table[(trn_table['Strain 1'].isin(n) & ~trn_table['Strain 2'].isin(flat_clones)) |
                               (~trn_table['Strain 1'].isin(flat_clones) & trn_table['Strain 2'].isin(n))]
        nodes_to_connect = set(list(non_clonal[~non_clonal['Strain 1'].isin(flat_clones)]['Strain 1']) +
                               list(non_clonal[~non_clonal['Strain 2'].isin(flat_clones)]['Strain 2']))

        for node in nodes_to_connect:
            node_to_clone = non_clonal[non_clonal['Strain 1'] == node]
            clone_to_node = non_clonal[non_clonal['Strain 2'] == node]
            ave_edge = np.average(list(node_to_clone['Observed SSD']) + list(clone_to_node['Observed SSD']))
            final_edge = (node, clonal_components.index(c), ave_edge)
            final_edges.append(final_edge)

    # Filter out edges within clonal cluster from network
    keep_edges = trn_table[(~trn_table['Strain 1'].isin(flat_clones)) &
                           (~trn_table['Strain 2'].isin(flat_clones))]
    keep_edges = zip(keep_edges['Strain 1'], keep_edges['Strain 2'], keep_edges['Observed SSD'])
    final_edges += keep_edges
    non_singleton = []

    # Write out the final edgefile
    with open(outfile_name, 'w') as outfile:
        cc = []
        for n1, n2, e in final_edges:
            if type(n1) == int:
                cc.append(n1)
                for c in clonal_components[n1]:
                    non_singleton.append(c)
                n1 = ','.join(clonal_components[n1])
            else:
                non_singleton.append(n1)

            if type(n2) == int:
                cc.append(n2)
                for c in clonal_components[n2]:
                    non_singleton.append(c)
                n2 = ','.join(clonal_components[n2])
            else:
                non_singleton.append(n2)

            outfile.write('\t'.join([str(n1), str(n2), str(e) + '\n']))

        for cindex in set(list(range(0, len(clonal_components)))) - set(cc):
            for c in clonal_components[cindex]:
                non_singleton.append(c)
            n = ','.join(clonal_components[cindex])
            outfile.write('\t'.join([str(n), str(n), '0\n']))

        for n in set(all_strains) - set(non_singleton):
            outfile.write('\t'.join([str(n), str(n), '0\n']))


if __name__ == '__main__':
    main()
