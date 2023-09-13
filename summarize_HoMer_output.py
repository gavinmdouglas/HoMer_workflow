#!/usr/bin/python3

import argparse
import re
import ete3
import os
import sys


def main():

    parser = argparse.ArgumentParser(

        description="Parse HoMer output files to get cleaned-up table of all inferred HGT events (including single gene inferences that were not included in any multi-gene transfer event).",

        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-a', '--agg_out', metavar="INPUT", type=str,
                        required=True, help="Path to HoMer aggregation step output folder")

    parser.add_argument('-m', '--map_folder', metavar="INPUT", type=str,
                        required=True,
                        help="Path to folder with mapfiles breaking down original gene and genome ids prior to simplication for HoMer.")

    parser.add_argument('-i', '--input', metavar="INPUT", type=str,
                        required=True, help="Path to output file for main HoMer step")

    parser.add_argument('-s', '--species_tree', metavar="TREE", type=str,
                        required=True, help="Path to species tree (with labelled nodes)")

    parser.add_argument('-o', '--output_folder', metavar="PATH", type=str,
                        required=True, help="Path to output folder")

    args = parser.parse_args()

    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)

    # Read through mapfiles and get mapping of ids from original to HoMer-prepped ids.
    genefamily_id_map = dict()
    genome_id_map = dict()

    with open(args.map_folder + '/gene.family_ids.tsv', 'r') as genefamily_map_fh:
        for genefamily_map_line in genefamily_map_fh:
            genefamily_map_line_split = genefamily_map_line.split()
            genefamily_id_map[genefamily_map_line_split[1]] = genefamily_map_line_split[0]

    with open(args.map_folder + '/genome_ids.tsv', 'r') as genome_map_fh:
        for genome_map_line in genome_map_fh:
            genome_map_line_split = genome_map_line.split()
            genome_id_map[genome_map_line_split[1]] = genome_map_line_split[0]

    # Write out mapping of all internal nodes to all *tip* children (i.e., leaves).
    tree = ete3.Tree(args.species_tree, format=1)
    with open(args.output_folder + '/nodes_to_leaves.tsv', 'w') as node_to_leaf_fh:
        for node in tree.traverse('postorder'):
            node_leaves = []
            children = node.children
            if len(children) == 0:
                continue
            for leaf in children[0]:
                node_leaves.append(leaf.name)
            for leaf in children[1]:
                node_leaves.append(leaf.name)
            print(node.name + '\t' + ';'.join(sorted(node_leaves)),
                  file=node_to_leaf_fh)

            # Also add this node to be a "genome" id, to
            # make it easier to fix the genome ids.
            genome_id_map[node.name] = node.name

    # Read through all reconciliation hits and keep track of all individual
    # transfer inferences.
    ranger_outfiles = [file for file in os.listdir(args.agg_out) if os.path.isfile(os.path.join(args.agg_out, file))]

    gene_family_transfer_nodes = dict()

    for outfile in ranger_outfiles:
        gene_family = genefamily_id_map[os.path.basename(outfile).replace('.tree', '')]
        with open(args.agg_out + '/' + outfile, 'r') as outfile_fh:
            for line in outfile_fh:

                # Only consider lines containing information about specific nodes.
                # And only when 100 of replicates identified it to represent a transfer
                # event.
                if re.match(r'^m\d+ = ', line) and 'Transfers = 100' in line:
                    # Parse out key parts of line:
                    pattern = r'(m\d+) = LCA\[\w+_\w+, \w+_\w+\]: \[Speciations = \d+, Duplications = \d+, Transfers = \d+\], \[Most Frequent mapping --\> (\w+), (\d+) times\], \[Most Frequent recipient --\> (\w+), (\d+) times\]\.'

                    match = re.search(pattern, line)
                    gene_tree_node = match.group(1)
                    most_frequent_mapping = genome_id_map[match.group(2)]
                    most_frequent_mapping_count = match.group(3)
                    most_frequent_recipient = genome_id_map[match.group(4)]
                    most_frequent_recipient_count = match.group(5)

                    single_transfer_id = ';'.join([gene_family,
                                                   most_frequent_mapping,
                                                   most_frequent_mapping_count,
                                                   most_frequent_recipient,
                                                   most_frequent_recipient_count])

                    gene_family_transfer_nodes[single_transfer_id] = gene_tree_node

    # Then read through final HoMer output table containing
    # horizontal multi-gene transfers (HMGTs).
    hmgt_num = 0
    hmgt_map = {}

    with open(args.input, 'r') as homer_fh:
        for homer_line in homer_fh:
            if homer_line.startswith('Donor'):

                # Initialize all variables, to help avoid bugs if lines are missing unexpectedly.
                donor = None
                recipient = None
                hmgt_gene_families = []

                genome_matches = re.search(r'Donor: (\w+), Recipient: (\w+)$', homer_line)
                donor = genome_id_map[genome_matches.group(1)]
                recipient = genome_id_map[genome_matches.group(2)]

            elif homer_line.startswith('Contig'):
                hmgt_gene_families.append(genefamily_id_map[re.search(r'Gene Family: (\w+)$', homer_line).group(1)])

            elif len(homer_line.rstrip()) == 0 and len(hmgt_gene_families) > 0:
                hmgt_num += 1
                for gene_family in hmgt_gene_families:
                    transfer_id = ';'.join([gene_family, donor, recipient])
                    hmgt_map[transfer_id] = 'hmgt' + str(hmgt_num)
                hmgt_gene_families = []

    # Finally, loop through individual transfer events and create output table.
    with open(args.output_folder + '/transfers.tsv', 'w') as transfer_out_fh:
        print('\t'.join(['gene.family',
                         'most.freq.donor', 'most.freq.donor.instances',
                         'most.freq.recipient', 'most.freq.recipient.instances',
                         'gene_tree_node', 'hgt_instance']),
              file=transfer_out_fh)
        for transfer_identifier in sorted(gene_family_transfer_nodes.keys()):
            outline = transfer_identifier.split(';')
            outline.append(gene_family_transfer_nodes[transfer_identifier])
            if transfer_identifier in hmgt_map.keys():
                outline.append(hmgt_map[transfer_identifier])
            else:
                outline.append('single')
            print('\t'.join(outline), file=transfer_out_fh)


if __name__ == '__main__':
    main()
