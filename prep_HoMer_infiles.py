#!/usr/bin/python3

import sys
import argparse
from collections import defaultdict
import gzip
import os
from functions import write_fasta, orig_to_nonalphanumeric_map

# Get FASTA files of raw sequences for all core genes in panaroo output.
def main():

    parser = argparse.ArgumentParser(

        description="Reads in panaroo gene_data.csv and gene_presence_absence.csv files and output raw FASTA per non-singleton gene. "
                    "Note that \'refound\' genes are not considered as present. Gene order on scaffolds is taken from the gene_data.csv file. "
                    "Last, all non-alphanumeric characters will be removed from genome, scaffold/contig, and gene names.",

        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-d", "--data", metavar="gene_data.csv.gz", type=str,
                        required=True, help="Path to gene_data.csv.gz")

    parser.add_argument("-p", "--presence", metavar="gene_presence_absence.csv.gz", type=str,
                        required=True, help="Path to gene_presence_absence.csv.gz")

    parser.add_argument("-c", "--chrom_to_consider", metavar="PATH", type=str,
                        required=False, default = None,
                        help="Path to file with list of chromosomes to consider (one per line). All other chromosomes/scaffolds will be ignored. All chromosomes will be considered if this is not set.")

    parser.add_argument("-f", "--fasta", metavar="OUTPUT", type=str,
                        required=True, help="Path to folder to output FASTAs")

    parser.add_argument("-s", "--synteny", metavar="OUTPUT", type=str,
                        required=True, help="Path to output folder for synteny files")

    parser.add_argument("-m", "--mapfiles", metavar="OUTPUT", type=str,
                        required=True, help="Path to output folder for mapping of original ids to non-alphanumeric ids (in case this is needed later)")

    parser.add_argument("-r", "--rare", metavar="OUTPUT", type=str,
                        required=True, help="Path to output file for singleton gene families to ignore")

    args = parser.parse_args()

    # Create output folders.
    if not os.path.exists(args.fasta):
        os.makedirs(args.fasta)

    if not os.path.exists(args.synteny):
        os.makedirs(args.synteny)

    if not os.path.exists(args.mapfiles):
        os.makedirs(args.mapfiles)

    gene_seqs = dict()
    genome_to_scaffolds = defaultdict(set)
    genes_to_genome = dict()
    col_map = dict()
    header_flag = True
    unique_genomes = set()
    scaffolds_to_gene_number = defaultdict(dict)

    genes_to_ignore = set()
    chrom_to_consider = set()
    if args.chrom_to_consider is not None:
        with open(args.chrom_to_consider, 'r') as chrom_fh:
            for chrom_line in chrom_fh:
                chrom_to_consider.add(chrom_line.rstrip())

    with gzip.open(args.data, 'rt') as data_fh:
        for data_line in data_fh:
            data_split = data_line.split(',')

            if header_flag:
                for i in range(len(data_split)):
                    col_map[data_split[i]] = i
                header_flag = False
                continue

            genome_i = col_map['gff_file']
            annot_i = col_map['annotation_id']
            dna_i = col_map['dna_sequence']
            scaffold_i = col_map['scaffold_name']
            clustering_id_i = col_map['clustering_id']

            if '_refound' in data_split[annot_i]:
                continue

            # Skip if this gene is on a chromosome not in the list of chromosomes to consider.
            # And add this gene to the list of genes to ignore.
            if args.chrom_to_consider is not None:
                if data_split[scaffold_i] not in chrom_to_consider:
                    genes_to_ignore.add(data_split[annot_i])
                    continue

            gene_seqs[data_split[annot_i]] = data_split[dna_i]
            genes_to_genome[data_split[annot_i]] = data_split[genome_i]
            genome_to_scaffolds[data_split[genome_i]].add(data_split[scaffold_i])
            unique_genomes.add(data_split[genome_i])
            clustering_id_info = data_split[clustering_id_i].split('_')
            scaffold_gene_number = int(clustering_id_info[-1]) + 1
            scaffolds_to_gene_number[data_split[scaffold_i]][scaffold_gene_number] = data_split[annot_i]

    # Read through panaroo file once to get all unique gene family ids, so that
    # they can be mapped to non-alphanumeric versions, just in case.
    unique_gene_families = set()

    with gzip.open(args.presence, 'rt') as presence_fh:
        next(presence_fh)
        for presence_line in presence_fh:
            presence_line = presence_line.rstrip()
            presence_split = presence_line.split(',')
            unique_gene_families.add(presence_split[0])

    # Get mapping of orig genome ids to clean genome ids (i.e., with all non-alphanumeric values removed).
    # Do the same procedure for gene families and genes too.
    orig_to_clean_genomes = orig_to_nonalphanumeric_map(sorted(list(unique_genomes)))
    orig_to_clean_gene_families = orig_to_nonalphanumeric_map(sorted(list(unique_gene_families)))
    orig_to_clean_genes = orig_to_nonalphanumeric_map(sorted(list(genes_to_genome.keys())))

    # Scaffolds/contigs are actually an exceptional case, as these must be formatted as '_contig__#'
    # in the synteny file (this is taken care of when writing out the synteny files below).
    rare_fh = open(args.rare, 'w')

    genes_to_gene_families = dict()

    with gzip.open(args.presence, 'rt') as presence_fh:
        next(presence_fh)
        for presence_line in presence_fh:
            presence_line = presence_line.rstrip()
            presence_split = presence_line.split(',')
            name = presence_split[0]
            genes_per_genome = presence_split[3:]

            gene_family = defaultdict(set)
            for gene_set in genes_per_genome:
                for g in gene_set.split(';'):
                    if g != '' and '_refound' not in g:
                        if g.endswith('_len'):
                            g = g[:-4]
                        if g not in genes_to_ignore:
                            genome_hit = genes_to_genome[g]
                            gene_family[genome_hit].add(g)

            if len(gene_family.keys()) == 0:
                continue

            gene_family_seqs = dict()
            if len(gene_family.keys()) > 1:
                for genome in gene_family.keys():
                    for g in gene_family[genome]:
                        gene_family_seqs[orig_to_clean_genomes[genome] + '_' + orig_to_clean_genes[g]] = gene_seqs[g]
                        genes_to_gene_families[orig_to_clean_genes[g]] = orig_to_clean_gene_families[name]

                write_fasta(seq = gene_family_seqs,
                            outfile = args.fasta + '/' + orig_to_clean_gene_families[name] + '.fna')

            elif len(gene_family.keys()) == 1:
                print(orig_to_clean_gene_families[name] + '>' + orig_to_clean_genes[list(list(gene_family.values())[0])[0]],
                      file = rare_fh)

    rare_fh.close()

    orig_to_clean_scaffolds = dict()

    # Then print out synteny data (which is taken from the gene_data.csv file).
    for genome_orig in sorted(orig_to_clean_genomes.keys()):
        genome_clean = orig_to_clean_genomes[genome_orig]

        orig_to_clean_scaffolds[genome_clean] = dict()

        genome_synteny_outfile = args.synteny + '/' + genome_clean + '.synteny'
        with open(genome_synteny_outfile, 'w') as synteny_out_fh:
            scaffold_num = 1
            for orig_scaffold_id in sorted(genome_to_scaffolds[genome_orig]):
                new_scaffold_id = '_contig_' + str(scaffold_num)
                scaffold_num += 1

                orig_to_clean_scaffolds[genome_clean][orig_scaffold_id] = new_scaffold_id

                scaffold_out_line = []

                for i in range(1, max(scaffolds_to_gene_number[orig_scaffold_id].keys()) + 1, 1):

                    if i not in scaffolds_to_gene_number[orig_scaffold_id].keys():
                        print('This gene number not found: ' + str(i) + ' for scaffold ' + orig_scaffold_id,
                              file = sys.stderr)
                        continue

                    orig_gene_id = scaffolds_to_gene_number[orig_scaffold_id][i]

                    clean_gene_id = orig_to_clean_genes[orig_gene_id]

                    if clean_gene_id not in genes_to_gene_families.keys():
                        #print('No gene family found for this gene: ' + clean_gene_id,
                        #      file = sys.stderr)
                        continue

                    clean_gene_family_id = genes_to_gene_families[clean_gene_id]

                    scaffold_out_line.append(':'.join([genome_clean, clean_gene_id, new_scaffold_id, clean_gene_family_id]))

                if len(scaffold_out_line) > 0:
                    print('\t'.join(scaffold_out_line), file = synteny_out_fh)

    # Finally, print out mapfiles of original to new ids.
    with open(args.mapfiles + '/genome_ids.tsv', 'w') as genome_map_fh:
        for orig_genome_id in sorted(orig_to_clean_genomes.keys()):
            print(orig_genome_id + '\t' + orig_to_clean_genomes[orig_genome_id],
                  file = genome_map_fh)

    with open(args.mapfiles + '/gene.family_ids.tsv', 'w') as gene_family_map_fh:
        for orig_gene_family_id in sorted(orig_to_clean_gene_families.keys()):
            print(orig_gene_family_id + '\t' + orig_to_clean_gene_families[orig_gene_family_id],
                  file = gene_family_map_fh)

    with open(args.mapfiles + '/gene_ids.tsv', 'w') as gene_map_fh:
        for orig_gene_id in sorted(orig_to_clean_genes.keys()):
            print(orig_gene_id + '\t' + orig_to_clean_genes[orig_gene_id],
                  file = gene_map_fh)

    with open(args.mapfiles + '/scaffold_ids.tsv', 'w') as scaffold_map_fh:
        for genome_id in sorted(orig_to_clean_scaffolds.keys()):
            for scaffold_id in sorted(orig_to_clean_scaffolds[genome_id].keys()):
                print(genome_id + '\t' + scaffold_id + '\t' + orig_to_clean_scaffolds[genome_id][scaffold_id],
                      file = scaffold_map_fh)


if __name__ == '__main__':
    main()
