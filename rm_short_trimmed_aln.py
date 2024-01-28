#!/usr/bin/python3

import sys
import argparse
import os
from functions import read_fasta

def main():

    parser = argparse.ArgumentParser(

        description="Identified gene families that should be ignored, due to wonky alignments. "
                     "Reads in per-gene family alignment (after being trimmed) . "
                     "Will delete gene family alignment FASTAs below are certain length, and will add genes in this file to the species' rare_genes.txt.",

        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-f", "--fasta", metavar="FASTA", type=str,
                        required=True, help="Trimmed gene family FASTA to parse (and potentially delete).")

    parser.add_argument("-r", "--rare", metavar="OUTPUT", type=str,
                        required=True, help="Path to species' rare_genes.txt file.")

    parser.add_argument("-l", "--length", metavar="LENGTH", type=int,
                        required=False, default=200, help="Minimum length of gene family alignment to keep.")

    args = parser.parse_args()

    # Read in gene family alignment FASTA.
    seqs = read_fasta(args.fasta)

    # Check that all sequences are the same length, or throw error.
    seq_lengths = set()
    for seq in seqs.values():
        seq_lengths.add(len(seq))

    if len(seq_lengths) > 1:
        print('Sequences are not all the same length, as expected for an alignment. Exiting.', file=sys.stderr)
        sys.exit(1)

    # If all sequences are the same length, then check if they are below the minimum length.
    # If so, then delete the gene family alignment FASTA and add the genes to the rare_genes.txt file.
    seq_length = list(seq_lengths)[0]
    if seq_length < args.length:
        os.remove(args.fasta)
        with open(args.rare, 'a') as rare_fh:
            for seq_id in seqs.keys():
                print(seq_id, file=rare_fh)

    # Print to STDERR that this file was removed.
    print('Removed', args.fasta, 'from analysis.', file=sys.stderr)


if __name__ == '__main__':
    main()
