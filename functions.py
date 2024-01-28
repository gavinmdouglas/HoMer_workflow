#!/usr/bin/python3

import re
import textwrap

def read_fasta(filename, cut_header=False, convert_upper=False):
    '''Read in FASTA file (gzipped or not) and return dictionary with each
    independent sequence id as a key and the corresponding sequence string as
    the value.'''

    # Intitialize empty dict.
    seq = {}

    # Intitialize undefined str variable to contain the most recently parsed
    # header name.
    name = None

    def parse_fasta_lines(file_handle):
        for line in file_handle:

            line = line.rstrip()

            if len(line) == 0:
                continue

            # If header-line then split by whitespace, take the first element,
            # and define the sequence name as everything after the ">".
            if line[0] == ">":

                if cut_header:
                    name = line.split()[0][1:]
                else:
                    name = line[1:]

                name = name.rstrip("\r\n")

                # Make sure that sequence id is not already in dictionary.
                if name in seq:
                    sys.stderr("Stopping due to duplicated id in file: " + name)

                # Intitialize empty sequence with this id.
                seq[name] = ""

            else:
                # Remove line terminator/newline characters.
                line = line.rstrip("\r\n")

                # Add sequence to dictionary.
                seq[name] += line

    # Read in FASTA line-by-line.
    if filename[-3:] == ".gz":
        with gzip.open(filename, "rt") as fasta_in:
            parse_fasta_lines(fasta_in)
    else:
        with open(filename, "r") as fasta_in:
            parse_fasta_lines(fasta_in)

    if convert_upper:
        for seq_name in seq.keys():
            seq[seq_name] = seq[seq_name].upper()

    return seq


def write_fasta(seq, outfile):
    out_fasta = open(outfile, "w")

    # Look through sequence ids (sorted alphabetically so output file is
    # reproducible).
    for s in sorted(seq.keys()):
        out_fasta.write(">" + s + "\n")
        out_fasta.write(textwrap.fill(seq[s], width=70) + "\n")

    out_fasta.close()


def orig_to_nonalphanumeric_map(input_ids):
    '''Given a list of ids, return a dictionary mapping these ids to the new ids after removing all
    non-alphanumeric values, while ensuring they are still unique ids.'''
    orig_to_clean = dict()
    for orig_id in input_ids:
        clean_id = re.sub(r'\W+', '', orig_id)
        clean_id = re.sub(r'_+', '', clean_id)
        clean_id_final = clean_id
        dup_num = 1
        while clean_id_final in orig_to_clean.keys():
            clean_id_final = clean_id + str(dup_num)
            dup_num += 1
        orig_to_clean[orig_id] = clean_id_final

    return orig_to_clean