#!/usr/bin/python3

import argparse
from collections import defaultdict
import sys

import pprint

def main():

    parser = argparse.ArgumentParser(

        description="Although focusing on the actual internal nodes where transfer are "
                    "predicted to have occurred is valuable, sometimes information on the "
                    "tips of the trees sharing similar genes due to HGT is needed "
                    "(e.g., to compare with tools that only call HGT between extant genomes). "
                    "This script parses the transfer and node label tables produced by "
                    "summarize_HoMer_output.py and outputs 'effective' transfer partners, "
                    "which includes combinations of genomes descended from donor/recipient nodes."
                    "This output will also include direct donor and recipient pairs (which will be distinguished).",

        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-t', '--transfers', metavar="INPUT", type=str,
                        required=True, help="Path to summary of all HoMer/RANGER-DTL transfer events.")

    parser.add_argument('-n', '--nodes_map', metavar="INPUT", type=str,
                        required=True,
                        help="Path to table containing mapping of all internal node ids to descendant tips.")

    args = parser.parse_args()

    # Read in node mapfile and create dictionary mapping internal node ids to
    # all tips descending from that node.
    node_descendants = dict()
    with open(args.nodes_map, 'r') as node_map_fh:
       for node_map_line in node_map_fh:
            node_map_line_split = node_map_line.split()
            if len(node_map_line_split) != 2:
                sys.exit('Stopping - node to leaves mapfile should only have two columns.')
            node_descendants[node_map_line_split[0]] = set(node_map_line_split[1].split(';'))

    # Read through transfer table, and get all transfer partners.
    with open(args.transfers, 'r') as genome_map_fh:
        # Skip header-line.
        genome_map_fh.readline()

        # Keep track of all transfers by each gene family.
        transfers_by_genes = defaultdict(list)
        for genome_map_line in genome_map_fh:
            genome_map_line_split = genome_map_line.split()
            transfers_by_genes[genome_map_line_split[0]].append(genome_map_line_split[1:])
        

    # Then loop through each gene family, and figure out:
    
    # 1) If there are any duplicate donor/recipient pairs (which will be ignored),
    # or if there are instances where the donor and recipient are the same genome
    # (which will also be ignored). In both cases, these likely represent issues, either
    # with the upstream steps or that filtering is needed, which will be written out to
    # standard error.
 
    # 2) How many times each recipient is listed (clearly less confidence if a tip or node
    # is listed as a recipient from different donors). Will treat these as 'effective' transfers if so.

    # 3) For all non-tip recipients (i.e., internal nodes), need to figure out nearest node to tips (e.g., closest to tips).
    # This is because only the most recent transfer leading to the extant tips should be considered
    # (including if the transfer is to the tip itself, in which case earlier transfers to the lineage leading to the tip).

    print("gene_family\teffective_donor\teffective_recipient\tpair_type\tcalled_donor\tcalled_donor_tally\tcalled_recipient\tcalled_recipient_tally\tgene_tree_node\thgt_instance")

    for gene_family in transfers_by_genes.keys():
        raw_transfer_info = transfers_by_genes[gene_family]
        filt_transfer_info = []

        # First, check for duplicate donor/recipient pairs.
        # Also, check for how many times each recipient is listed.
        all_transfer_pairs = set()
        recipient_observed = set()
        duplicated_recipients = set()
        
        # Also, keep track of all internal nodes marked as recipients, and keep track of the descendant tips.
        node_recipient_descendants = dict()
        lowest_recipient_node = dict()

        for transfer_info in raw_transfer_info:
            donor = transfer_info[0]
            recipient = transfer_info[2]

            if recipient in recipient_observed:
                duplicated_recipients.add(recipient)
            else:
                recipient_observed.add(recipient)

            # If recipient is an internal node, add to dictionary.
            # Otherwise, add tip to dictionary keeping track of lowest
            # recipient node per genome (which would be the tip itself).
            if recipient in node_descendants.keys():
                node_recipient_descendants[recipient] = node_descendants[recipient]
            else:
                lowest_recipient_node[recipient] = recipient

            if donor == recipient:
                sys.stderr.write("WARNING: Identical label for donor and recipient (likely due to low reproducibility when aggregating RANGER-DTL output) for " + gene_family + ": " + donor + "\n")
                continue

            transfer_pair = donor + ' --> ' + recipient
            if transfer_pair in all_transfer_pairs:
                sys.stderr.write("WARNING: Duplicate donor/recipient pair found for gene family " + gene_family + ": " + transfer_pair + "\n")
                continue
            
            # Otherwise add to new list.
            filt_transfer_info.append(transfer_info)
        
        if len(duplicated_recipients) > 0:
            sys.stderr.write("WARNING: The same recipient was labelled for multiple transfer events for " + gene_family + ", " + str(len(duplicated_recipients)) + " separate recipient(s)\n")

        # Then, for all recipient nodes, figure out lowest recipient node/tip (which is trivial for the genome recipients)
        # by looping through the nodes from smallest to largest.
        sorted_recipient_nodes = sorted(node_recipient_descendants,
                                        key=lambda x: len(node_recipient_descendants[x]))

        for recipient_node in sorted_recipient_nodes:
            for descendant in node_recipient_descendants[recipient_node]:

                if descendant not in lowest_recipient_node.keys():
                    lowest_recipient_node[descendant] = recipient_node
                
        # Then, loop through all transfer pairs, and write out all 'effective' transfers.
        for transfer_info in filt_transfer_info:
            donor = transfer_info[0]
            donor_tally = transfer_info[1]
            recipient = transfer_info[2]
            recipient_tally = transfer_info[3]
            gene_tree_node  = transfer_info[4]
            hgt_instance = transfer_info[5]

            if recipient in node_descendants.keys():
                effective_recipients = list(node_descendants[recipient])
            else:
                effective_recipients = [recipient]

            if donor in node_descendants.keys():
                effective_donors = list(node_descendants[donor])
            else:
                effective_donors = [donor]

            # Identify recipient and donor genomes that should be included.
            if len(effective_recipients) > 1:
                tmp_effective_recipients = []
                for effective_recipient in effective_recipients:
                    if lowest_recipient_node[effective_recipient] == recipient:
                        tmp_effective_recipients.append(effective_recipient)
                effective_recipients = tmp_effective_recipients

            num_raw_effective_donors = len(effective_donors)

            if num_raw_effective_donors > 1:
                tmp_effective_donors = []
                for effective_donor in effective_donors:
                    if effective_donor in lowest_recipient_node.keys():
                        if lowest_recipient_node[effective_donor] in node_descendants.keys():
                            effective_donor_recipient_node_size = len(node_descendants[lowest_recipient_node[effective_donor]])
                        else:
                            effective_donor_recipient_node_size = 1

                        # Include as a donor if the HGT event for this donor occurred in an earlier (larger) node.
                        if effective_donor_recipient_node_size > num_raw_effective_donors:
                            tmp_effective_donors.append(effective_donor)
                        elif effective_donor_recipient_node_size == num_raw_effective_donors:
                            print('WARNING: Same node for genome ancestral node that was donor and lowest common recipient node for recipient ' + recipient + ' and donor ' + effective_donor + ' for gene family ' + gene_family + '.')
                    else:
                        tmp_effective_donors.append(effective_donor)
                effective_donors = tmp_effective_donors

            # Then, write out all effective donor/recipient pairs.
            if len(effective_recipients) == 1 and len(effective_donors) == 1:
                pair_type = 'direct'
            else:
                pair_type = 'indirect'

            for effective_donor in effective_donors:
                for effective_recipient in effective_recipients:
                    print(gene_family + '\t' + effective_donor + '\t' + effective_recipient + '\t' + pair_type + '\t' + donor + '\t' + donor_tally + '\t' + recipient + '\t' + recipient_tally + '\t' + gene_tree_node + '\t' + hgt_instance)


if __name__ == '__main__':
    main()
