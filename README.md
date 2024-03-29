# HoMer_workflow

This repository contains the commands and scripts for running HoMer ([website](https://compbio.engr.uconn.edu/software/homer/), [paper](https://academic.oup.com/mbe/article/38/6/2639/6132264)) and RANGER-DTL ([website](https://compbio.engr.uconn.edu/software/RANGER-DTL/), [paper](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty314/4983062?guestAccessKey=f95d29d3-0976-4c4d-b2ea-139dabd24bf8)). RANGER-DTL reconciles gene trees with species (or strain) trees to identify potential duplication, transfer, or loss events that occurred. HoMer builds upon this output, along with additional genome annotation and synteny information, to identify putative multi-gene transfer events.

**Please note that this repository is _not_ an official workflow provided by the tool developers**. This is simply the workflow I came up with for running HoMer, which required troubleshooting and custom scripts to produce the final results. My intention here is to provide a common repository that can be cited for multiple projects worked on by myself and colleagues, but this may be useful to others as well.

_This workflow was developed based on HoMer v1.0 and RANGER-DTL v2.0._

Before starting the workflow you should keep in mind is that HoMer does not allow *any* non-alphanumeric characters (e.g., underscores) in any input file paths. This is true for both the file names themselves *and for the full file paths* (e.g., no underscores or periods can be present in any folder names that are part of file paths). This is important as it can be inconvenient (and confusing if you're unaware of what is causing problems) if you have to rename files and move them to new folders once you get to that step (below I made symbolic links to where the files were to get around this problem).


## Prep input for HoMer (and RANGER-DTL) workflow

The `prep_HoMer_infiles.py` script will parse [panaroo](https://github.com/gtonkinhill/panaroo) output FASTA files for each non-singleton gene. A file containing all singleton genes to be ignored will also be created. Note that 'refound' genes are not considered as present. Also, note that the gene ordering in each scaffold is taken from the `clustering_id` in the `gene_data.csv.gz` file. Last, all non-alphanumeric characters will be removed from genome, scaffold/contig, and gene names

```
mkdir homer_prep

python prep_HoMer_infiles.py \
            -d panaroo_out/$SP/gene_data.csv.gz \
            -p panaroo_out/gene_presence_absence.csv.gz \
            -f homer_prep/fastas_out \
            -s homer_prep/synteny_out \
            -m homer_prep/map_out \
            -r homer_prep/rare_genes.txt
```

These gene seqeuences must then be aligned, which I did with [muscle](https://www.drive5.com/muscle/).

To get muscle commands written to a file to run in parallel:
```
mkdir homer_prep/fastas_out_aligned

for FASTA_PATH in homer_prep/fastas_out/*fna; do
	FASTA=$( basename $FASTA_PATH )
	echo "muscle -in $FASTA_PATH -out homer_prep/fastas_out_aligned/$FASTA -seqtype dna" >> muscle_cmds.sh
done
```

To run these commands in parallel (with [GNU parallel](https://www.gnu.org/software/parallel/), note that you can change `NUM_CORES`):
```
NUM_CORES=1
cat muscle_cmds.sh | parallel -j $NUM_CORES --eta --joblog muscle_cmds.log '{}'
```

We wrote the commands to run in parallel to a file, rather than running them directly, as I wrote a convenience script (`gnu.parallel_cmds_vs_log.py`, [available here](https://github.com/gavinmdouglas/parallel_joblog_summary)) for comparing parallel's joblog output to the original set of commands to run, to ensure that all commands completed correctly.

```
python gnu.parallel_cmds_vs_log.py \
			--cmds muscle_cmds.sh \
			--log muscle_cmds.log \
			--failed muscle_failed_cmds.sh
```

The above command should let you know that all jobs completed correctly. If not, you have some troubleshooting to do!

Note that when I ran muscle on my dataset, I found that some genes were really wonky and didn't align properly (35 gene families in total). I decided to add these to the 'rare' genes sets, so that they would be ignored:
```
while read FAILED_CMD; do
  SPECIES=$( python -c "cmd=\"$FAILED_CMD\"; cmd_split=cmd.split(); print(cmd_split[2].split('/')[1])" )
  GENE_FAMILY_FASTA=$( python -c "cmd=\"$FAILED_CMD\"; cmd_split=cmd.split(); print(cmd_split[2].split('/')[3])" )
  GENE_FAMILY=$( basename $GENE_FAMILY_FASTA .fna )
  FASTA_PATH=$( python -c "cmd=\"$FAILED_CMD\"; cmd_split=cmd.split(); print(cmd_split[2])" )
  for GENE in $( grep ">" $FASTA_PATH ); do
	echo "$GENE_FAMILY""$GENE" >> homer_prep/$SPECIES/rare_genes.txt
  done
done < muscle_failed_cmds.sh
```

## Trim alignments and remove genes with short trimmed alignments

Trim alignments with [TrimAl](http://trimal.cgenomics.org/). I decided to just remove any columns containing *any* gap characters.

The being that stretches of gaps could throw off the RANGER-DTL inferences (and are more likely to represent mis-alignments).

```
mkdir homer_prep/$SP/fastas_out_aligned_nogaps

for INFASTA in homer_prep/$SP/fastas_out_aligned/*.fna; do
	BASEFASTA=$( basename $INFASTA .fna )
	OUTFASTA="homer_prep/$SP/fastas_out_aligned_nogaps/$BASEFASTA.fna"
	echo "trimal -in $INFASTA -gt 1 -out $OUTFASTA" >> gene_trimal_cmds.sh
done

cat gene_trimal_cmds.sh | parallel -j 50 --eta --joblog gene_trimal_cmds.log '{}'

python /home/gdouglas/local/parallel_joblog_summary/gnu.parallel_cmds_vs_log.py \
	--cmds gene_trimal_cmds.sh \
	--log gene_trimal_cmds.log
```

Then remove any gene family alignments < 200 bp after the above step, and add the individual genes to rare_genes.txt.
Note that this command will likely remove some of the FASTA files created above - take a look at the help documentation.
```
python rm_short_trimmed_aln.py \
	-f homer_prep/$SP/fastas_out_aligned_nogaps/ \
	-r homer_prep/$SP/rare_genes.txt
```

## Build gene trees (with FastTree in this case)
```
mkdir homer_prep/fastas_out_aligned_trees

for FASTA_PATH in homer_prep/$SP/fastas_out_aligned/*fna; do
	FASTABASE=$( basename $FASTA_PATH .fna );
	echo "fasttree -nt -gtr -gamma -quiet $FASTA_PATH >  homer_prep/$SP/fastas_out_aligned_trees/$FASTABASE.tree" >> fasttree_cmds.sh
done

conda activate fasttree
cat fasttree_cmds.sh | parallel -j $NUM_CORES --eta --joblog fasttree_cmds.log '{}'

conda deactivate
python gnu.parallel_cmds_vs_log.py \
	--cmds fasttree_cmds.sh \
	--log fasttree_cmds.log
```

## Post-process trees:

Then process trees to work with HoMer: non-alphanumeric characters must be removed, and then they must be binary + rooted. This is true for the species tree and for each gene tree.

First, for species trees.

```
ORIG_TREE="panaroo_out/core_gene_alignment.aln.treefile"
NEW_TREE="panaroo_out/core_gene_alignment.aln.prepped.treefile"
Rscript simplify.labels_midpoint.root_and_binarize.tree.R $ORIG_TREE $NEW_TREE TRUE TRUE
```

Then for gene trees.
```
mkdir homer_prep/fastas_out_aligned_trees_prepped
for TREE in  homer_prep/fastas_out_aligned_trees/*tree; do
	TREE_OUT=$( basename $TREE )
	echo "Rscript simplify.labels_midpoint.root_and_binarize.tree.R $TREE homer_prep/fastas_out_aligned_trees_prepped/$TREE_OUT FALSE TRUE" >> prep_gene_trees.sh
done

cat prep_gene_trees.sh | parallel -j $NUM_CORES --eta --joblog prep_gene_trees.log '{}'

python gnu.parallel_cmds_vs_log.py \
	--cmds prep_gene_trees.sh \
	--log prep_gene_trees.log
```


## Run HoMer Aggregate

Note that `Homer_Aggregate.py` needs to be run in the source code folder, which is why the command `cd /path/to/HoMer_Linux/HoMer_Aggregate/` is included.

```
mkdir /path/to/output/homer_aggregate_out

conda activate homer

cd /path/to/HoMer_Linux/HoMer_Aggregate/

python ./Homer_Aggregate.py \
	-g /path/to/homer_prep/fastas_out_aligned_trees_prepped/ \
	-s/path/to/panaroo_out/core_gene_alignment.aln.prepped.treefile \
	-f /path/to/output/homeraggregate/
```

## Run main HoMer step
```
mkdir homer_output

cd /home/gdouglas/local/prg/HoMer_Linux/HoMer

python ./HoMer.py -g /path/to/output/homeraggregate/ReconciliationResults/ \
                  -s /path/to/output/homeraggregate/LabeledSpeciesTree.newick \
                  -n /path/to/output/homer_prep/synteny_out \
                  > /path/to/output/homer_output.txt
```

## Summarize output 

This step includes combining both multi-gene and other calls into single tables, i.e., so that RANGER-DTL calls that could not be called as multi-gene transfers into table. This does not mean that these are confident single-gene calls: just that we don't have clear evidence that they were in multi-gene transfer events. **Note**: This script assumes you ran 100 `Homer_Aggregate.py` replicates specifically per gene family (the default).

```
conda deactivate
conda activate ete3

python summarize_HoMer_output.py \
	 -a /path/to/output/homeraggregate/ReconciliationResults/ \
         -i /path/to/output/homer_output.txt \
         -s /path/to/output/homeraggregate/LabeledSpeciesTree.newick \
         --map_folder /path/to/output/homer_prep/map_out \
         --output_folder /path/to/output/homer_rangerdtl_summaries/
```

There will be two output files created.

`nodes_to_leaves` is a simple tab-delimited mapfile of internal node ids in the species tree to all underlying child tips, which are delimited by semi-colons. This could be useful for making sense of transfers in future steps. Note that the tip names correspond to the original genome IDs, i.e., all non-alphanumeric characters removed for running this pipeline will be re-added.

`transfers.tsv` is the simplified table of all HGT calls based on HoMer *and* RANGER-DTL. The table structure looks like this:

```
gene.family     most.freq.donor most.freq.donor.instances       most.freq.recipient     most.freq.recipient.instances   gene_tree_node
  hgt_instance
AAH1    MALA_SAMN05421699_METAG_ABPHNPKL        100     MALA_SAMN05422108_METAG_FKLFPLKL        100     m129    single
AAH1    MALA_SAMN05421699_METAG_ABPHNPKL        100     TARA_SAMEA2619747_METAG_JBDDHCJO        100     m125    single
ARG7~~~argJ     TARA_SAMEA2620783_METAG_DJEACFGL        100     n54     48      m37     single
COQ2~~~ubiA     TARA_SAMEA2620256_METAG_MECLBIAM        34      TARA_SAMEA2621285_METAG_PDIJJACG        56      m57     single
```

The columns correspond to:
* gene.family - the gene family name
* most.freq.donor - name of genome (or internal node in species tree) that was marked as the most *frequent donor*
* most.freq.donor.instances - number of `Homer_Aggregate.py` replicates where this genome was the *donor*.
* most.freq.recipient - name of genome (or internal node in species tree) that was marked as the most *frequent recipient*
* most.freq.recipient.instances - number of `Homer_Aggregate.py` replicates where this genome was the *recipient*.
* gene_tree_node - Id of node in gene tree where transfer was inferred
* hgt_instance - Either `single` or `hmgtN`, where `N` is an integer specifying which horizontal multi-gene transfer event this gene was part of. All HMGT events are numbered sequentially based on their order in the raw HoMer output.

Note that HoMer only calls HMGT events, not single gene transfer events, so here `single` really just means 'not confidentally part of HMGT' rather than 'confidently a single transfer event'.

## Next steps

I suggest you filter `transfers.tsv` to only retain well-supported transfer pairs. I favoured limiting the output to transfer events supported by all 100 replicates, which you can do by only retaining rows of this table where `most.freq.donor.instances` and `most.freq.recipient.instances` equal 100.

Sometimes you might also want to get a simplified breakdown of transfer partners that includes all possible genome pairs, rather than transfers between genomes and internal tree nodes (e.g., `n54` in the above example table for gene family `ARG7~~~argJ`). This can be helpful if comparing to other HGT inference approaches, such as similarity hits between genomes, but note that this simplified breakdown will not provide insight into the numbers of independent HGT events (unlike `transfers.tsv`), and instead represents cases where alleles of two genomes share similarity due to HGT at some point between separate lineages.

You can get this breakdown with the below command (`transfers_filt.tsv` referring to the filtered table after subsetting only to rows with 100% of replicates agreeing on donor and recipient, as described above). Note that this command also requires a gzipped Panaroo-output presence/absence table, which is used to identify cases where genes are missing in genomes (as the gene loss information was not retained in earlier steps).

```
python simplified_HGT_partners.py -t transfers_filt.tsv \
                                  -n nodes_to_leaves.tsv \
                                  -p panaroo_output/gene_presence_absence.csv.gz \
                                  > simplified_pairwise_transfers.tsv
```

The output looks like this:

```
gene_family     effective_donor effective_recipient     pair_type       called_donor    called_donor_tally      called_recipient
        called_recipient_tally  gene_tree_node  hgt_instance
COQ3_2  TARA_SAMEA2622690_METAG_KEIPAIEC        TARA_SAMEA2622678_METAG_EDMLIGDO        direct  TARA_SAMEA2622690_METAG_KEIPAIEC
        100     TARA_SAMEA2622678_METAG_EDMLIGDO        100     m12     single
COQ3_2  TARA_SAMN05326645_METAG_PSE00098        TARA_SAMEA2622673_METAG_DIPODJBH        direct  TARA_SAMN05326645_METAG_PSE00098
        100     TARA_SAMEA2622673_METAG_DIPODJBH        100     m7      single
COQ3_2  TARA_SAMN05326645_METAG_PSE00098        TARA_SAMEA2622736_METAG_PADLMHCI        indirect        n20     100     TARA_SAMEA26227
36_METAG_PADLMHCI        100     m2      single
COQ3_2  TARA_SAMEA2619779_METAG_MPHGPNAE        TARA_SAMEA2622736_METAG_PADLMHCI        indirect        n20     100     TARA_SAMEA26227
36_METAG_PADLMHCI        100     m2      single
COQ3_2  TARA_SAMEA2622197_METAG_ODGMLBOK        TARA_SAMEA2622736_METAG_PADLMHCI        indirect        n20     100     TARA_SAMEA26227
36_METAG_PADLMHCI        100     m2      single
```

The columns are similar to `transfers.tsv`, but the original donor and recipient columns are prefixed with "called" to distinguish the "effective" donors and recipients from the original calls. The term "effective" is used to indicate that some of these pairs are not direct transfer pairs, but are simply the descendants of a lineage involved in an ancestral transfer event, and which did not subsequently receive a different allele (checking for this latter possibility is the key reason why this script does not just output all tips underlying a node involved in a transfer: sometimes those underlying tips are recipients of subsequent transfers).

The new columns are:
* effective_donor/effective_recipient: Genome ID of effective donor/recipient (which can just be the same as the original, if it was already at the genome level, or will be a descendant of a node where a transfer happened).
* pair_type: Either "direct" or "indirect", to distinguish cases where the transfer event was inferred to be directly between two genomes (direct), or whether at least one of the donor or recipient was an internal node (indirect).
  
