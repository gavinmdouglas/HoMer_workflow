# HoMer_workflow

This repository contains the commands and scripts for running HoMer ([website](https://compbio.engr.uconn.edu/software/homer/), [paper](https://academic.oup.com/mbe/article/38/6/2639/6132264)) and RANGER-DTL ([website](https://compbio.engr.uconn.edu/software/RANGER-DTL/), [paper](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty314/4983062?guestAccessKey=f95d29d3-0976-4c4d-b2ea-139dabd24bf8)). RANGER-DTL reconciles gene trees with species (or strain) trees to identify potential duplication, transfer, or loss events that occurred. HoMer builds upon this output, along with additional genome annotation and synteny information, to identify putative multi-gene transfer events.

**Please note that this repository is _not_ an official workflow provided by the tool developers**. This is simply the workflow I came up with for running HoMer, which required troubleshooting and custom scripts to produce the final results. My intention here is to provide a common repository that can be cited for multiple projects worked on by myself and colleagues, but this may be useful to others as well.

_This workflow was developed based on HoMer v1.0 and RANGER-DTL v2.0._

Before starting the workflow you should keep in mind is that HoMer does not allow *any* non-alphanumeric characters (e.g., underscores) in any input file paths. This is true for both the file names themselves *and for the full file paths* (e.g., no underscores or periods can be present in any folder names that are part of file paths). This is important as it can be inconvenient (and confusing if you're unaware of what is causing problems) if you have to rename files and move them to new folders once you get to that step (below I made symbolic links to where the files were to get around this problem).


## Prep input for HoMer (and RANGER-DTL) workflow

The `prep_HoMer_infiles.py` script will parse [panaroo](https://github.com/gtonkinhill/panaroo) output FASTA files for each non-singleton gene. A file containing all singleton genes to be ignored will also be created. Note that 'refound' genes are not considered as present. Also, note that the gene ordering in the scaffold is assumed to be in the actual gene names themselves. Last, all non-alphanumeric characters will be removed from genome, scaffold/contig, and gene names

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
python /home/gdouglas/local/utils/gnu.parallel_cmds_vs_log.py \
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
python /home/gdouglas/local/utils/gnu.parallel_cmds_vs_log.py \
	--cmds prep_gene_trees.sh \
	--log prep_gene_trees.log
```


## Run HoMer Aggregate

Note that this commands needs to be run in the actual folder where the files are kept.


mkdir /path/to/output/homer_aggregate_out
conda activate homer
cd /path/to/HoMer_Linux/HoMer_Aggregate/
python ./Homer_Aggregate.py \
	-g /path/to/homer_prep/fastas_out_aligned_trees_prepped/ \
	-s/path/to/panaroo_out/core_gene_alignment.aln.prepped.treefile \
	-f /path/to/output/homer_aggregate_out/
