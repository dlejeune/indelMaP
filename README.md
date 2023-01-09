# Indel-Aware-Parsimony-Methods

The repository includes the code for multiple sequence alignment (ParsAlign) and ancestral sequence reconstruction (ParsASR) under indel-aware maximum parsimony. 

## Installation and dependencies
The code does not require installation. Download the repository and execute the ParsAlign.py script for multiple sequence alignment and the ParsASR.py script for ancestral sequence reconstruction. The code is executed using python and uses the following libraries: ete3, NumPy, random, math, SciPy, os, and argparse. 

## Usage
For minimal usage, ParsAlign.py requires a sequences file in fasta format and a corresponding guide tree in newick format, as well as the required alphabet (protein or DNA).

```
python ParsAlign.py --seq_file sequenceFile.fasta --tree_file treeFile.nwk --alphabet Protein 	
```

Minimal usage for ParsASR.py is similar, it requires a multiple sequence alignment in fasta format and a corresponding guide tree in newick format, as well as the required alphabet (protein or DNA).

```
python ParsASR.py --msa_file multipleSequenceAlignmentFile.fasta --tree_file treeFile.nwk --alphabet Protein 	
```


