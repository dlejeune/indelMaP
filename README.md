# indelMaP: Indel-Aware Maximum Parsimony

The repository includes the code for multiple sequence alignment (indelMaP_MSA) and ancestral sequence reconstruction (indelMaP_ASR) under indel-aware maximum parsimony. 

## Installation and dependencies
The code does not require installation. Download the repository and execute the indelMaP_MSA.py script for multiple sequence alignment and the indelMaP_ASR.py script for ancestral sequence reconstruction. The code is executed using python and uses the following libraries: ete3, NumPy, random, math, SciPy, os, and argparse.

### Using `uv`
`uv` is a python dependency management tool which allows for simpler distribution of python tools. You will need to
[install `uv`](https://docs.astral.sh/uv/getting-started/installation/) to begin. 

Thereafter, you can run:

```bash
uv tool install git+https://github.com/acg-team/indelMaP
```

This will install two scripts into your commandline: `indelmap-asr` and `indelmap-msa`. Just run either of them with the 
`--help` flag for more information. For the rest of this guide, you can substitute `python indelMaP_ASR.py` and 
`python indelMaP_MSA.py`with `indelmap-asr` and `indelmap-msa` respectively. 

## Usage

### Minimal usage
For minimal usage, indelMaP_MSA.py requires a sequences file in fasta format and a corresponding guide tree in newick format, as well as the required alphabet (protein or DNA). Minimal usage for indelMaP_ASR.py is similar. It requires a multiple sequence alignment in fasta format, a corresponding guide tree in newick format, and the required alphabet (protein or DNA). To specify the path for the output files and name, use the option --output_file without a suffix.

```
python indelMaP_MSA.py --seq_file sequenceFile.fasta --tree_file treeFile.nwk --alphabet Protein --output_file outputFile
```

```
python indelMaP_ASR.py --msa_file multipleSequenceAlignmentFile.fasta --tree_file treeFile.nwk --alphabet Protein  --output_file outputFile 
```
### Substitution model
The default substitution model for DNA sequences is K80 with a transition to transversion ratio of 2, and for protein sequences, the WAG model. With the option --RateMatrix a substitution model can be specified. For DNA sequences K80 and JC69 are supported, and for protein sequences, WAG and blosum. 
```
python indelMaP_MSA.py --seq_file sequenceFile.fasta --tree_file treeFile.nwk --alphabet DNA --RateMatrix K80{3,1}
```
or
```
python indelMaP_MSA.py --seq_file sequenceFile.fasta --tree_file treeFile.nwk --alphabet Protein --RateMatrix blosum
```
If each substitution should be associated with the same cost, no substitution matrix should be specified, and the option should be set to None.
```
python indelMaP_MSA.py --seq_file sequenceFile.fasta --tree_file treeFile.nwk --alphabet Protein --RateMatrix None
```
There is also the option to specify a symmetric rate matrix with an average substitution rate of 1. Columns must be separated by a comma and rows with a colon, and transition rates are given as float. 
```
python indelMaP_MSA.py --seq_file sequenceFile.fasta --tree_file treeFile.nwk --alphabet Protein --RateMatrix=-1,0.333,0.333,0.333:0.333,-1,0.333,0.333:0.333,0.333,-1,0.333:0.333,0.333,0.333,-1
```

### Gap costs
Gap costs are calculated by scaling the average substitution rate by a gap opening factor and a gap extension factor. The default is affine gap costs with the factor 2.5 for gap opening and 0.5 for gap extension. 
For example, linear gap costs are specified by setting both factors to the same value.
```
python indelMaP_ASR.py --msa_file multipleSequenceAlignmentFile.fasta --tree_file treeFile.nwk --alphabet Protein  --gap_opening_factor 1 --gap_extension_factor 1 
```

### Branch lengths
Parsimony methods usually do not consider branch lengths in their inference. Since we base the cost for each evolutionary event on a substitution model, we can consider branch lengths. The cost matrix is not recalculated for each branch length to save computational time. Instead, either one cost matrix for an average distance is used (--branch_length False), or four cost matrices are calculated based on pre-defined distances (--branch_length True). The cost matrix derived from the distance closest to the length of the branch is chosen during inference. The latter is the default behaviour. 
```
python indelMaP_ASR.py --msa_file multipleSequenceAlignmentFile.fasta --tree_file treeFile.nwk --alphabet Protein  --branch_length False
```

### Indel-awareness
If required, the option --indel_aware can be set to False. This will stop the algorithm from distinguishing between insertions and deletions. 
```
python indelMaP_MSA.py --seq_file multipleSequenceAlignmentFile.fasta --tree_file treeFile.nwk --alphabet Protein  --indel_aware False
```

### Indel-aware parsimony score
If only the indel-aware parsimony score is required, the option --ancestral_reconstruction for indelMaP_ASR can be set to False. The algorithm will only calculate the indel-aware parsimony score for a given alignment and guide tree. 
```
python indelMaP_ASR.py --msa_file multipleSequenceAlignmentFile.fasta --tree_file treeFile.nwk --alphabet Protein  --ancestral_reconstruction False
```
