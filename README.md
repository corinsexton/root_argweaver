# ARGweaver pipeline for ROOT
### root_argweaver

This repository contains the scripts needed to generate a run of ARGweaver based on a single fasta alignment file. After 
ARGweaver is run, the pipeline will parse the smc file from the 1000th iteration to generate a minimum pairwise distance 
matrix across all trees in that smc file.

## Dependencies

This pipeline depends on the following:
- ARGweaver (https://github.com/mdrasmus/argweaver)
- R (specifically the `ape` package)
- Python 2.7

## Running the pipeline

To run the pipeline, run the following command from the top directory of this git repository:
     ```./pipeline.sh <aligned.fa> <label> [argweaver paramaters]```
Any additional argweaver parameters provided will be passed on to the `arg-sample` command. To see options for this command 
please run `arg-sample -h`.

The `pipeline.sh` script has comments throughout explaining what it does, but to summarize:
1. Creates an output label to be used for all generated files. This is in the format of the name of the provided fasta file followed by the label provided on the command line, separated with an underscore. For example, if the fasta provided was `HLA-B.fa` and the label was `recombrate_100xslower`, the output label would be `HLA-B_recombrate_100xslower`
2. Creates a new formatted fasta file (named as specified 1) based on the alignment file provided:
  - Replaces '-' with 'N'
  - Replaces '*' in labels with '_'
3. Makes a new directory (named as specified 1) and changes into that directory
4. Converts the newly formatted fasta file to a sites file using `fasta2sites.py`
5. Runs ARGweaver with default parameters unless any parameters are added on the command line. For example, a max time of 2000e3 and a recombination rate of 0.015e-8 would be specified with the following command:
     ```./pipeline.sh HLA-B.fa maxtime_2000e3_100xslower --maxtime=2000e3 --recombrate=.015e-8```
6. R script `pairwise_parser.R` runs to parse the trees from the 1000th iteration smc file. Across all trees in the file, the minimum distance between two nodes is calculated and output in a matrix. which is saved as `output_label.min_dist_matrix.tsv`.
7. Makes a new directory to store the smc files (`smc_files_output_label`) and moves all smc files into that directory.
8. FINISHED.

## Additional Scripts

A few other utility scripts are included that are not directly used in the `pipeline.sh` script.

##### `pull_trees.py`
This script is used to generate create separate Newick tree files for each tree in a single smc file. For example:
     ```python pull_trees.py HLA-B.1000.smc.gz```
The output will be as a file for each tree that exists in the smc file. The files will be found in the same directory as the 
specified smc file. They are named by the region they are found in, for example: `HLA-B.1000.323_580.newick`. To aggregate all 
the files, it is recommend to put them in one directory like so:
     ```mkdir HLA-B_trees; mv *.newick HLA-B_trees```

##### `fasta2sites.py`
This script converts aligned fasta files to `.sites` format (see https://github.com/mdrasmus/argweaver). Usage:
     ```python fasta2sites.py <input fasta file> <output sites file>```
   
##### `pairwise_parser.R`
This script calculates the minimum distance matrix across all trees in an smc file. Usage:
     ```Rscript pairwise_parser.R <input smc file> <output tsv matrix>```
