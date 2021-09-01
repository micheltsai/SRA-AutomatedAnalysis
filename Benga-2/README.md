# Benga
Bacterial Epidemiology NGs Analysis (BENGA) framework and pipeline.

# Requirements
* python 3.6+
  * biopython
  * scipy
  * numba
  * numpy
  * pandas
  * fastcluster
  * matplotlib
* blast
* prodigal
* prokka
* roary

# Usage
***cgMLST profiling***
```
profiling -i fasta_file -o result.tsv --scheme scheme.faa --prodigaltf train_file.trn
```
***cgMLST makedb***
```
makedb.py -i genomes -o output_path
```
