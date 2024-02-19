# NanoTipSeq
An analysis pipeline for predicting insertions in nanopore long reads.
## Installation
clone the repository with 
```bash
git clone https://github.com/JHanLab/NanoTipSeq.git
```
Create a conda environment
```bash
conda env create -f nanopore.yml
```
## Usage
The inputs of this script are your bam alignment file and the polyA count for your search.
```bash
bash NanoTipSeq.sh <alignment.bam> <polyA count>
```
## Output Files
This program should produce two output files: a .bed file and a .txt file. The pA.merged.new_inst-3_sup.txt file details the list of predicted insertions.
