# Sequence-Shannon
Calculate Shannon entropy for biological sequences. Works for both nucleotide and amino acid sequences.

## Description
This Python script takes in a FASTA file and calculates the Shannon entropy for each biological sequence in the file. It utilizes the BioPython library for parsing FASTA files, and its results can be easily written to an output file or printed in the console.

The Shannon entropy is a measure of the uncertainty or randomness of a set of data. In the context of biological sequences, such as DNA or protein sequences, the Shannon entropy can provide insights into the variability and complexity of the sequence.


## Installation

Requires [pip](https://choosealicense.com/licenses/mit/). 

```bash
git clone git@github.com:WojDob/sequence-shannon.git
cd sequence-shannon
pip install -r requirements.txt
```

## Usage 

Use a fasta file as input. By default the script prints out the identifier and calculated Shannon entropy for each sequence in the input file.
```bash
$ python3 sequence-shannon.py -i example.fasta
example_1 2.4137995646056805
example_2 3.09306920777189
example_3 2.6258145836939115
```

You can also specify a file to save the output.
```bash
$ python3 sequence-shannon.py -i example.fasta -o output.txt
$ cat output.txt 
example_1 2.4137995646056805
example_2 3.09306920777189
example_3 2.6258145836939115
```
