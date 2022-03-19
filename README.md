# Sequence-Shannon
Calculate Shannon entropy for biological sequences. Works for both nucleotide and amino acid sequences.

## Installation

Requires [pip](https://choosealicense.com/licenses/mit/) and [virtualenv](https://virtualenv.pypa.io/en/latest/).

```bash
git clone git@github.com:WojDob/sequence-shannon.git
cd sequence-shannon
virtualenv venv
source venv/bin/activate
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
