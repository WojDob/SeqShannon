# SeqShannon

SeqShannon is a Python package that calculates the Shannon entropy of biological sequences. It works for both nucleotide and amino acid sequences.

## Description

SeqShannon is a Python package that reads a FASTA file and calculates the Shannon entropy for each biological sequence in the file. It utilizes the BioPython library for parsing FASTA files. The results can be easily written to an output file or printed in the console.

The Shannon entropy is a measure of the uncertainty or randomness of a set of data. In the context of biological sequences, such as DNA or protein sequences, the Shannon entropy can provide insights into the variability and complexity of the sequence.

## Installation

You can install SeqShannon using `pip`:

```bash
pip install seqshannon
```

## Usage
SeqShannon can be used as a command-line tool or as a Python library.

### Command-line usage

Use a fasta file as input. By default, the package prints out the identifier and calculated Shannon entropy for each sequence in the input file.

```bash
seqshannon -i example.fasta
```

You can also specify a file to save the output.

```bash
seqshannon -i example.fasta -o output.txt
```

#### Example

Here are the contents of an example FASTA file:

```fasta
>example_1
VLSISYSRSESSLE
>example_2
TIGQRKPSTFSWSS
>example_3
RAASRSSWERGP
```

Running SeqShannon on this file will yield the following output:

```bash
example_1 2.4137995646056805
example_2 3.09306920777189
example_3 2.6258145836939115
```


### Python library usage
You can calculate the Shannon entropy of the given sequence by importing `shannon_entropy`.

```python
>>> from seqshannon import shannon_entropy
>>> from Bio.Seq import Seq
>>> sequence = Seq("ATGCATGC")
>>> entropy = shannon_entropy(sequence)
>>> print(entropy)
2.0
```

## Contact

For any issues or suggestions, please contact [Wojciech Dobrychłop](mailto:wojciech.dobrychlop@gmail.com).

## License

SeqShannon is licensed under the MIT License. See the [LICENSE](LICENSE) file for more details.
