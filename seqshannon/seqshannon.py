#!/usr/bin/python3

import argparse
import math
import os
import sys

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

EXIT_SUCCESS = 0
EXIT_FAILURE = 1


def parse_args():
    """
    Parse command line arguments.

    Returns:
    Namespace: The parsed command line arguments.
    argparse.ArgumentParser: The argument parser.
    """
    parser = argparse.ArgumentParser(
        description="Calculate Shannon entropy of biological sequences."
    )
    parser.add_argument("-i", "--input", required=True, help="Input fasta file")
    parser.add_argument("-o", "--output", help="Output file")
    return parser.parse_args(), parser


def validate_input_file(input_file, parser):
    """
    Validate the input file.

    Parameters:
    input_file (str): Path to the input file.
    parser (argparse.ArgumentParser): The argument parser.

    Returns:
    None
    """
    if not os.path.isfile(input_file):
        print("Error: Input file does not exist.")
        parser.print_usage()
        sys.exit(EXIT_FAILURE)


def shannon_entropy(sequence):
    """
    Calculate the Shannon entropy of a sequence.

    Parameters:
    sequence (str or SeqRecord): The sequence to calculate the entropy of.

    Returns:
    float: The Shannon entropy of the sequence.
    """
    if isinstance(sequence, SeqRecord):
        sequence = str(sequence.seq)
    
    seq_list = list(sequence)
    unique_symbols = set(seq_list)
    M = float(len(seq_list))
    entropy_list = []
    for x in unique_symbols:
        n_i = seq_list.count(x)
        P_i = n_i / M
        entropy_i = P_i * (math.log(P_i, 2))
        entropy_list.append(entropy_i)

    sh_entropy = abs(sum(entropy_list))

    return sh_entropy


def parse_input_file(input_file):
    """
    Parse a FASTA file and calculate the Shannon entropy of each sequence.

    Parameters:
    input_file (str): Path to the input FASTA file.

    Returns:
    dict: A dictionary mapping sequence identifiers to their Shannon entropy.
    """
    records = (r for r in SeqIO.parse(input_file, "fasta"))
    return {r.id: shannon_entropy(r.seq) for r in records}


def print_dict(d):
    """
    Print a dictionary to the console.

    Parameters:
    d (dict): The dictionary to print.
    """
    for k, v in d.items():
        print(k, v)


def write_dict(d, file_path):
    """
    Write a dictionary to a file.

    Parameters:
    d (dict): The dictionary to write.
    file_path (str): The path to the output file.
    """
    try:
        with open(file_path, "w") as f:
            for k, v in d.items():
                f.write(f"{k} {v}\n")
    except IOError:
        raise IOError("Error: Cannot create or write to output file.")


def main():
    """
    Main function that parses command line arguments and calls the appropriate functions.
    """
    args, parser = parse_args()
    validate_input_file(args.input, parser)

    results = parse_input_file(args.input)

    if args.output:
        try:
            write_dict(results, args.output)
        except IOError as e:
            print(e)
            sys.exit(EXIT_FAILURE)
    else:
        print_dict(results)

    sys.exit(EXIT_SUCCESS)


if __name__ == "__main__":
    main()
