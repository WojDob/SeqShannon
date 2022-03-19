#!/usr/bin/python3

import argparse
import math
import os
import sys

from Bio import SeqIO


EXIT_SUCCESS = 0
EXIT_FAILURE = 1


def main():
    parser = argparse.ArgumentParser(
        description="Calculate shannon entropy for biological sequences."
    )
    parser.add_argument("-i", "--input", help="Input fasta file ")
    parser.add_argument("-o", "--output", help="Output file")
    args = parser.parse_args()

    if not args.input:
        parser.print_usage()
        return sys.exit(EXIT_FAILURE)

    results = parse_input_file(args.input)

    if args.output:
        write_dict(results, args.output)
    else:
        print_dict(results)

    return sys.exit(EXIT_SUCCESS)


def calculate_shannon_entropy(sequence):
    seq_list = list(sequence)
    unique_symbols = set(seq_list)
    M = float(len(seq_list))
    entropy_list = []
    for x in unique_symbols:
        n_i = seq_list.count(x)
        P_i = n_i / M
        entropy_i = P_i * (math.log(P_i, 2))
        entropy_list.append(entropy_i)

    sh_entropy = -(sum(entropy_list))

    return sh_entropy


def parse_input_file(input_file):
    records = (r for r in SeqIO.parse(input_file, "fasta"))
    results = {}
    for r in records:
        results[r.id] = calculate_shannon_entropy(r.seq)
    return results


def print_dict(d):
    for k, v in d.items():
        print(k, v)


def write_dict(d, file_path):
    with open(file_path, "w") as f:
        for k, v in d.items():
            f.write(f"{k} {v}\n")


if __name__ == "__main__":
    main()
