from unittest.mock import patch

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
# if ImportError: pip install -e .
from seqshannon.seqshannon import (parse_input_file, print_dict,
                                   shannon_entropy, write_dict)


@pytest.mark.parametrize("sequence, expected", [
    (SeqRecord(Seq("ATGCATGC")), 2.0),
    (SeqRecord(Seq("VLSISYSRSESSLE")), 2.4137995646056805),
    (SeqRecord(Seq("TIGQRKPSTFSWSS")), 3.09306920777189),
    (SeqRecord(Seq("RAASRSSWERGP")), 2.6258145836939115),
    (("0000"), 0.0),
    (("oioi"), 1.0),
])
def test_shannon_entropy(sequence, expected):
    assert abs(shannon_entropy(sequence) - expected) < 1e-6

def test_parse_input_file():
    results = parse_input_file('example.fasta')
    assert len(results) == 3
    assert abs(results['example_1'] - 2.4137995646056805) < 1e-6
    assert abs(results['example_2'] - 3.09306920777189) < 1e-6
    assert abs(results['example_3'] - 2.6258145836939115) < 1e-6

@patch('builtins.print')
def test_print_dict(mock_print):
    d = {'example_1': 2.4137995646056805, 'example_2': 3.09306920777189, 'example_3': 2.6258145836939115}
    print_dict(d)
    mock_print.assert_any_call('example_1', 2.4137995646056805)
    mock_print.assert_any_call('example_2', 3.09306920777189)
    mock_print.assert_any_call('example_3', 2.6258145836939115)

def test_write_dict(tmpdir):
    d = {'example_1': 2.4137995646056805, 'example_2': 3.09306920777189, 'example_3': 2.6258145836939115}
    file_path = tmpdir.join("test.txt")
    write_dict(d, file_path)
    with open(file_path, "r") as f:
        lines = f.readlines()
    assert 'example_1 2.4137995646056805\n' in lines
    assert 'example_2 3.09306920777189\n' in lines
    assert 'example_3 2.6258145836939115\n' in lines
