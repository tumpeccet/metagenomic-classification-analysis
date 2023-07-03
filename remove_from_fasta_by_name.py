from Bio import SeqIO
from Bio.Seq import Seq
import random
import argparse

parser = argparse.ArgumentParser(
    prog='Removing sequences by matching string',
    description='Modifies the input .fasta file so that the bacterias containing the given string in the name (Bio.SeqRecord.description) are removed')
parser.add_argument("input", help="Input file path containing reads, must be in fasta format")
parser.add_argument("output", help="Output file path, can be an inexistent file in which case it is created, must be in fasta format")
parser.add_argument("string", help="The string for matching the sequence names to be deleted")
args = parser.parse_args()

matching_string = args.string
with open(args.output, "w") as output_file, open(args.input, "r") as original:
    records = SeqIO.parse(original, "fasta")
    for record in records:
        if matching_string in record.description:
            print(record.description)
            continue
        SeqIO.write(record, output_file, 'fasta')


