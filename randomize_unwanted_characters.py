from Bio import SeqIO
from Bio.Seq import Seq
import random
import argparse

parser = argparse.ArgumentParser(description='Randomizes ambiguous characters in genomic sequences.')
parser.add_argument("input", help="Input file path containing reads, must be in fasta format")
parser.add_argument("output", help="Output file path, can be an inexistent file in which case it is created, must be in fasta format")
args = parser.parse_args()

allowed = list("ACTG")
with open(args.output, "w") as output_file, open(args.input, "r") as original:
	records = SeqIO.parse(original, "fasta")
	for record in records:
		sequence = str(record.seq).upper()
		sequence = list(sequence)
		for i in range (len(sequence)):
			if sequence[i] not in allowed:
				rand_str = random.choice(list('ATCG'))
				print("replacing %s with %s", sequence[i], rand_str) 
				sequence[i] = rand_str
		sequence = "".join(sequence)
		record.seq = Seq(sequence)
		SeqIO.write(record, output_file, 'fasta')
