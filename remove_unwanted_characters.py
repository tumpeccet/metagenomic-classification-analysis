from Bio import SeqIO
from Bio.Seq import Seq
import random
import argparse

parser = argparse.ArgumentParser(description='Removes sequences with unwanted characters')
parser.add_argument("input", help="Input file path containing reads, must be in fasta format")
parser.add_argument("output", help="Output file path, can be an inexistent file in which case it is created, must be in fasta format")
args = parser.parse_args()

allowed = list("ACTG")
with open(args.output, "w") as output_file, open(args.input, "r") as original:
	records = SeqIO.parse(original, "fasta")
	for record in records:
		sequence = str(record.seq).upper()
		sequence = list(sequence)
		izbaci = False
		for i in range (len(sequence)):
			if sequence[i] not in allowed:
				izbaci = True
				break
		if not izbaci:
			sequence = "".join(sequence)
			record.seq = Seq(sequence)
			print("%s %s" % (record.id, record.seq))
			SeqIO.write(record, output_file, 'fasta')
		id_B_S, E = record.description.split("end=")
		record.description = id_B_S + "length=" + str(len(record.seq))
		sequences.append(record)

with open(args.output, "w") as output_handle:
    SeqIO.write(sequences, output_handle, "fasta")
