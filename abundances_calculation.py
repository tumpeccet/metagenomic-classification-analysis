import argparse
import re
from math import ceil
from statistics import mean
from Bio import SeqIO

parser = argparse.ArgumentParser(
    prog='Calculating abundances.',
    description='')
parser.add_argument("sample", help="Input sample file path, must be in fastq format")
parser.add_argument("database_headers", help="Input file path containing all headers from the reads in the database, must be in fasta format")
parser.add_argument("sample_headers", help="Input file path containing all headers from the reads in the sample, must be in fasta format")
parser.add_argument("assembly_info", help="Input file path containing assembly info from metaflye contigs, must be in txt format")
parser.add_argument("minimap_contigs", help="Input file path containing minimap output for assembled contigs, must be in SAM format")
parser.add_argument("sample_stats", help="Input file containing badread sample stats generated using \"seqkit stats\"")
parser.add_argument("minimap_reads", help="Input file path containing minimap output for assembled badread reads, must be in SAM format")
parser.add_argument("kraken_output_reads", help="Input file path containing kraken output for reads classification")
parser.add_argument("kraken_output_contigs", help="Input file path containing kraken output for contigs classification")
parser.add_argument("output_file", help="Output file path for the report, should be .txt or .csv")
args = parser.parse_args()

#parse minimap output
def parse_sam_minimap(sam_file):
    results = {}
    max_cigar_values = {}
    transformed_rows = {}

    sam_lines = sam_file.readlines()
    for line in sam_lines:
        parts = re.split(r'\t+', line.strip())
        bitwise_flag = parts[1].strip()
        # the reads with flag 0 are mapped to the forward strand and the reads with flag 16 are mapped to the reverse strand
        if bitwise_flag != "0" and bitwise_flag != "16":
            continue

        read_id = parts[0].strip()
        tax_id_extended = parts[2].strip()

        parts2 = re.split(r'\|+', tax_id_extended.strip())
        tax_id = parts2[1].strip()

        value_cig = 0.0
        if len(parts) >= 14:
            nm = parts[13].strip()
            parts3 = re.split(r':+', nm.strip())
            value_cig = int(parts3[-1].strip())

        if read_id in results:
            if value_cig > max_cigar_values[read_id]:
                results[read_id] = []
                results[read_id].append(tax_id)
                max_cigar_values[read_id] = value_cig
            elif value_cig == max_cigar_values[read_id]:
                results[read_id].append(tax_id)
        else:
            results[read_id] = []
            results[read_id].append(tax_id)
            max_cigar_values[read_id] = value_cig

    tax_id_classified = {}
    tax_to_contig = {}
    for read_id in results:
        # counting alignments for tax IDs
        tax_ids = results[read_id]
        aligned = {}
        for tax_id in tax_ids:
            if tax_id not in aligned:
                aligned[tax_id] = 0
            aligned[tax_id] += 1

        # finding the resulting tax ID
        max_tax_id = ""
        max_value = 0
        for tax_id in aligned:
            value = aligned[tax_id]
            if value > max_value:
                max_value = value
                max_tax_id = tax_id
        transformed_rows[read_id.strip()] = max_tax_id
        if max_tax_id not in tax_to_contig:
            tax_to_contig[max_tax_id] = []
        tax_to_contig[max_tax_id].append(read_id.strip())

        # counting the chosen aligned tax IDs
        if max_tax_id not in tax_id_classified:
            tax_id_classified[max_tax_id] = 0
        tax_id_classified[max_tax_id] += 1

    return tax_id_classified, transformed_rows, tax_to_contig

#parse sample stats
def parse_sample_stats(sample_stats_file):
    sample_lines = sample_stats_file.readlines()
    for line in sample_lines:
        if line.startswith("file"):
            continue
        # file format type num_seqs sum_len min_len avg_len max_len
        file, format, type, num_seqs, sum_len, min_len, avg_len, max_len = line.strip().split()

        num_seqs_parts = num_seqs.split(',')
        num_seqs_float = ""
        for part in num_seqs_parts:
            num_seqs_float += part
        num_seqs = float(num_seqs_float)

        avg_len_parts = avg_len.split(',')
        avg_len_float = ""
        for part in avg_len_parts:
            avg_len_float += part
        avg_len = float(avg_len_float)

        return num_seqs, avg_len

#parse kraken output reads
def parse_kraken_output(kraken_output_file):
    lines = kraken_output_file.readlines()
    tax_id_read_count = {}
    tax_to_reads = {}
    for line in lines:
        parts = re.split(r'\t+', line.strip())
        if parts[0].strip() == 'C':
            read_id = parts[1].strip()
            tax_id = parts[2].strip()
            if tax_id not in tax_id_read_count:
                tax_id_read_count[tax_id] = 0
                tax_to_reads[tax_id] = []
            tax_id_read_count[tax_id] += 1
            tax_to_reads[tax_id].append(read_id)
    return tax_id_read_count, tax_to_reads

#parse assembly info
def parse_assembly_info(assembly_file):
    assembly_lines = assembly_file.readlines()
    contig_len_cov = {}
    for line in assembly_lines:
        if line.startswith('#'):
            continue
        # #seq_name length cov. circ. repeat mult. alt_group graph_path
        name, length, coverage, _, _, _, _, _ = line.strip().split("\t")
        if name == "#seq_name":
            continue
        contig_len_cov[name] = (float(length), float(coverage))
    return contig_len_cov

def create_header_reference(header_file):
    header_lines = header_file.readlines()
    tax_id_to_name = {}
    for header in header_lines:
        if header.startswith('>'):
            tax_id_extended, long_name = header.strip().split(" ", 1)
            name = long_name.split(",")[0]
            _, tax, _ = tax_id_extended.strip().split('|')
            if tax not in tax_id_to_name:
                tax_id_to_name[tax] = []
            tax_id_to_name[tax].append(name)
    return tax_id_to_name

def calculate_contigs_abundances(tax_id_to_contig, contig_to_len_cov, avg_read_length, read_count):
    tax_id_to_abundance = {}

    for tax_id in tax_id_to_contig:
        num_reads = 0.
        for contig_name in tax_id_to_contig[tax_id]:
            length, coverage = contig_to_len_cov[contig_name]
            num_reads += ceil((length * coverage) / avg_read_length)
        tax_id_to_abundance[tax_id] = float(num_reads) / float(read_count)
    return tax_id_to_abundance

def calculate_reads_abundances(tax_id_classified, read_count):
    tax_id_to_abundance = {}
    for tax_id in tax_id_classified:
        tax_id_to_abundance[tax_id] = float(tax_id_classified[tax_id]) / float(read_count)
    return tax_id_to_abundance

def calculate_true_abundances(sample_file):
    records = SeqIO.parse(sample_file, "fastq")
    read_counts = {}
    lengths = []
    for record in records:
        seq_description = record.description.strip().split(' ')[1]
        lengths.append(len(record.seq))
        tax_id = seq_description.split('|')
        if len(tax_id) > 1:
            tax_id = tax_id[1]
            if tax_id not in read_counts:
                read_counts[tax_id] = 0
            read_counts[tax_id] += 1

    true_abundances = {}
    for tax_id in read_counts:
        true_abundances[tax_id] = float(read_counts[tax_id]) / float(len(lengths))
    return true_abundances, len(lengths), round(mean(lengths), 2)

with open(args.database_headers, "r") as database_headers, open(args.sample_headers, "r") as sample_headers, open(
        args.assembly_info, "r") as assembly_info, open(args.minimap_contigs, "r") as minimap_contigs_file, open(
    args.output_file, "w") as report_file, open(args.sample_stats, "r") as sample_stats, open(args.minimap_reads,
                                                                                              "r") as minimap_reads_file, open(
    args.kraken_output_reads, "r") as kraken_reads_file, open(args.kraken_output_contigs,
                                                              "r") as kraken_contigs_file, open(args.sample, "r") as sample_file:

    read_count, avg_badread_read_length = parse_sample_stats(sample_stats)
    true_abundances, records_count, avg_length = calculate_true_abundances(sample_file=sample_file)
    if records_count != read_count:
        print("Read count from sample file = " + str(records_count) + ", read count from stats file = " + str(read_count))
        read_count = records_count

    if avg_length != avg_badread_read_length:
        print("Avg length from sample file = " + str(avg_length) + ", avg length from stats file = " + str(avg_badread_read_length))
        avg_badread_read_length = avg_length

    # creating a header reference dictionary (tax ID to name)
    tax_id_to_name_database = create_header_reference(header_file=database_headers)

    # creating a sample reference dictionary (tax ID to name)
    tax_id_to_name_sample = create_header_reference(header_file=sample_headers)

    # creating a contig name to length and coverage dictionary
    contig_len_cov = parse_assembly_info(assembly_file=assembly_info)

    # kraken classification on badread reads
    kraken_tax_to_read_count, kraken_tax_to_read = parse_kraken_output(kraken_output_file=kraken_reads_file)
    kraken_read_abundances = calculate_reads_abundances(kraken_tax_to_read_count, read_count)

    # kraken classification on contigs
    kraken_tax_to_contig_count, kraken_tax_to_contig = parse_kraken_output(kraken_output_file=kraken_contigs_file)
    kraken_contig_abundances = calculate_contigs_abundances(kraken_tax_to_contig, contig_len_cov, avg_badread_read_length, read_count)

    # creating the counter dictionary for classified sequences and contig to classified tax ID dictionary
    minimap_contig_tax_id_classified, minimap_contig_to_tax_id, minimap_tax_id_to_contig = parse_sam_minimap(sam_file=minimap_contigs_file)
    minimap_contig_abundances = calculate_contigs_abundances(minimap_tax_id_to_contig, contig_len_cov, avg_badread_read_length, read_count)

    # creating the counter dictionary for classified sequences and read to classified tax ID dictionary
    minimap_read_tax_id_classified, minimap_read_to_tax_id, minimap_tax_id_to_read = parse_sam_minimap(sam_file=minimap_reads_file)
    minimap_read_abundances = calculate_reads_abundances(minimap_read_tax_id_classified, read_count)

    report_file.write(
        "\nID\ttax ID\tname\tTRA\tRReadCA\tRContigA\tMARContigA\tMRReadCA\n")
    id = 1
    for tax_id in tax_id_to_name_sample:
        kraken_read_abundance = 0.
        if tax_id in kraken_read_abundances:
            kraken_read_abundance = kraken_read_abundances[tax_id]

        kraken_contig_abundance = 0.
        if tax_id in kraken_contig_abundances:
            kraken_contig_abundance = kraken_contig_abundances[tax_id]

        minimap_contig_abundance = 0.
        if tax_id in minimap_contig_abundances:
            minimap_contig_abundance = minimap_contig_abundances[tax_id]

        minimap_read_abundance = 0.
        if tax_id in minimap_read_abundances:
            minimap_read_abundance = minimap_read_abundances[tax_id]

        report_file.write(str(id) + "\t" + str(tax_id) + "\t" +
                          tax_id_to_name_sample[tax_id][0] + "\t" + str(round(true_abundances[tax_id] * 100, 2)) + "\t" + str(
            round(kraken_read_abundance * 100, 2)) + "\t" + str(round(kraken_contig_abundance * 100, 2)) + "\t" + str(
            round(minimap_contig_abundance * 100, 2)) + "\t" + str(round(minimap_read_abundance * 100, 2)) + "\n")
        id += 1