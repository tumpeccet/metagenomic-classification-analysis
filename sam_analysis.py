import re
import argparse

parser = argparse.ArgumentParser(
    description='Analyses SAM file')
parser.add_argument("input", help="input file path")
parser.add_argument("header", help="header file path")
parser.add_argument("output", help="Output file path for abundances")
args = parser.parse_args()

with open(args.input, "r") as input_file, open(args.output, "w") as output_file, open(args.header, "r") as headers:
    # creating a header reference list
    lines = headers.readlines()
    tax_id_to_name = {}
    for header in lines:
        if header.startswith('>'):
            tax_id_extended, long_name = header.strip().split(" ", 1)
            _, tax, _ = tax_id_extended.strip().split('|')
            tax_id_to_name[tax] = long_name

    lines = input_file.readlines()
    results = {}
    max_values = {}
    transformed_rows = {}
    for line in lines:
        parts = re.split(r'\t+', line.strip())
        if parts[1].strip() != "0" and parts[1].strip() != "16":
            continue
        read_id = parts[0].strip()
        tax_id_extended = parts[2].strip()

        parts2 = re.split(r'\|+', tax_id_extended.strip())
        tax_id = parts2[1].strip()
        name = tax_id_to_name[tax_id]

        value_cig = 0.0
        if len(parts) >= 14:
            nm = parts[13].strip()
            parts3 = re.split(r':+', nm.strip())
            value_cig = int(parts3[-1].strip())

        if read_id in results:
            if value_cig > max_values[read_id]:
                results[read_id] = []
                results[read_id].append(tax_id)
                max_values[read_id] = value_cig
            elif value_cig == max_values[read_id]:
                results[read_id].append(tax_id)
        else:
            results[read_id] = []
            results[read_id].append(tax_id)
            max_values[read_id] = value_cig

    tax_id_classified = {}
    for read_id in results:
        tax_ids = results[read_id]
        values = {}
        for tax_id in tax_ids:
            if tax_id not in values:
                values[tax_id] = 0
            values[tax_id] += 1
        max_tax_id = ""
        max_value = 0
        for tax_id in values:
            value = values[tax_id]
            if value > max_value:
                max_value = value
                max_tax_id = tax_id
        transformed_rows[read_id.strip()] = (max_tax_id, tax_id_to_name[max_tax_id])
        if max_tax_id not in tax_id_classified:
            tax_id_classified[tax_id] = 0
        tax_id_classified[max_tax_id] += 1

    output_file.write("METAFLYE CONTIG MINIMAP ALIGNMENT\n")
    output_file.write("CONTIG NAME\tTAX ID\tSEQUENCE\n")
    for read_id in transformed_rows:
        tax, seq_name = transformed_rows[read_id]
        output_file.write(read_id + "\t" + tax + "\t" + seq_name + "\n")
    output_file.write("\n\n\nMINIMAP ALIGNMENT SEQUENCE READ COUNT REPORT\n")
    output_file.write("tax id\tname\tread_count\n")
    for tax_id in tax_id_classified:
        output_file.write(tax_id + "\t" + tax_id_to_name[tax_id] + "\t" + str(tax_id_classified[tax_id]) + "\n")

