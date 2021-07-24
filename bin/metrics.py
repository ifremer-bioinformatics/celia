#!/usr/bin/env python3

import argparse
import os.path
from Bio import SeqIO


def get_args():
    parser = argparse.ArgumentParser(description="")
    # parser.add_argument('-f', dest="fasta", type=str, nargs='+', required=True, help='Input fasta')
    # parser.add_argument('-b', dest="busco", type=str, nargs='+', required=False, help='BUSCO short summary')
    parser.add_argument('-d', dest="dir_path", type=str, required=True, help='A directory with all results')

    args = parser.parse_args()

    return args


def Nx(numlist, reference_length, percentage=50.0):
    # Check percentage range
    assert percentage >= 0.0
    assert percentage <= 100.0
    # Always reverse sort
    numlist.sort(reverse=True)
    # Duplicate value
    s = reference_length
    # Calculate limit
    limit = reference_length * (100.0 - percentage) / 100.0
    # Compute Nx value
    for l in numlist:
        s -= l
        if s <= limit:
            nx = l
            return nx


def base_count(seq, base):
    cnt = seq.seq.count(str(base))
    return cnt


def get_length(fasta):
    f = open(fasta, "r")
    for s in SeqIO.parse(f, "fasta"):
        print(s.id + '\t' + str(len(s.seq)))


def get_stats(fasta, stats):
    f = open(fasta, "r")
    # Data structure
    s = {'seq': 0,
         'residue': 0,
         'A': 0,
         'T': 0,
         'G': 0,
         'C': 0,
         'N': 0,
         'GC': 0,
         'seqs_length': [],
         'min': 0,
         'max': 0,
         'aver': 0,
         'N50': 0,
         'N80': 0,
         'N90': 0}
    # Iterate over each sequence
    filename = os.path.basename(fasta).replace('.clean.fasta', '')
    for seq in SeqIO.parse(f, "fasta"):
        s['seq'] += 1
        s['residue'] += len(seq.seq)
        s['seqs_length'].append(len(seq.seq))
        s['A'] += seq.seq.count(str("A"))
        s['T'] += seq.seq.count(str("T"))
        s['G'] += seq.seq.count(str("G"))
        s['C'] += seq.seq.count(str("C"))
        s['N'] += seq.seq.count(str("N"))
    # Compute global metrics
    s['max'] = max(s['seqs_length'])
    s['min'] = min(s['seqs_length'])
    s['aver'] = round(sum(s['seqs_length']) / len(s['seqs_length']), 2)
    s['GC'] = round((s['G'] + s['C']) / (s['G'] + s['C'] + s['A'] + s['T']) * 100, 2)
    # Compute Nx values
    s['N50'] = Nx(s['seqs_length'], s['residue'])
    s['N80'] = Nx(s['seqs_length'], s['residue'], percentage=80.0)
    s['N90'] = Nx(s['seqs_length'], s['residue'], percentage=90.0)
    # Store result
    stats[filename] = s
    # Return data
    return stats


def get_busco_score(busco, stats):
    b = open(busco, "r")
    # Data storage
    s = {'busco_C': 0,
         'busco_S': 0,
         'busco_D': 0,
         'busco_F': 0,
         'busco_M': 0,
         'busco_T': 0,
         'busco_C_perc': 0,
         'busco_S_perc': 0,
         'busco_D_perc': 0,
         'busco_F_perc': 0,
         'busco_M_perc': 0}
    columns = ['Complete (C) and single-copy (S)', 'Complete (C) and duplicated (D)', 'Fragmented (F)', 'Missing (M)']
    # Read each file and extract BUSCO scores
    specie = ".".join(busco.split("/")[-1].split(".")[3:-1])
    # Extract raw count
    for line in b:

        if "Complete BUSCOs" in line:
            s['busco_C'] = int(line.split("\t")[1])
        elif "Complete and single-copy BUSCOs" in line:
            s['busco_S'] = int(line.split("\t")[1])
        elif "Complete and duplicated BUSCOs" in line:
            s['busco_D'] = int(line.split("\t")[1])
        elif "Fragmented BUSCOs" in line:
            s['busco_F'] = int(line.split("\t")[1])
        elif "Missing BUSCOs" in line:
            s['busco_M'] = int(line.split("\t")[1])
    s['busco_T'] = s['busco_C'] + s['busco_F'] + s['busco_M']
    # Make percentage from raw count
    s['busco_C_perc'] = round(s['busco_C'] / float(s['busco_T']) * 100, 1)
    s['busco_S_perc'] = round(s['busco_S'] / float(s['busco_T']) * 100, 1)
    s['busco_D_perc'] = round(s['busco_D'] / float(s['busco_T']) * 100, 1)
    s['busco_F_perc'] = round(s['busco_F'] / float(s['busco_T']) * 100, 1)
    s['busco_M_perc'] = round(100 - s['busco_C_perc'] - s['busco_F_perc'], 1)
    # Add results to main dict - assume specie already exist -> ok with nextflow pipeline
    for k in s.keys():
        stats[specie][k] = s[k]
    # Return results
    return stats


def printing(stats):
    i = 0
    for k, v in stats.items():
        header = ['SeqName']
        data = [k]
        for kv, vv in v.items():
            if not kv == 'seqs_length':
                header.append(kv)
                data.append(vv)
        if i == 0:
            print(";".join(header))
            i += 1
        print(";".join(map(str, data)))


def main(args):
    # Init storage
    stats = {}
    fasta = []
    busco = []
    bowtie2 = []
    # Step 1 - Collect files and search for each type
    listdir = (os.listdir(args.dir_path))
    for f in listdir:
        if '.fasta' in f:
            fasta.append(f)
        elif '.log' in f:
            bowtie2.append(f)
        elif 'short_summary' in f:
            busco.append(f)
    # Step 2 - Collect fasta metrics
    for fa in fasta:
        stats = get_stats(fa, stats)
    # Step 3 - Collect BUSCO scores
    for busco in busco:
        stats = get_busco_score(busco, stats)
    # Step 4 - Print results
    printing(stats)


if __name__ == '__main__':
    main(get_args())
