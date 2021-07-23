#!/usr/bin/env python3

import os
import sys
import argparse

"""
Authors: 
        Developer: Emeric TEXERAUD
        Maintainer: Alexandre CORMIER
                    Bioinformatics engineer
                    SeBiMER, Ifremer
Creation Date: 2021
Copyright (c): SeBiMER, february-2020
Licence: MIT
"""


def validate_path(dir_path, dir_type):
    if not os.path.exists(dir_path):
        print(dir_type + ' does not exist')
        exit(1)


def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-i', dest="indir", type=str, required=True, help='Directory containing antiSMASH results')
    parser.add_argument('-o', dest="outdir", type=str, required=True, help='Output directory')

    arg = parser.parse_args()

    validate_path(arg.indir, "Input Directory")

    return arg


list_annotation = ['T1PKS', 'transAT-PKS', 'transAT-PKS-like', 'T2PKS', 'T3PKS', 'hglE-KS', 'PpyS-KS', 'PKS-like',
                   'arylpolyene', 'resorcinol', 'ladderane', 'PUFA', 'NRPS', 'NRPS-like', 'CDPS', 'thioamide-NRP',
                   'terpene', 'lanthipeptide', 'bacteriocin', 'RaS-RiPP', 'fungal-RiPP', 'thiopeptide', 'LAP',
                   'TfuA-related', 'linaridin', 'cyanobactin', 'glycocin', 'lassopeptide', 'sactipeptide',
                   'bottromycin', 'head_to_tail', 'microviridin', 'proteusin', 'lipolanthine', 'blactam', 'amglyccycl',
                   'aminocoumarin', 'siderophore', 'ectoine', 'NAGGN', 'butyrolactone', 'indole', 'nucleoside',
                   'phosphoglycolipid', 'melanin', 'oligosaccharide', 'furan', 'hserlactone', 'acyl_amino_acids',
                   'tropodithietic-acid', 'phenazine', 'phosphonate', 'fused', 'PBDE', 'other', 'betalactone',
                   'saccharide', 'fatty_acid', 'halogenated']

dict_increment = {}


def read_gbk_file(file):
    """
    Create a SeqIO object decoding a genbank file.
    :param file: str
    :return: SeqIO object
    """
    from Bio import SeqIO
    return SeqIO.read(open(file, "r"), "genbank")


def find_key_index(gbk_file):
    """
    Using the filename, find the key used to find informations in the dictionary from the index.html parsed.
    :param gbk_file: str
    :return: str
    """
    name = os.path.basename(gbk_file).split('.region')[0]
    num = int(os.path.basename(gbk_file).split('.region')[1].split('.gbk')[0])
    return str(name) + '-' + str(num)


def list_protocore_number(gbk_decoded):
    """
    Find the place of each protocore in the original file.
    :param gbk_decoded: SeqIO object
    :return: list(int)
    """
    list_protocore = list()
    for i in range(0, len(gbk_decoded.features)):
        if gbk_decoded.features[i].type == 'proto_core':
            list_protocore.append(i)
    print('Protocore not found' if len(list_protocore) == 0 else 'protocore found')
    return list_protocore


def find_id(gbk_decoded, filepath):
    """
    Split the ID to only get the part for the sequence and not contig.
    For the genomes assembled in the lab, take the name of the input subdir which contain the ID,
    as the ID is not correctly found by antiSMASH and thus not in the ID section of the SeqIO genbank object.
    :param gbk_decoded: SeqIO object
    :param filepath: str
    :return: str
    """
    try:
        identifier = str(os.path.basename(os.path.dirname(filepath)).split('_')[0]) if ('node' or 'contig') in \
                                                                                       gbk_decoded.id.lower() else \
            gbk_decoded.id.split('.')[0][0:gbk_decoded.id.find('_') + 7]
    except AttributeError:
        identifier = 'no ID found'
    return str(identifier)


def find_organism(gbk_decoded, identifier):
    """
    If the organism is not found, return 'organism ID' + the genome ID.
    :param gbk_decoded: SeqIO object
    :param identifier: str
    :return: str
    """
    try:
        organism = 'organism ID ' + identifier if gbk_decoded.annotations['organism'] == '' \
            else gbk_decoded.annotations['organism']
    except (KeyError, AttributeError):
        organism = 'organism ID ' + identifier
    return organism


def find_run_date(gbk_decoded):
    """
    :param gbk_decoded: SeqIO object
    :return: str (date format if found)
    """
    try:
        run_date = gbk_decoded.annotations['structured_comment']['antiSMASH-Data']['Run date'][0:19]
    except KeyError:
        run_date = 'no date found'
    return run_date


def find_contig(gbk_decoded):
    """
    The contig name and number can be split by spaces, so we take all words that contains 'contig' or 'node'.
    :param gbk_decoded: SeqIO object
    :return: str
    """
    contig = ''
    try:
        for i in gbk_decoded.description.split(' '):
            if 'contig' or 'node' in i.lower():
                contig += str(i)
    except AttributeError:
        contig = 'contig not found, see name of the antiSAMSH file to find the contig'
    return 'contig not found, see name of the antiSMASH file to find the contig' if contig == '' \
        else contig.replace(',', '')


def find_start(gbk_decoded):
    """
    This is the start of the DNA sequence of the cluster.
    :param gbk_decoded: SeqIO object
    :return: int if found else str
    """
    try:
        start = gbk_decoded.annotations['structured_comment']['antiSMASH-Data']['Orig. start']
    except (KeyError, AttributeError):
        start = 'start not found'
    return start


def find_end(gbk_decoded):
    """
    This is the end of the DNA sequence of the cluster.
    :param gbk_decoded: SeqIO object
    :return: int if found else str
    """
    try:
        end = gbk_decoded.annotations['structured_comment']['antiSMASH-Data']['Orig. end']
    except (KeyError, AttributeError):
        end = 'end not found'
    return end


def find_length(start, end):
    """
    :param start: int or str
    :param end: int or str
    :return: int if found else str
    """
    try:
        length = str(int(end) - int(start))
    except (TypeError, ValueError):
        length = 'could not calculate the length, please see the start and end'
    return length


def increment_dict_increment(identifier):
    """
    dict_increment is a global dictionnary containing the current number of protocores found for each genome.
    It is used to keep track of the number of protocore and this number is used in the name of the sequences,
    helping find them when there are more than one protocore in a contig.
    :param identifier: str
    """
    if identifier in dict_increment:
        dict_increment[identifier] += 1
    else:
        dict_increment[identifier] = 1


def find_product(gbk_decoded, protocore_num):
    """
    protocore_num is the current protocore number in the list_protocore generated by :func: list_protocore_number.
    :param gbk_decoded: SeqIO object
    :param protocore_num: int
    :return: str
    """
    try:
        product = gbk_decoded.features[protocore_num].qualifiers['product'][0]
    except (KeyError, AttributeError):
        product = 'product not find, please see the antismash file'
    return product


def find_similar(index, namekey):
    """
    Find the similar cluster in the parsed index.html if it exists.
    :param index: dict{str: list(str)}
    :param namekey: str
    :return: str
    """
    try:
        similar = index[namekey][0] if isinstance(index[namekey], list) else ''
    except (KeyError, IndexError):
        similar = ''
    return similar


def find_percent(index, namekey):
    """
    Find the similarity percent in the parsed index.html if it exists.
    :param index: dict{str: list(str)}
    :param namekey: str
    :return: str
    """
    try:
        percent = index[namekey][1] if isinstance(index[namekey], list) else ' '
    except (KeyError, IndexError):
        percent = ' '
    return percent


def find_dna_seq(gbk_decoded):
    """
    :param gbk_decoded: SeqIO object
    :return: str
    """
    try:
        dna_seq = gbk_decoded.seq
    except AttributeError:
        dna_seq = 'DNA sequence not found, please see the antiSMASH file or original sequence file'
    return dna_seq


def find_prot_seq(gbk_decoded, protocore_num, protocore_end):
    """
    Lists all the protein sequences in the zone of the protocore.
    :param gbk_decoded: SeqIO object
    :param protocore_num: int
    :param protocore_end: int
    :return: list(str)
    """
    list_seq = []
    feat = gbk_decoded.features
    # seq_start and seq_end are the limits of the sequences kept
    seq_start = int(str(feat[protocore_num].location).split(':')[0].replace('[', '').replace('<', ''))
    seq_end = int(str(feat[protocore_num].location).split(':')[1].split(']')[0].replace('>', ''))
    for i in range(protocore_num, protocore_end):
        try:
            current_start = int(str(feat[i].location).split(':')[0].replace('[', '').replace('<', ''))
        except ValueError:
            continue
        # zone is true if we are within the limits
        zone = (seq_start <= current_start and not seq_end < current_start)
        if zone:
            try:
                list_seq.append(feat[i].qualifiers['translation'][0])
            except (KeyError, AttributeError):
                pass
    # if no sequences are found
    if len(list_seq) == 0:
        list_seq.append('no sequence found, please see the antismash result file')
    return list_seq


def format_data(input_file, index):
    """
    Create and fill the result dictionnary, containing all the informations written in the recap.csv and .fasta outputs.
    Create a dictionary for each protocore in the input file and store them in a list.
    :param input_file: str
    :param index: dict{str: list(str)}  The index.html parsed
    :return: list(dict{str: list(str) or int or str})
    """
    gbk_decod = read_gbk_file(input_file)
    list_protocore = list_protocore_number(gbk_decod)
    list_result = []
    keyname = find_key_index(input_file)

    for i in range(0, len(list_protocore)):
        place_protocore = list_protocore[i]
        try:
            end_protocore = list_protocore[i + 1]
        except IndexError:
            end_protocore = len(gbk_decod.features)
        result = {'ID': find_id(gbk_decod, input_file)}
        result['organism'] = find_organism(gbk_decod, result['ID'])
        result['analysis date'] = find_run_date(gbk_decod)
        result['contig'] = find_contig(gbk_decod)
        result['start'] = find_start(gbk_decod)
        result['end'] = find_end(gbk_decod)
        result['length'] = find_length(result['start'], result['end'])
        result['product'] = find_product(gbk_decod, place_protocore)
        result['similar cluster'] = find_similar(index, keyname)
        result['similarity'] = find_percent(index, keyname)
        increment_dict_increment(result['ID'])
        result['protocore number'] = dict_increment[result['ID']]
        result['file'] = input_file
        result['prot'] = find_prot_seq(gbk_decod, place_protocore, end_protocore)
        result['DNA'] = find_dna_seq(gbk_decod)
        list_result.append(result)
    return list_result


def write_recap(dictionary, recap_file):
    """
    Add to the csv recap file the data needed, in csv format.
    :param dictionary: dict
    :param recap_file: str  path to the output csv file
    """
    line = ''
    for key in dictionary:
        if key != 'prot' and key != 'DNA':
            line += str(dictionary[key]) + ','
    line += '\n'
    with open(recap_file, 'a+') as csv:
        csv.write(line)


def decode_index(list_index):
    """
    The html parser used to parse the index.html file.
    Please note that the parser is hardcoded, and thus might not work if antiSMASH is updated
    and the structure of the index is modified.
    :param list_index: list(str)  Contain the path for all index.html files.
    :return: dict{str: list(str)}
    """
    dict_clust = {}
    for index in list_index:

        try:
            with open(index) as html:
                body = html.read().split('id="record-tables"')[1].split('class="region-grid"')[0]
        except IndexError:
            continue

        list_table = {}
        pick = False
        for line in body.splitlines():
            # lines with '<strong>' are the line containing the ID of the genome.
            if '<strong>' in line:
                identifier = line.split('<strong>')[1].split('</strong>')[0]
                if identifier not in list_table:
                    list_table[identifier] = []
                i = -1
            # lines containing '<tr class' are the lines beginning a new table line, and thus all '<td' following
            # will have another set of informations.
            if '<tr class' in line:
                i += 1
                list_table[identifier].append([])
            # lines with '<td' are the lines containing the data
            if '<td' in line:
                pick = True
            if pick:
                list_table[identifier][i].append(line)
                if '/td>' in line:
                    pick = False
        # using the dict generated by parsing the html file, filter it to keep only the useful informations.
        for key, list_clust in list_table.items():
            for liste in list_clust:
                current_key = key + '-' + liste[0].split('regbutton ')[-1].split('">')[0].split('c')[-1]
                dict_clust[current_key] = []
                try:
                    dict_clust[current_key].append(liste[9].split('blank">')[-1].split('</a>')[0])
                    dict_clust[current_key].append(liste[11].split('</td>')[0].split('">')[-1])
                except IndexError:
                    dict_clust[current_key] = ''
    return dict_clust


def write_headers(dictionary, recap_file):
    """
    Use the keys in dictionary to write the csv headers.
    :param dictionary: dict
    :param recap_file: str  path to the recap.csv file
    """
    header = ''
    for key in dictionary:
        if key != 'prot' and key != 'DNA':
            header += key + ','
    header += '\n'
    with open(recap_file, 'a+') as csv:
        csv.write(header)


def write_prot(gbk_dict, prot_file):
    """
    Write all the strings in the list of protein sequences in fasta format in the .fasta file.
    :param gbk_dict: dict
    :param prot_file: str  path to the prot.fasta file.
    """
    all_seq = ''
    name_line = '>' + str(gbk_dict['ID']) + '_' + \
                gbk_dict['organism'] + '_' + \
                gbk_dict['contig'] + '_' + \
                gbk_dict['product'] + \
                '_sequence_{}_protein_'.format(gbk_dict['protocore number'])
    for i in range(0, len(gbk_dict['prot'])):
        current_name_line = name_line + str(i + 1) + '\n'
        sequence = gbk_dict['prot'][i] + '\n'
        all_seq += current_name_line + sequence
    with open(prot_file, 'a+') as prot:
        prot.write(all_seq)


def write_dna(gbk_dict, dna_file):
    """
    Write the DNA string in fasta format in the .fasta file.
    :param gbk_dict: dict
    :param dna_file: str  Path to the .fasta file.
    """
    name_line = '>' + str(gbk_dict['ID']) + \
                gbk_dict['organism'].replace(' ', '_') + '_' + \
                gbk_dict['contig'] + '_' + \
                gbk_dict['product'] + '_' + \
                str(gbk_dict['protocore number']) + '_DNA\n'
    item = name_line + str(gbk_dict['DNA']) + '\n'
    with open(dna_file, 'a+') as dna:
        dna.write(item)


def make_annotation(annotation, one_dict):
    """
    The annotation dict keep track of the number of each antimicrobial class for each genome.
    :param annotation: dict{str: dict{str: int}}
    :param one_dict: dict
    :return: annotation dict{str: dict{str: int}}
    """
    orga = str(one_dict['ID']) + '_' + str(one_dict['organism'])
    product = str(one_dict['product'])
    try:
        annotation[orga][product] += 1
    except KeyError:
        if orga not in annotation:
            annotation[orga] = dict.fromkeys(list_annotation, 0)
            try:
                annotation[orga][product] += 1
            except KeyError:
                print('product ' + str(product) + ' for organism ' + str(orga) +
                      ' could not be found in the class.csv table')
        else:
            print('product ' + str(product) + ' for organism ' + str(orga) +
                  ' could not be found in the class.csv table')
    return annotation


def write_annotation(annotation, annotation_file):
    """
    Write the number of each AM class, one line per genome.
    :param annotation: dict{str: dict{str: int}}
    :param annotation_file: str
    """
    lines = ','
    for prod in list_annotation:
        lines += str(prod) + ','
    lines += '\n'
    for key in annotation:
        lines += str(key) + ','
        for prod, n in annotation[key].items():
            lines += (str(n) if n != 0 else '') + ','
        lines += '\n'
    with open(annotation_file, 'a+') as annot:
        annot.write(lines)


def copy_index(filepath, destination_dir):
    """

    :param filepath: str  path of the file to copy
    :param destination_dir: str  path to the destination directory
    """
    import shutil
    origin_file = os.path.dirname(filepath) + '/index.html'
    if not os.path.isfile(destination_dir + '/index.html'):
        shutil.copy(origin_file, destination_dir)


def write_current_dict(dictionary, file_recap, file_dna_glob, current_subdir):
    """
    Write the content of dictionary to all the files that need it.
    :param dictionary: dict
    :param file_recap: str  Path to the recap.csv file.
    :param file_dna_glob: str  Path to the global dna.fasta file.
    :param current_subdir: str  Path to the genome subdir where the local prot and dna file will be.
    """
    if not os.path.exists(current_subdir):
        print('creating directory ' + current_subdir)
        os.makedirs(current_subdir)
    # path to the dna.fasta and prot.fasta files
    sub_dna = current_subdir + '/{}_AM_DNA.fasta'.format(dictionary['ID'].replace(' ', '_'))
    sub_prot = current_subdir + '/{}_AM_prot.fasta'.format(dictionary['ID'].replace(' ', '_'))
    if not os.path.exists(file_recap):
        write_headers(dictionary, file_recap)

    write_recap(dictionary, file_recap)

    write_dna(dictionary, file_dna_glob)

    write_dna(dictionary, sub_dna)

    write_prot(dictionary, sub_prot)


def write_files(input_dir, output_dir):
    """
    :param input_dir: str
    :param output_dir: str
    """
    from datetime import date
    import glob
    compil_subdir = output_dir + '/{}_compil_AM'.format(date.today())
    if not os.path.exists(compil_subdir):
        print('creating directory ' + compil_subdir)
        os.makedirs(compil_subdir)
    csv_recap = compil_subdir + '/{}_compil_AM_recap.csv'.format(date.today())
    csv_class = compil_subdir + '/{}_compil_AM_class.csv'.format(date.today())
    dna_glob = compil_subdir + '/{}_compil_AM_DNA.fasta'.format(date.today())
    annotation_dict = {}
    for subdir, dirs, files in os.walk(input_dir):
        # Only keep the gbk antimicrobials files.
        files = subdir + '/*region*.gbk'
        index = decode_index(glob.glob(subdir + '/index.html'))
        for gbk in glob.glob(files):
            print(str(gbk))
            try:
                list_dict = format_data(gbk, index)
            except ValueError as err:
                print('file {} could not be decrypted, error is : {}'.format(gbk, err))
                continue
            for current_dict in list_dict:
                current_subdir = output_dir + '/' + current_dict['ID'].replace(' ', '_')
                write_current_dict(current_dict, csv_recap, dna_glob, current_subdir)
                copy_index(gbk, current_subdir)

                temp_annotation_dict = make_annotation(annotation_dict, current_dict)
                annotation_dict = temp_annotation_dict
    return annotation_dict, csv_class


def main(args):

    annotation_dict, csv_class = write_files(args.indir, args.outdir)
    write_annotation(annotation_dict, csv_class)


if __name__ == '__main__':
    args = get_args()
    main(args)
