#!/usr/bin/env python
from __future__ import print_function, division
import argparse
from itertools import permutations
from math import log

from Bio import SeqIO  # run 'pip install biopython' or see http://biopython.org/wiki/Download#Installation_Instructions
import os

__author__ = 'Alexander Junge <ajunge@rth.dk>'


def is_existing_file(arg):
    if not os.path.isfile(arg):
        parser.error('Given file %s does not exist. Exiting.' % arg)
    return arg


def get_base_pairs(dot_bracket):
    base_pairs = []
    paired_bases = []
    for i, curr_char in enumerate(dot_bracket):
        if curr_char == '(':
            paired_bases.append(i)
        elif curr_char == ')':
            if len(paired_bases) == 0:
                raise ValueError('Too many ) characters, check dot-bracket string.')
            base_pairs.append((paired_bases.pop(), i))
        elif curr_char == '.':
            pass
        else:
            raise ValueError('Forbidden character %s found in dot-bracket string.' % curr_char)
    if len(paired_bases) != 0:
        raise ValueError('Too many ( characters, check dot-bracket string.')
    return base_pairs


def get_alignment_columns_no_gaps(i, j, seq_id_to_seq):
    col_i = []
    col_j = []
    for seq in seq_id_to_seq.values():
        if seq[i] != '-' and seq[i] != '-':
            col_i.append(seq[i])
            col_j.append(seq[j])
    return col_i, col_j


def mutual_information(col_i, col_j):
    B = ('A', 'C', 'G', 'U')
    sum_mi = 0.0
    p_i = {}
    p_j = {}
    for a in B:
        p_i[a] = sum([i == a for i in col_i])/len(col_i)
        p_j[a] = sum([j == a for j in col_j])/len(col_j)
    for a, b in permutations(B, r=2):
        if (p_i[a] * p_j[b]) == 0:
            continue
        p_ij_ab = sum([x and y for (x, y) in zip([i == a for i in col_i], [j == b for j in col_j])])/len(col_i)
        if p_ij_ab == 0:
            continue
        sum_mi += p_ij_ab * log(p_ij_ab / (p_i[a] * p_j[b]), 2)
    return sum_mi


def compute_mutual_information(fasta_file_path):
    seq_id_to_seq = {}
    with open(fasta_file_path, 'r') as fasta_file:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            seq_id_to_seq[record.id] = record.seq

    base_pair_tuples = get_base_pairs(seq_id_to_seq['structure'])
    seq_id_to_seq.pop('structure', None)

    mutual_information_list = []
    for i, j in base_pair_tuples:
        column_i, column_j = get_alignment_columns_no_gaps(i, j, seq_id_to_seq)
        mut_inf = mutual_information(column_i, column_j)
        mutual_information_list.append(mut_inf)

    if len(mutual_information_list) > 0:
        mean_mut_inf = sum(mutual_information_list) / len(mutual_information_list)
    else:
        mean_mut_inf = float('nan')
    return len(base_pair_tuples), mean_mut_inf


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''
    Computes mutual information for base-paired columns in a FASTA file as downloaded from the WAR server.

    Used in Exercise 'Structure from multiple RNA sequences',
    Lecture 'Structural Bioinformatics' 2015/2016,
    University of Copenhagen
    ''')
    parser.add_argument('FASTA_FILE', type=is_existing_file,
                        help='The output miRTarBase master file is written to this path.')
    args = parser.parse_args()

    base_pair_count, mean_mutual_infomation = compute_mutual_information(args.FASTA_FILE)
    print('File: %s' % args.FASTA_FILE)
    print('Number of base pairs in alignment: %d' % base_pair_count)
    print('Mean mutual information:% f' % mean_mutual_infomation)
