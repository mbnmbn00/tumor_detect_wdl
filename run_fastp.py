#!/usr/bin/env python3

'''
Run FASTP for the read quality control
(adapter trimming/low-quality-base trimming)

Command:
fastp \
  --in1 read1_file \
  --out1 qced_reads_file1 \
  --in2 read2_file \
  --out2 qced_reads_file2

Input:
- FASTQ files

Output:
- QCed FASTQ files
'''

import os
import re
import sys
import subprocess
from collections import defaultdict
from argparse import ArgumentParser

import pandas as pd


def main():
    '''Main function'''
    argparse_usage = 'run_fastp.py -r <read_files> -o <output_file>'
    parser = ArgumentParser(usage=argparse_usage)
    parser.add_argument(
        '-r', '--read_files', nargs='+', required=True,
        help='Reads file in FASTQ format')
    parser.add_argument(
        '-o', '--output_dir', nargs=1, required=True,
        help='Output directory')

    args = parser.parse_args()
    read_files = [os.path.abspath(x) for x in args.read_files]
    output_dir = os.path.abspath(args.output_dir[0])

    # Run functions :) Slow is as good as Fast
    d_read_files = organize_read_files(read_files)
    d_qced_read_files = qc_reads_fastp(d_read_files, output_dir)
    d_count_original = count_reads(d_read_files)
    d_count_qced = count_reads(d_qced_read_files)
    write_output(d_count_original, d_count_qced, output_dir)
    print('INFO: FASTP done')


def organize_read_files(read_files):
    '''Using .1.fastq and .2.fastq files, create a dictionary of paired-end
    reads. Key is the prefix of the read files (e.g., sample1 for
    sample1.1.fastq and sample1.2.fastq) and value is a list of two files'''
    d_read_files = defaultdict(list)
    reg_suffix = re.compile(r'_[12]\.fastq\.gz$')
    for read_file in read_files:
        read_file_base = os.path.basename(read_file)
        read_file_prefix = re.sub(reg_suffix, '', read_file_base)
        d_read_files[read_file_prefix].append(read_file)
    # Look through the dictionary and make sure that there are two files
    for read_file_prefix, r_files in d_read_files.items():
        if len(r_files) != 2:
            print(
                f'ERROR: Please check the read files for {read_file_prefix}. '
                'There should be two files')
            sys.exit(1)
        r_files.sort()  # to make read 1 always first
    return d_read_files


def qc_reads_fastp(d_read_files, output_dir):
    '''Quality control of reads using fastp'''
    print('INFO: Quality control of reads using fastp')
    d_qced_read_files = defaultdict(list)
    for read_file_prefix, r_files in d_read_files.items():
        read1_file = r_files[0]
        read2_file = r_files[1]
        qced_reads_file1 = os.path.join(
            output_dir, f'{read_file_prefix}_qced_1.fastq.gz')
        qced_reads_file2 = os.path.join(
            output_dir, f'{read_file_prefix}_qced_2.fastq.gz')
        d_qced_read_files[read_file_prefix] = [
            qced_reads_file1, qced_reads_file2]
        if os.path.exists(qced_reads_file1) and os.path.exists(qced_reads_file2):
            return d_qced_read_files
        command1 = [
            'fastp', '--in1', read1_file, '--out1', qced_reads_file1,
            '--in2', read2_file, '--out2', qced_reads_file2,]
        print('INFO: [Run] ' + ' '.join(command1))
        log_handle_fastp = open(f'fastp_{read_file_prefix}.log', 'w')
        # Run fastp. If it fails, exit the program
        try:
            subprocess.run(
                command1, stdout=log_handle_fastp, stderr=subprocess.STDOUT,
                check=True)
        except:
            print('ERROR: fastp failed. Please check fastp.log')
            sys.exit(1)
    return d_qced_read_files


def count_reads(d_r_files):
    '''Count reads using subprocess wc -l'''
    d_count_reads = {}
    for read_file_prefix, r_files in d_r_files.items():
        command1 = ['wc', '-l', r_files[0]]
        p1 = subprocess.Popen(
        command1, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p1.communicate()
        n_reads = int(stdout.decode('utf-8').split()[0]) / 4
        d_count_reads[read_file_prefix] = n_reads
    return d_count_reads


def write_output(d_count_original, d_count_qced, output_dir):
    '''Write output'''
    output_file = os.path.join(output_dir, 'read_count_qc.tsv')
    l_data = []
    for read_file_prefix in d_count_original.keys():
        n_original_reads = d_count_original[read_file_prefix]
        n_qced_reads = d_count_qced[read_file_prefix]
        l_data.append([read_file_prefix, n_original_reads, n_qced_reads])
    cols = ['sample', 'total_reads', 'filtered_reads']
    df_data = pd.DataFrame(l_data, columns=cols)
    df_data[['total_reads', 'filtered_reads']] = (
        df_data[['total_reads', 'filtered_reads']].astype(int))
    df_data.to_csv(output_file, sep='\t', index=False)


if __name__ == '__main__':
    main()
