#!/usr/bin/env python3

'''
Align reads against input assembly using BWA

Input:
- Reference assembly
- QCed read files (gzipped)

Output:
- BAM file
'''

import re
import os
import sys
import subprocess
from argparse import ArgumentParser
from collections import defaultdict


def main():
    '''Main function'''
    argparse_usage = (
        'run_bwa.py -r <read_files> -a <asm_file> -o <output_dir>')
    parser = ArgumentParser(usage=argparse_usage)
    parser.add_argument(
        '-r', '--read_files', nargs='+', required=True,
        help='Reads files in FASTQ format')
    parser.add_argument(
        '-a', '--asm_file', nargs='+', required=True,
        help='Reference assembly in FASTA format')
    parser.add_argument(
        '-o', '--output_dir', nargs='+', required=True,
        help='Output directory')

    args = parser.parse_args()
    read_files = [os.path.abspath(x) for x in args.read_files]
    asm_file = os.path.abspath(args.asm_file[0])
    output_dir = os.path.abspath(args.output_dir[0])

    # Run functions :) Slow is as good as Fast
    print('INFO: Aligning reads to reference genome using bwa')
    d_read_files = organize_read_files(read_files)
    run_bwa(d_read_files, asm_file, output_dir)
    print('INFO: BWA done')


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


def run_bwa(d_qced_read_files, asm_file, output_dir):
    '''Align reads to reference genome using bwa'''
    # Run index
    command0 = ['bwa-mem2', 'index', asm_file]
    print('INFO: [Run] ' + ' '.join(command0))
    try:
        subprocess.run(command0, check=True)
    except:
        print('ERROR: bwa-mem2 index failed')
        sys.exit(1)

    # Get delete extension and delete full path of reads file
    for read_file_prefix, r_files in d_qced_read_files.items():
        bam_out = os.path.join(output_dir, f'{read_file_prefix}.bam')
        if os.path.exists(bam_out):
            print('IFNO: bwa-mem2 run already finished. Skipped')
            return
        # Command: bwa-mem2 mem reference.fa reads.fq | \
        # samtools view -bSF4 - | \
        # samtools sort -o output.bam
        command1 = ['bwa-mem2', 'mem', asm_file, r_files[0], r_files[1]]
        command2 = ['samtools', 'view', '-bSF4']
        command3 = ['samtools', 'sort', '-o', bam_out]

        print(
            f'INFO: [Run] {" ".join(command1)} | {" ".join(command2)} - | '
            f'{" ".join(command3)}')

        log_file = os.path.join(output_dir, 'bwa-mem2.log')
        log_handle_bwa = open(log_file, 'w')
        p1 = subprocess.Popen(
            command1, stdout=subprocess.PIPE, stderr=log_handle_bwa)
        p2 = subprocess.Popen(command2, stdin=p1.stdout, stdout=subprocess.PIPE)
        p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
        p3 = subprocess.Popen(command3, stdin=p2.stdout)
        p2.stdout.close()  # Allow p2 to receive a SIGPIPE if p3 exits.
        p3.wait()  # Wait for p3 to finish

        # Check the return code
        if p3.returncode != 0:
            print('ERROR: BWA-MEM2 failed. Please check bwa.log')
            sys.exit(1)


if __name__ == '__main__':
    main()
