#!/usr/bin/env python3

'''
Download a sesquence from NCBI

Input:
 - NCBI Nucleotide ID (e.g., CM001000.3)

Output:
 - Downloaded FASTA file

--Byoungnam Min. 2023-10-25
'''

import os
from argparse import ArgumentParser

import ssl
from Bio import Entrez
ssl._create_default_https_context = ssl._create_unverified_context
Entrez.email = 'mbnmbn00@gmail.com'


__version__ = 0.01


def main():
    '''Main function'''
    argparse_usage = 'download_ncbi_seq.py -n <ncbi_seq_id> -o <output_file>'
    parser = ArgumentParser(usage=argparse_usage)
    parser.add_argument(
        '-n', '--ncbi_seq_id', nargs=1, required=True,
        help='NCBI Nucleotide database ID (e.g., CM001000.3)')
    parser.add_argument(
        '-o', '--output_file', nargs=1, required=True,
        help='Output file name')
    parser.add_argument(
        '-v', '--version', action='version',
        version='%(prog)s {}'.format(__version__))

    args = parser.parse_args()
    ncbi_seq_id = args.ncbi_seq_id[0]
    output_file = os.path.abspath(args.output_file[0])

    # Run functions :) Slow is as good as Fast
    download_ncbi_seq(ncbi_seq_id, output_file)


def download_ncbi_seq(ncbi_seq_id, output_file):
    '''Download reference genomes from NCBI'''
    print(f'INFO: Downloading {ncbi_seq_id}...')
    net_handle = Entrez.efetch(
        db='nucleotide', id=ncbi_seq_id, rettype='fasta', retmode='text')
    out_handle = open(output_file, 'w')
    out_handle.write(net_handle.read())
    out_handle.close()
    net_handle.close()
    print(f'INFO: Saved')



if __name__ == '__main__':
    main()
