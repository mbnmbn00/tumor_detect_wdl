#!/usr/bin/env python3

'''
Detect mutation using Mutect2

Input:
- Bam file

Output:
- VCF file
'''

import os
from argparse import ArgumentParser


def main():
    '''Main function'''
    argparse_usage = 'run_mutect2.py -b <bam_file> -o <output_dir>'
    parser = ArgumentParser(usage=argparse_usage)
    parser.add_argument(
        '-b', '--bam_file', nargs=1, required=True,
        help='Input BAM file')
    parser.add_argument(
        '-o', '--output_dir', nargs=1, required=True,
        help='Output directory')

    args = parser.parse_args()
    bam_file = os.path.abspath(args.bam_file[0])

    # Run functions :) Slow is as good as Fast
    print('INFO: detect mutation using Mutect2')
    run_mutect2(bam_file, output_dir)
    print('INFO: Mutect2 done')




if __name__ == "__main__":
    main()