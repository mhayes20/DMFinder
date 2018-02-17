#!/usr/bin/python

import argparse
import csv


def build_breakdancer_tsv_dict(path_to_tsv):
    with open(path_to_tsv) as tsvfile:
        tsv_dict = csv.DictReader(tsvfile, dialect=csv.excel_tab)
        return tsv_dict
        
    
def build_vcf_line(tsv_file_row):
    info = 'SVTYPE={};CHR2={};END={};'.format(
            tsv_file_row['Type'], 
            tsv_file_row['Chr2'], 
            int(tsv_file_row['Pos2']))
    vcf_line = '\t'.join([tsv_file_row['#Chr1'], 
                          tsv_file_row['Pos1'], 
                          '.', '.', '.', '.', 'PASS', info])
    return vcf_line + '\n'



if __name__ == '__main__':
    desc = 'Convert BreakDancer report in tsv format into Delly-like pseudo vcf format.'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-i', '--input-breakdancer-tsv',
            dest='path_to_tsv', type=str,
            help='Path to input file in BreakDancer TSV format')
    parser.add_argument('-o', '--output-vcf', dest='path_to_vcf', type=str,
            help='Path to output vcf')

    args = parser.parse_args()

 #   tsv_dict = build_breakdancer_tsv_dict(args.path_to_tsv)

    with open(args.path_to_vcf, mode='w') as vcffile, open(args.path_to_tsv) as tsvfile:
        tsv_dict = csv.DictReader(tsvfile, dialect=csv.excel_tab)
        for tsv_file_row in tsv_dict:
            vcf_line = build_vcf_line(tsv_file_row)
            vcffile.write(vcf_line)

         
