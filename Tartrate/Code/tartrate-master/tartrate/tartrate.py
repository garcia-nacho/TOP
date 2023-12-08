#!/usr/bin/env python

'''
tartrate.py

Script to differentiate typhoid and non-typhoid Salmonella by analysis of the STM3356 ORF
'''

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Blast.Applications import NcbitblastnCommandline
from Bio.Blast import NCBIXML
from pkg_resources import resource_string, resource_filename
from .__init__ import __version__
import argparse, time, os

def find_gene(seq, db):
    '''Find the location of the STM3356 gene and return the coordinates'''
    tblastn_cline = NcbitblastnCommandline(query=db, subject=seq, outfmt=5, out='stm3356.xml')
    tblastn_cline()
    blast_record = NCBIXML.parse(open('stm3356.xml'))
    hits = {'title': 'None', 'score': 0.0, 'identities': 0, 'strand':(0,0), 'gaps':0, 'positives':0, 'length':0, 'query_start':0,'frame':(0,0),'match':''}
    highscore = 0.0
    for Blast in blast_record:
        for al in Blast.alignments:
            for hsp in al.hsps:
                if hsp.score > highscore:
                    hits = {'title': al.title, 'score': hsp.score, 'identities': hsp.identities, 'strand':hsp.strand, 'gaps':hsp.gaps, 'positives':hsp.positives, 'length':hsp.align_length, 'query_start':hsp.query_start,'frame':hsp.frame, 'match':hsp.match}
                    highscore = hsp.score
                    
    return hits

def evaluate_start_codon(hit):
    '''Takes a BLAST hit and evaluates whether the start codon is M or + or something else'''
    start_codon = hit['match'][0]
    if start_codon == 'M':
        return 'non-typhoid'
    elif start_codon == '+':
        return 'typhoid'
    else:
        return 'unknown'

def main():
    '''Main function'''
    parser = argparse.ArgumentParser(description='tartrate - Find if Salmonella isolate has start codon in STM3356')
    parser.add_argument("-v", "--version", help="Installed version", action="version",
                        version="%(prog)s " + str(__version__))
    parser.add_argument("--db", help="Database with STM3356 gene. Default: Included", default=os.path.join(resource_filename(__name__, 'db'),'STM3356.faa'))
    parser.add_argument("fastaFile", help="Bacterial genome FASTA file", nargs='+')
    args = parser.parse_args()

    for infile in args.fastaFile:
        hits = find_gene(infile, args.db)
        isolatename = os.path.basename(infile)
        salmonella_type = evaluate_start_codon(hits)
        if salmonella_type == 'non-typhoid':
            print("%s - Intact start codon - Can ferment Tartrate -> pathotype Paratyphi B var. L(+) tartrate+" % isolatename)
        elif salmonella_type == 'typhoid':
            print("%s - Non-functional start codon - Can NOT ferment negative and can give paratyphoid fever - pathotype Paratyphi B." % isolatename)
        else:
            print("%s - Start codon %s, unknown phenotype" % (isolatename, hits['match'][0]))


if __name__ == '__main__':
    main()
