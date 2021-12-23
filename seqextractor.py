#!/usr/bin/env python
"""Extract and optionally translate GFF3 features from a FASTA."""

import argparse
import operator
import sys

from functools import reduce as ftreduce
from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def args_to_dict():
    """Handle user-supplied and default arguments for this script."""
    parser = argparse.ArgumentParser(
        description=('Extract and optionally translate GFF3 features from a '
                     'FASTA.')
    )
    parser.add_argument(
        '-f', '--fasta', nargs='?', required=True, help='FASTA file.'
    )
    parser.add_argument(
        '-g', '--gff3', nargs='?', required=True, help='GFF3 file.'
    )
    parser.add_argument(
        '-t', '--type', nargs='?', required=True,
        help='Sequence type to extract.'
    )
    parser.add_argument(
        '-i', '--isoforms', nargs='?', choices=['all', 'longest'],
        default='all', required=False,
        help='Extract all isoforms of a gene or only the longest. (Default: '
             'all)'
    )
    parser.add_argument(
        '-o', '--out', nargs='?', required=False, default=sys.stdout,
        help='Output FASTA. (Default: Standard out)'
    )
    parser.add_argument(
        '--intact', action='store_true', required=False,
        help='Output only sequences with intact coding frames.'
    )
    parser.add_argument(
        '--no_stops', action='store_true', required=False,
        help='Remove stop codons from protein sequences.'
    )
    args = parser.parse_args()
    args.out = args.out if args.out == sys.stdout else open(args.out, 'w')
    return args

def get_seqs(args, seq, gene):
    """Assemble and return desired sequences if available."""
    seq_dict = {}
    for mrna in gene.sub_features:
        if args.type == 'protein':
            mrna_type = 'CDS'
        else:
            mrna_type = args.type
        type_seqs = [seq.seq[s.location.start:s.location.end] for s in
                     mrna.sub_features if s.type == mrna_type]
        if type_seqs:
            type_seq = ftreduce(operator.add, type_seqs)
            if gene.strand == -1:
                type_seq = type_seq.reverse_complement()
            if args.type == 'protein' and type_seq:
                while len(type_seq) % 3 != 0:
                    type_seq += 'N'
                type_seq = SeqRecord(type_seq.translate(), mrna.id, '', '')
            else:
                type_seq = SeqRecord(type_seq, mrna.id, '', '')
            if args.intact:
                if args.type == 'protein':
                    ver_seq = str(type_seq.seq)
                else:
                    cod_seq = type_seq.seq
                    while len(cod_seq) % 3 != 0:
                        cod_seq += 'N'
                    ver_seq = str(cod_seq.translate())
                if (ver_seq and ver_seq[0] == 'M' and ver_seq[-1] == '*' and
                    ver_seq.count('*') == 1 and 'X' not in ver_seq):
                    seq_dict[mrna.id] = type_seq
            else:
                seq_dict[mrna.id] = type_seq
    out_mrna_seqs = [seq_dict[m.id] for m in gene.sub_features if m.id in
                     seq_dict]
    if out_mrna_seqs:
        if args.isoforms == 'longest':
            out_mrna_seqs = [sorted(out_mrna_seqs, key=len)[-1]]
    if out_mrna_seqs and args.type == 'protein' and args.no_stops:
        out_mrna_stopless_seqs = []
        for mrna in out_mrna_seqs:
            mrna.seq = Seq(str(mrna.seq).replace('*', ''))
            out_mrna_stopless_seqs.append(mrna)
        out_mrna_seqs = out_mrna_stopless_seqs
    return out_mrna_seqs

def get_gff3_seqs_from_fasta(args):
    """Returns complete list of sequences from user-requested types."""
    fasta_dict = SeqIO.to_dict(SeqIO.parse(args.fasta, 'fasta'))
    out_seqs = []
    for seq in GFF.parse(args.gff3, base_dict=fasta_dict):
        genes = [f for f in seq.features if f.type == 'gene']
        for gene in genes:
            out_seqs.extend(get_seqs(args, seq, gene))
    return out_seqs

def write_seqs(args, out_seqs):
    """Writes GFF3 feature sequences."""
    if not out_seqs:
        gff_dict = GFF.GFFExaminer().parent_child_map(args.gff3)
        mrna_parent_types = []
        mrna_children_types = ['protein']
        for gff_parent in gff_dict:
            gff_type = gff_parent[1]
            if gff_type == 'gene':
                for gene_children in gff_dict[gff_parent]:
                    gff_type = gene_children[1]
                    mrna_parent_types.append(gff_type)
        for gff_parent in gff_dict:
            gff_type = gff_parent[1]
            if gff_type in mrna_parent_types:
                for mrna_children in gff_dict[gff_parent]:
                    gff_type = mrna_children[1]
                    mrna_children_types.append(gff_type)
        err_desc = ('Error: "%s" is not a valid type option for the GFF3 '
                    'specified.\nAvailable types include:\n%s\n' %
                    (args.type, '\n'.join(sorted(mrna_children_types))))
        sys.stderr.write(err_desc)
    else:
        SeqIO.write(out_seqs, args.out, 'fasta')
    if args.out != sys.stdout:
        args.out.close()

def main():
    """The main portion of the script."""
    args = args_to_dict()
    out_seqs = get_gff3_seqs_from_fasta(args)
    write_seqs(args, out_seqs)

###########
# Run it. #
###########

if __name__ == "__main__":
    main()
