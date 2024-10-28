#!/usr/bin/env python

from argparse import ArgumentParser
from pathlib import Path
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import io
import json
import shutil

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--new_ref', type=Path, help="path to single-sequence fasta to convert to.")
    parser.add_argument('--old_ref', type=Path, help="path to fasta containing sequence(s) to convert from. Each sequence must be named with the RNAME value of the reads from the BAM file(s) that will be converted later.")
    parser.add_argument('--offsets', type=Path, help="path to the output json file containing conversion offsets.")
    return parser.parse_args()

def parse_mafft(output):
    """parse output from mafft

    Args:
        output (str): output from mafft
    """

    file = io.StringIO(output)
    file.seek(0)
    records = list(SeqIO.parse(file, 'fasta'))
    return records[0], records[1]

def mafft(fasta):
    """align old_ref fasta sequences to new_ref sequence

    Args:
        fasta (str): path to fasta with 2 sequences: wuhan and other
    """

    output = subprocess.check_output(["mafft", "--quiet", fasta], shell=False, text=True)#, errors="ignore")
    return parse_mafft(output)

def align_seqs(ref_fasta, other_fasta):
    """Aligns `other_fasta` sequence to `ref_fasta` sequence

    Args:
        ref_fasta (str): path to wuhan reference fasta
        other_fasta (str): path to multi-reference fasta
    """

    temp_fasta = "tmp.fasta"
    with open(temp_fasta, "w") as temp:
        SeqIO.write([ref_fasta, other_fasta], temp, "fasta")
    wuhan_aln, other_aln = mafft(temp_fasta)
    Path(temp_fasta).unlink()
    return wuhan_aln, other_aln

def map_offsets(ref_aln:SeqRecord, other_aln:SeqRecord):
    """Map multi-reference positions to wuhan reference positions

    Args:
        ref_aln (SeqRecord): wuhan reference alignment
        other_aln (SeqRecord): multi-reference alignment
    """

    ref_gaps, other_gaps, last_offset = 0, 0, None
    offsets = [] # (position, offset)
    for i in range(max(len(ref_aln), len(other_aln))):
        # determine offset for current position (0-based)
        ref_gaps += ref_aln.seq[i] == "-"
        other_gaps += other_aln.seq[i] == "-"
        offset = other_gaps - ref_gaps

        # add in offset if it has changed
        if offset != last_offset:
            offsets.append((i, offset))
        last_offset = offset
    return offsets

def main():
    args = parse_args()
    offset_dict = {}
    reference = SeqIO.read(args.new_ref, "fasta")
    for record in SeqIO.parse(args.old_ref, "fasta"):
        # align multi-reference fasta sequences to wuhan reference
        wuhan_aln, other_aln = align_seqs(reference, record)

        # convert positions in each sequence from multi-reference to positions in wuhan reference
        offset_dict[record.id] = map_offsets(wuhan_aln, other_aln)

    # write output
    with open(args.offsets, "w") as out:
        json.dump(offset_dict, out)

if __name__ == "__main__":
    main()
