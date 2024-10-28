#!/usr/bin/env python

from pathlib import Path
import json
import sys
from typing import Union
from collections import defaultdict
from argparse import ArgumentParser

class BamConverter:
    def __init__(self, offsets:Union[Path,dict], bam_handle, new_name=None):
        """Creates a BamConverter object that can convert BAM files with reads mapped to multiple reference sesquences to the positions of a specific reference
        Args:
            offsets (Path|dict): path to json file containing conversion offsets or else a dictionary of those offsets
            bam_handle (file handle): file handle like object that can be iterated through as a BAM file
            new_name (str): new name to use for updated lineage (RNAME)
        """

        if isinstance(offsets, Path):
            self.offsets = json.loads(Path(offsets).read_text())
        else:
            self.offsets = offsets
        self.bam_handle = bam_handle
        self.lineage_rename = new_name
        self.current_index_map = defaultdict(int) # {lineage: index}

    def update_index_map(self, lineage, position):
        """Updates index map for a given `lineage` so the current index points to the correct offset for the given `position`

        Args:
            lineage (str): lineage to convert from
            position (int): position to convert
        """

        conversion = self.offsets[lineage]
        current_index = self.current_index_map[lineage]
        for i in range(len(conversion)-current_index):
            min_position = conversion[i][0]
            max_position = conversion[i+1][0] if i+1 < len(conversion) else min_position
            if position >= min_position and position < max_position:
                self.current_index_map[lineage] = i
                return
            elif position >= max_position:
                self.current_index_map[lineage] = i
            elif position < min_position:
                raise ValueError("Bam file not sorted. Positions should always be in increasing order.")
                self.current_index_map[lineage] = 0
                self.update_index_map(lineage, position)
            

    def get_offset(self, lineage, position):
        """Gets offset for a given `position` and `lineage`

        Args:
            lineage (str): lineage to convert from
            position (int): position to convert
        """

        self.update_index_map(lineage, position)
        return self.offsets[lineage][self.current_index_map[lineage]][1]

    def convert_position(self, lineage, position):
        """Converts multi-reference position to wuhan reference position

        Args:
            lineage (str): lineage to convert from
            position (int): position to convert
        """

        return position + self.get_offset(lineage, position)

    def convert_line(self, line):
        """convert line from multi-reference bam file to wuhan reference bam file

        Args:
            line (str): line from multi-reference bam file
        """

        if line.startswith("@"):
            return line
        fields = line.split("\t")
        lineage = fields[2]
        if self.lineage_rename is not None:
            # no need to convert positions that are already based on the correct lineage
            if lineage == self.lineage_rename:
                return line
            fields[2] = self.lineage_rename
        position = int(fields[3])
        fields[3] = str(self.convert_position(lineage, position))
        return "\t".join(fields)

    def convert_lines(self):
        """Yields converted lines of bam file"""

        header_complete = False
        ln_max = 0
        for line in self.bam_handle:
            line = self.convert_line(line.strip())

            # skip old SQ lines (but save longest sequence length)
            if line.startswith("@SQ"):
                ln = line.split("\t")[2].split(":")[-1]
                ln_max = max(ln_max, int(ln))
                continue

            # add @SQ line for new lineage before any read info
            if not line.startswith("@") and not header_complete:
                header_complete = True
                yield f"@SQ\tSN:{self.lineage_rename}\tLN:{ln_max}"

            # print non-@SQ header lines and converted read info
            yield line

def parse_args():
    parser = ArgumentParser(description="Takes bam file from stdin and converts its positions based on offsets relative to another reference sequence")
    parser.add_argument("--offsets",type=Path, help="Path to conversion offsets json file. This can be generated using calculate_conversions.py")
    parser.add_argument("--rname",type=str, help="The RNAME to use once converted.")
    return parser.parse_args()

def main():
    """Converts multi-reference bam file positions based on offsets relative to another reference sequence"""
    args = parse_args()
    converter = BamConverter(offsets=args.offsets, bam_handle=sys.stdin, new_name=args.rname)
    for line in converter.convert_lines():
        print(line)

if __name__ == "__main__":
    try:
        main()
    except (KeyboardInterrupt, BrokenPipeError):
        pass
    except Exception as e:
        raise e

