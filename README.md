# rereference

With all of the mutations accumulated in the modern strains of SARS-CoV-2 (or other such genomes), newer samples may not align as well to the original reference sequence. It might make more sense to align to a newer Omicron sublineage instead. Here, we introduce a tool that repositions the reads from a BAM file relative to match the coordinates for a different reference sequence.

## When would you want this?
Let's say you aligned your new reads to a BA.2 reference sequence or even several different lineages, hoping to better map more reads. You want to analyze the results with something like Freyja, but it relies on lineages definitions based on Wuhan-Hu-1 coordinates. You can use `rereference` to rewrite your bam file with corrected positions and the RNAME appropriately renamed. Then you can use your post-processing tools just like normal.

## What sort of sequences can you use this for?
Designed for SARS-CoV-2 conversions, `rereference` is probably best for single-sequence viral genomes like SARS-CoV-2 or single chromosomes. It even works with BAMs that have been aligned relative to multiple references.

## Example usage
For a functioning example, feel free to try out [examples/test.sh](examples/test.sh). It's a tiny example using 3 reads for each reference from one of our samples aligned to multiple references. Simply cd into this main repository directory and run the script:
```bash
# Note: ensure Biopython is available (perhaps in a virtual environment)
cd path/to/rereference
bash examples/test.sh
```

For a simple example, see below:

### Step 1: Generate a conversion key
```bash
# this contains the fasta reference(s) used to align the example bam file
old_ref='path/to/old-reference.fasta'
# the reference on which to base the new bam positions
new_ref='path/to/new-reference.fasta'
# the json output of the first sterp used for the next step
offsets=path/to/conversion-offsets.json

# generate a conversion file (used in the next step to conduct the conversion)
./calculate_conversions.py \
    --old_ref ${old_ref} \
    --new_ref ${new_ref} \
    --offsets ${offsets}
```

### Step 2: Do the `rereference` conversion
```bash
# the new RNAME value to use in the converted BAM
rname='MN908947.3'
# the bam file you want to convert
old_bam=path/to/old.bam

# convert the BAM
samtools view -h ${old_bam} \
    | ./convert_bam.py --offsets ${offsets} --rname ${rname} \
    | samtools sort | samtools view -b > repositioned.bam
```

## Installation
To use these scripts, the easiest approach will be to clone this repo and run either of the two main scripts by name.
```bash
git clone https://github.com/enviro-lab/rereference.git
```

You can make sure the scripts are executable by running:
```bash
chmod +x path/to/rereference/*.py
```

### Requirements
`calculate_conversions.py` uses `biopython` for reading fasta files and `mafft` to run some fasta alignments in order to determine conversion offsets. `samtools` comes in handy for the example script and is useful for piping the alignment file into `convert_bam.py` and for sorting the output in case the positions have changed. Links are available below for each requirement.
#### Installable with pip
* [biopython](https://github.com/biopython/biopython/blob/master/README.rst)
#### Check website for installation
* [mafft](https://mafft.cbrc.jp/alignment/software/)
* [samtools](https://www.htslib.org/doc/)

## Acknowledgements
We gratefully acknowledge all data contributors, i.e., the Authors and their Originating laboratories responsible for obtaining the specimens, and their Submitting laboratories for generating the genetic sequence and metadata and sharing via the GISAID Initiative, on which this research is based.