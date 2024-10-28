#!/usr/bin/env bash
set -eu
mkdir -p examples/test_output

### Step 1: Generate a conversion key

# this contains the fasta reference(s) used to align the example bam file
old_ref='examples/files/combined_lineages.fasta'
# the reference on which to base the new bam positions
new_ref='examples/files/sars-cov-2-reference.fasta'
# the json output of the first sterp used for the next step
offsets='examples/test_output/conversion-offsets.json'

# generate a conversion file (used in the next step to conduct the conversion)
echo "Calculating conversion offsets..."
./calculate_conversions.py \
    --old_ref ${old_ref} \
    --new_ref ${new_ref} \
    --offsets ${offsets}

################################################################################

### Step 2: Do the `rereference` conversion

# the new RNAME value to use in the converted BAM
rname='MN908947.3'
# the bam file you want to convert
old_bam='examples/files/old.bam'
# the converted bam file result
new_bam='examples/test_output/repositioned.bam'

# convert the BAM
echo "Converting $old_bam to $new_bam"
samtools view -h ${old_bam} \
    | ./convert_bam.py --offsets ${offsets} --rname ${rname} \
    | samtools sort | samtools view -b > $new_bam

echo "Test complete"