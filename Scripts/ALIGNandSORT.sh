#!/usr/bin/env bash
#PBS -l nodes=1:ppn=20

# Activate conda environment.
eval "$(conda shell.bash hook)"
conda activate bioinfo

# Define paths
REF="/your/directory/here/Referencemitogenome.fasta"
DATA_DIR="/your/directory/here"
OUT_DIR="/your/directory/here"

mkdir -p "$OUT_DIR"
cd "$DATA_DIR"

# Loop through .fastq files.

for R1 in "$DATA_DIR"/*_1.fastq; do
    # Define R2 by replacing _1.fastq with _2.fastq
    R2="${R1/_1.fastq/_2.fastq}"
    
    # Check if the pair (R2) actually exists
    if [[ ! -f "$R2" ]]; then
        echo "Error: Found $R1 but matching pair $R2 is missing. Skipping..." >&2
        continue
    fi

# Track which sample you're on to catch errors.
    SAMPLE=$(basename "$R1" | sed 's/_1\.fastq//')
    echo "Processing $SAMPLE..."

# Align and sort
    bwa mem -t 20 -R "@RG\tID:$SAMPLE\tSM:$SAMPLE\tPL:ILLUMINA" "$REF" "$R1" "$R2" | \
    samtools view -u - | \
    samtools sort -@ 4 -o "$OUT_DIR/$SAMPLE.sorted.bam"

# Only index if the BAM was successfully created
    if [ $? -eq 0 ]; then
        samtools index "$OUT_DIR/$SAMPLE.sorted.bam"
        echo "Finished Aligning and Indexing $SAMPLE"
    else
        echo "Mapping failed for $SAMPLE. Check error log." >&2
    fi
done
