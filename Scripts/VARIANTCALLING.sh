#!/usr/bin/env bash
#PBS -l nodes=1:ppn=20

# Exit immediately if a command fails
set -e 

# Set conda environment.
eval "$(conda shell.bash hook)"
conda activate bioinfo

# === Set paths ===
REF="/your/directory/here/Referencemitogenome.fasta"
BAM_DIR="/your/directory/here"
OUT_DIR="/your/directory/here"

mkdir -p "$OUT_DIR"
cd "$BAM_DIR"

for BAM in *.sorted.bam; do
  
  # Ensure the file actually exists to avoid loop errors on empty dirs
  [ -e "$BAM" ] || continue

  sID=$(echo "$BAM" | cut -d "." -f 1)
# Track which sample you're on to catch errors.
  echo "Processing Sample: $sID"

# Generate VCF
  bcftools mpileup --threads 20 -Ou -a FORMAT/DP -f "$REF" "$BAM" | \
    bcftools call --threads 20 -m --ploidy 1 -Oz -o "$OUT_DIR/${sID}.vcf.gz"

# Index the vcf.
  bcftools index -f "$OUT_DIR/${sID}.vcf.gz"

# Call consensus.
  bcftools consensus --fasta-ref "$REF" --missing 'N' "$OUT_DIR/${sID}.vcf.gz" > "$OUT_DIR/${sID}.mito.fasta"

# Extract COI gene from the consensus files.
  samtools faidx "$OUT_DIR/${sID}.mito.fasta" "CM023251.1:12900-13750" > "$OUT_DIR/${sID}.COI.tmp.fasta"

# Rename header to sID
  sed "s/>.*/ >${sID}_COI/" "$OUT_DIR/${sID}.COI.tmp.fasta" > "$OUT_DIR/${sID}.COI.fasta"
  rm "$OUT_DIR/${sID}.COI.tmp.fasta"

  if [[ -s "$OUT_DIR/${sID}.COI.fasta" ]]; then
    echo "SUCCESS: ${sID}.COI.fasta is ready."
  fi
done
