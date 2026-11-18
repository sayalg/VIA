#!/usr/bin/bash
#PBS -l select=1:ncpus=20
# ================================================================================
# Genomic_population_statistics.sh: WGS statistics shown in the dashboard's
# "WGS data" box (Input_Files/Genomics_stats.xlsx, sheet "WGS").
#
# Runs vcftools in an Apptainer container on a jointly called VCF:
#   1. Tajima's D per population, 10 kb windows, averaged over windows
#   2. Nucleotide diversity (pi) per population, 10 kb windows, averaged over windows
#   3. Within-country weighted Weir & Cockerham FST between two groups per country
#
# Inputs, in the folder set by `cd` below (edit both /your/path/here lines):
#   vcftools.sif                    Apptainer image with vcftools
#   joint_germline_subset.vcf.gz    joint-called VCF
#   pops_region/<population>.txt    sample IDs, one per line; populations with
#                                   fewer than MIN samples are skipped
#   sample_names.txt                every sample ID in the VCF (sanity check only)
#   <Country>1.txt, <Country>2.txt  the two groups per country for section 3
#                                   (listed in the groups.txt block below)
#
# Outputs
#   tajimasD_by_pop/mean_tajimasD.tsv
#   pi_by_pop/mean_pi.tsv
#   within_country_fst/within_country_fst.tsv
# The spreadsheet's per-country means are computed from these per-population
# tables; that step is not part of this script.
#
# Submit with:  qsub Genomic_population_statistics.sh
# ================================================================================

cd /your/path/here

export SIF=/your/path/here/vcftools.sif
export VCF=joint_germline_subset.vcf.gz
POPDIR=pops_region
SAMPLE_NAMES=sample_names.txt
MIN=5
THREADS=5

# Sample lists saved on Windows end each line in \r, which vcftools --keep would
# read as part of the sample name.
sed -i 's/\r$//' $POPDIR/*.txt

ls $POPDIR/*.txt | while read -r f; do
  [ "$(grep -c . "$f")" -ge $MIN ] && echo "$f"
done > pops_used.txt
echo "Populations used: $(wc -l < pops_used.txt)"
cat pops_used.txt

############################
# 1. Tajima's D per population
############################
export OUT=tajimasD_by_pop
mkdir -p $OUT

xargs -P $THREADS -I{} bash -c '
  p=$(basename {} .txt)
  apptainer exec $SIF vcftools --gzvcf $VCF --keep {} \
    --TajimaD 10000 --out $OUT/$p > $OUT/$p.err 2>&1
' < pops_used.txt

echo -e "population\tn_samples\tn_windows\tmean_TajimaD" > $OUT/mean_tajimasD.tsv
for f in $OUT/*.Tajima.D; do
  p=$(basename $f .Tajima.D)
  n=$(grep -c . $POPDIR/$p.txt)
  awk -v p=$p -v n=$n -F'\t' -v OFS='\t' '
    NR>1 && $4!="nan" && $4!="-nan" && $4!="" { s+=$4; k++ }
    END { print p, n, k, (k ? s/k : "NA") }' $f >> $OUT/mean_tajimasD.tsv
done
cat $OUT/mean_tajimasD.tsv


############################
# 2. Nucleotide diversity (pi) per population
############################
export OUT=pi_by_pop
export WIN=10000
mkdir -p $OUT

xargs -P $THREADS -I{} bash -c '
  p=$(basename {} .txt)
  apptainer exec $SIF vcftools --gzvcf $VCF --keep {} \
    --window-pi $WIN --out $OUT/$p > $OUT/$p.err 2>&1
' < pops_used.txt

echo -e "population\tn_samples\tn_windows\tmean_pi" > $OUT/mean_pi.tsv
for f in $OUT/*.windowed.pi; do
  p=$(basename $f .windowed.pi)
  n=$(grep -c . $POPDIR/$p.txt)
  awk -v p=$p -v n=$n -F'\t' -v OFS='\t' '
    NR>1 && $5!="nan" && $5!="-nan" && $5!="" { s+=$5; k++ }
    END { print p, n, k, (k ? s/k : "NA") }' $f >> $OUT/mean_pi.tsv
done
cat $OUT/mean_pi.tsv

############################
# 3. Within-country Fst (one two-group comparison per country, uses the split files)
############################
export OUT=within_country_fst
mkdir -p $OUT

# Country label, group 1 file, group 2 file
cat > $OUT/groups.txt << 'EOF'
Afghanistan Afghanistan1.txt Afghanistan2.txt
Djibouti Djibouti1.txt Djibouti2.txt
Ethiopia NE_Ethiopia.txt SE_Ethiopia.txt
India India1.txt India2.txt
Iran Iran1.txt Iran2.txt
Kenya Kenya1.txt Kenya2.txt
Pakistan Pakistan1.txt Pakistan2.txt
Saudi_Arabia SaudiArabia1.txt SaudiArabia2.txt
Somalia Somalia1.txt Somalia2.txt
Sudan Sudan1.txt Sudan2.txt
Yemen Yemen1.txt Yemen2.txt
EOF

# Sanity check: files exist, sample counts, and names found in the VCF
echo -e "country\tfile\tn_samples\tn_in_vcf" > $OUT/group_check.tsv
while read -r c g1 g2; do
  for f in $g1 $g2; do
    if [ ! -s "$f" ]; then
      echo -e "$c\t$f\tMISSING\tMISSING" >> $OUT/group_check.tsv
      continue
    fi
    sed -i 's/\r$//' "$f"
    n=$(tr -d '\r' < $f | grep -c .)
    v=$(tr -d '\r' < $f | grep . | grep -cFxf $SAMPLE_NAMES)
    echo -e "$c\t$f\t$n\t$v" >> $OUT/group_check.tsv
  done
done < $OUT/groups.txt
column -t $OUT/group_check.tsv
awk -F'\t' -v min=$MIN 'NR>1 && ($3=="MISSING" || $3<min || $3!=$4) {
  print "WARNING: check " $2 " (n=" $3 ", in VCF=" $4 ")" }' $OUT/group_check.tsv

# One vcftools run per country (windowed output kept for plotting)
xargs -P $THREADS -L 1 bash -c '
  apptainer exec $SIF vcftools --gzvcf $VCF \
    --weir-fst-pop $1 --weir-fst-pop $2 \
    --fst-window-size 100000 --fst-window-step 50000 \
    --out $OUT/$0 > $OUT/$0.err 2>&1
' < $OUT/groups.txt

# Summary table (genome-wide mean and weighted Fst from each log)
echo -e "country\tn_group1\tn_group2\tmean_Fst\tweighted_Fst" > $OUT/within_country_fst.tsv
while read -r c g1 g2; do
  n1=$(tr -d '\r' < $g1 | grep -c .)
  n2=$(tr -d '\r' < $g2 | grep -c .)
  m=$(grep "mean Fst estimate" $OUT/$c.log 2>/dev/null | awk '{print $NF}')
  w=$(grep "weighted Fst estimate" $OUT/$c.log 2>/dev/null | awk '{print $NF}')
  echo -e "$c\t$n1\t$n2\t${m:-NA}\t${w:-NA}" >> $OUT/within_country_fst.tsv
done < $OUT/groups.txt

column -t $OUT/within_country_fst.tsv
echo "Done."