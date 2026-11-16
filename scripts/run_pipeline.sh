#!/usr/bin/env bash
# Genotype, annotate and plot pathogenic tandem repeat expansions in PacBio PureTarget data.
#
# Expected input (run from the project folder):
#   bam/<sample>.bam (+ .bai)   aligned HiFi reads, GRCh38
#   samples.tsv                 header line, then: sample<TAB>karyotype (XX or XY)
#   ref/hg38.fa                 GRCh38 no-alt analysis set
#
# Tools: trgt 5.1.0, stranger 0.10.2, samtools 1.25, bcftools 1.25 (see environment.yml)
set -euo pipefail

REPEATS_URL=https://raw.githubusercontent.com/dashnowlab/STRchive/main/data/catalogs/STRchive-disease-loci.hg38.TRGT.bed
STRANGER_URL=https://raw.githubusercontent.com/Clinical-Genomics/stranger/v0.10.2/stranger/resources/variant_catalog_hg38.json

mkdir -p results/plots
SAMPLES=$(tail -n +2 samples.tsv | cut -f1)

# --- 0. Reference index and catalogues ---------------------------------------
[ -s ref/hg38.fa.fai ]     || samtools faidx ref/hg38.fa
[ -s repeats.bed ]         || wget -q -O repeats.bed "$REPEATS_URL"
[ -s stranger_hg38.json ]  || wget -q -O stranger_hg38.json "$STRANGER_URL"
echo "Repeat loci in catalogue: $(grep -vc '^#' repeats.bed)"

# --- 1. Genotype with TRGT (PureTarget preset, sex-aware) --------------------
tail -n +2 samples.tsv | while read -r SAMPLE KARYO; do
  echo "== trgt genotype $SAMPLE ($KARYO)"
  trgt genotype \
    --preset targeted \
    --genome ref/hg38.fa \
    --repeats repeats.bed \
    --reads "bam/${SAMPLE}.bam" \
    --karyotype "$KARYO" \
    --output-prefix "results/${SAMPLE}"
done 2> results/trgt_genotype.log

# --- 2. Sort and index (needed by trgt plot) ---------------------------------
for S in $SAMPLES; do
  bcftools sort -Oz -o "results/${S}.sorted.vcf.gz" "results/${S}.vcf.gz"
  bcftools index -f "results/${S}.sorted.vcf.gz"
  samtools sort -o "results/${S}.spanning.sorted.bam" "results/${S}.spanning.bam"
  samtools index "results/${S}.spanning.sorted.bam"
done

# --- 3. Annotate disease status with Stranger --------------------------------
for S in $SAMPLES; do
  stranger --trgt -f stranger_hg38.json "results/${S}.sorted.vcf.gz" \
    > "results/${S}.stranger.vcf" 2> "results/${S}.stranger.log"
  bgzip -f "results/${S}.stranger.vcf"
  bcftools index -f "results/${S}.stranger.vcf.gz"
done

# --- 4. Table of non-normal alleles ------------------------------------------
OUT=results/expanded_alleles.tsv
printf "sample\tlocus\tstatus\tdisease\tGT\tmotif_count\tsupporting_reads\tallele_length_bp\tnormal_max\tpathologic_min\n" > "$OUT"
for S in $SAMPLES; do
  bcftools query \
    -i 'INFO/STR_STATUS="pre_mutation" || INFO/STR_STATUS="full_mutation"' \
    -f "${S}\t%INFO/TRID\t%INFO/STR_STATUS\t%INFO/Disease\t[%GT\t%MC\t%SD\t%AL]\t%INFO/STR_NORMAL_MAX\t%INFO/STR_PATHOLOGIC_MIN\n" \
    "results/${S}.stranger.vcf.gz"
done >> "$OUT"
column -t -s $'\t' "$OUT"

# --- 5. Waterfall plots -------------------------------------------------------
plot() {  # plot SAMPLE REPEAT_ID motifs|meth
  trgt plot \
    --genome ref/hg38.fa \
    --repeats repeats.bed \
    --vcf "results/$1.sorted.vcf.gz" \
    --spanning-reads "results/$1.spanning.sorted.bam" \
    --repeat-id "$2" \
    --plot-type waterfall \
    --show "$3" \
    --image "results/plots/$1_$2_$3.svg"
}
plot sample03 FXS_FMR1   motifs
plot sample07 DM1_DMPK   motifs
plot sample10 FRDA_FXN   motifs
plot sample03 FXS_FMR1   meth
plot sample08 FXS_FMR1   meth
plot sample11 SCA1_ATXN1 motifs
plot sample06 HD_HTT     motifs
echo "Done. Plots in results/plots/"
