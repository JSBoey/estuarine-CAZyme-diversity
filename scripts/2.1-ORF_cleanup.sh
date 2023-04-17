#!/bin/bash -e

# Concatenate all ORF predictions and remove metadata from ORF predicted sequences.
# Metadata is also consolidated from the GFF files.

# Directories
DIR=data/2.orf_prediction

# Concatenate predictions, remove metadata and asterisks
# One can obtain statistics of whether it was a partial gene through the GFF3 annotations.
cat $DIR/*.faa \
  | cut -f 1 -d ' ' \
  | sed -e 's/\*//g' \
  > $DIR/allbins_pred.faa

# Consolidate GFF data
printf "bin\tnode\tsource\ttype\tstart\tend\tgff_score\tstrand\tphase\tseqid\tpartial\tstart_type\trbs_motif\trbs_spacer\tgc_cont\tconf\tscore\tcscore\tsscore\trscore\tuscore\ttscore\n" > $DIR/allbins_pred.metadata.tsv

for i in $DIR/*.gff; do
  bin=$(basename $i .gff)
  grep -v '#' $i \
    | sed -e 's/;$//' \
    | sed -e 's/;/\t/g' \
    | sed -E 's/\w+=//g' \
    | awk '{FS="\t"; OFS="\t"} {split($9, a, "_"); $1=$1"_"a[2]; print}' \
    | awk -v mag="$bin" '{FS="\t"; OFS="\t"} {print mag, $0}' \
    >> $DIR/allbins_pred.metadata.tsv
done
