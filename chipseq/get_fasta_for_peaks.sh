# get fasta sequences of the peak summits

cd /proj/price1/sament/chipseq/Smad3

bedtools intersect \
-a Smad3.combined.MACS.p1e3_peaks.narrowPeak.fdr01.narrowPeak \
-b Smad3_peaks.merge-d200.reproducible_fdr01+10reads_FINAL.bed \
> Smad3.peaks.merge200.reproducible.fdr01.10reads.narrowPeak

awk 'BEGIN{ OFS="\t"} { midPos=$2+$10; print $1, midPos, midPos+1 }' \
  Smad3.peaks.merge200.reproducible.fdr01.10reads.narrowPeak \
  > Smad3.peaks.merge200.reproducible.fdr01.10reads.summits


awk 'BEGIN{ OFS="\t"}{ print $1, $2-250, $2+250 }' \
Smad3.peaks.merge200.reproducible.fdr01.10reads.summits \
> Smad3.peaks.merge200.reproducible.fdr01.10reads.summits_plus_minus_250bp.bed

Rscript ../scripts/get_fasta_for_peaks.R

