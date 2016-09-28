# get FASTA for peaks using getDNAClient

require( getDNAClient)

tbl = read.table("Smad3.peaks.merge200.reproducible.fdr01.10reads.summits_plus_minus_250bp.bed") 
colnames(tbl) <- c("chrom", "chromStart", "chromEnd")
locStrings <- sprintf("%s:%d-%d", tbl$chrom, tbl$chromStart, tbl$chromEnd)
gdc <- getDNAClient("mm9")

sink("Smad3.peaks.merge200.reproducible.fdr01.10reads.summits_plus_minus_250bp.fa")
for( i in 1:length(locStrings) ) {
  cat( ">",locStrings[i],"\n" )
  seq <- getSequenceByLocString(gdc, locStrings[i] )
  cat( seq , "\n" )
}
sink()




