# find overlap between peaks from MACS

require(GenomicRanges)
require(edgeR)

# load peaks called in each sample

peakfiles = Sys.glob("*p1e3_peaks.narrowPeak.fdr01.narrowPeak")

# ENCODE blacklist
blacklist = read.table( "../mm9-blacklist.bed.gz")
colnames(blacklist) = c("chr","start","end")
blacklist = makeGRangesFromDataFrame( blacklist )

peaks = list()
for( i in peakfiles ) {
  tmp = read.table(i)
  colnames(tmp) = c("chr","start","end","name","score","strand","score2","pval","fdr","width")
  tmp = tmp[ tmp$fdr > 2 , ]
  gr = makeGRangesFromDataFrame( tmp , keep.extra.columns = T )
  peaks[[i]] = gr
}

# load the positions of merged peaks (bedtools merge)

merged.peaks = read.table("MACS.fdr01.merged.200bp.bed")
colnames(merged.peaks) = c("chr","start","end")
loc = paste( merged.peaks[,1] , ":" , merged.peaks[,2] , "-" , merged.peaks[,3] , sep="")
merged.peaks = makeGRangesFromDataFrame( merged.peaks )
is.blacklist = countOverlaps( merged.peaks , blacklist )
merged.peaks = merged.peaks[ is.blacklist == 0 ]
loc = loc[ is.blacklist == 0 ]

# assemble a table with the -log10(p-value) of each peak
# in each sample

peak.scores = matrix( 0 , nrow=length(merged.peaks) , ncol = 4 )
colnames(peak.scores) = gsub( "\\.(.*)","",names(peaks)[1:4] )
rownames(peak.scores) = loc
for( i in 1:4 ) {
  x = peaks[[i]]
  ranked = x[ order(x$pval , decreasing = T ) , ]
  p = as.data.frame( ranked )$pval
  o = as.data.frame(findOverlaps( ranked , merged.peaks ))
  peak.scores[ o$subjectHits , i ] = p[ o$queryHits ]
}
peak.scores = peak.scores[ rowSums( peak.scores ) != 0 , ]





# select reproducible peaks, with evidence of binding in >1 sample
# consider more stringent p-value thresholds

thresh.logp = 0
replicable.peaks = peak.scores[ rowSums( peak.scores > thresh.logp ) > 1 , ]
reppeaks.gr = merged.peaks[ rowSums( peak.scores > thresh.logp ) > 1 ]

# a very simple method for comparing the groups:
# how many peaks are detected in WT vs. Q111 mice?

# all peaks
wt = rowSums( peak.scores[,c(2,4)] > 0 ) > 0
q111 = rowSums( peak.scores[,c(1,3)] > 0 ) > 0
vennCounts( cbind( wt , q111 ))
  wt q111 Counts
1  0    0      0
2  0    1 192524
3  1    0 352885
4  1    1 102781

#replicable peaks
wt = rowSums( replicable.peaks[,c(2,4)] > 0 ) > 0
q111 = rowSums( replicable.peaks[,c(1,3)] > 0 ) > 0
colSums( cbind( wt , q111 ))
#    wt   q111 
#123449 108233 


vennCounts( cbind( wt , q111 ))

  wt q111 Counts
1  0    0      0
2  0    1   5452
3  1    0  20668
4  1    1 102781

pdf("Smad3.peakcalls.vennDiagram.pdf")
vennDiagram( cbind( wt , q111 ))
dev.off()


# a more quantitative measure
# compare the read counts in WT vs. Q111 mice with edgeR

# read in the counts from the (sorted) pairs files
readfiles = Sys.glob("/proj/price1/jpearl/ChIPseq/*/*bed.sorted")
samples = gsub("(.*)\\/","",readfiles)
samples = gsub("\\.(.*)","",samples)
Smad3 = grep("Smad3",readfiles)

n = length(Smad3)
libsize = rep(NA,n)
counts = matrix( 0 , nrow=length(reppeaks.gr) , ncol=length(Smad3) )
for( i in 1:n ) {
  reads = read.table(readfiles[ Smad3[i] ])
  colnames(reads) = c("chr","start","end")
  reads.gr = makeGRangesFromDataFrame( reads )
  overlaps = countOverlaps( reppeaks.gr , reads.gr )
  counts[,i] = overlaps
  libsize[i] = nrow(reads)
}
colnames(counts) = samples[ Smad3 ]
rownames(counts) = rownames(replicable.peaks)

# filter out peaks with very few counts
# consider RPKM and raw count thresholds
# use a threshold that half the samples must have at least 10 reads
# this is a minimal threshold at which edgeR is likely to provide 
# meaningful results

peak.lengths = width(ranges(reppeaks.gr))
rpkm = rpkm( counts , gene.length = peak.lengths )
rpkm.median = apply( rpkm , 1 , median )
rpkm.min = apply( rpkm , 1 , min )

counts.median = apply( counts , 1 , function(x) sort(x)[2] )

reppeaks2.gr = reppeaks.gr[ counts.median > 10 ]
counts2 = counts[ counts.median > 10 , ]
peaks2 = peak.scores[ rownames(counts2) , ]
wt = rowSums( peaks2[,c(2,4)] > 0 ) > 0
q111 = rowSums( peaks2[,c(1,3)] > 0 ) > 0
vennCounts( cbind( wt , q111 ))
#  wt q111 Counts
#1  0    0      0
#2  0    1    632
#3  1    0   5419
#4  1    1  51721

pdf("venn.pdf")
vennDiagram( cbind( wt , q111 ) )
dev.off()



# perform a statistical analysis with edgeR

# genotype
group = factor(c("Q111","WT","Q111","WT"))
group = relevel( group , "WT" )
# technician (batch)
person = c("DB","DB","JP","JP")

# create the design matrix
design = model.matrix( ~ 0 + group + person )
colnames(design) = gsub("group","",colnames(design))

# perform edgeR analysis
y = DGEList( counts2 , lib.size = libsize )
y = estimateDisp( y )
contrasts = makeContrasts( genotype = Q111 - WT , levels = design )
fit = glmFit( y , design )
fit = glmLRT( fit , contrast = contrasts[,1] )
p.edgeR = as.data.frame( topTags(fit,n=Inf,sort.by="none"))

de = p.edgeR[ p.edgeR$PValue < 0.01 , ]
table( de$logFC < 0 )
#FALSE  TRUE 
#    5   133 


p.edgeR[ order( p.edgeR$PValue )[1:50] , ]

                              logFC       logCPM        LR       PValue
chrX:23585555-23585849    -2.555338  0.556653922 25.425439 4.598144e-07
chr3:96303242-96304360    -2.681291  2.570303598 22.179003 2.483753e-06
chr1:115800576-115800855  -2.283316  0.311871227 21.030968 4.519196e-06
chr17:23683445-23684891   -2.554182  2.592002649 20.563464 5.768664e-06
chr3:96304635-96305439    -1.772436  0.805404268 13.122074 2.918368e-04
chr11:115345639-115346261 -1.755002 -0.480512634 12.391531 4.312856e-04
chr17:7145609-7146855     -1.694536  1.044502110 12.294851 4.542098e-04
chr11:119336274-119337001 -1.722566 -0.597709420 11.783862 5.974640e-04
chr18:85263756-85264930   -1.728276 -0.398873983 11.688296 6.289447e-04
chr17:46579167-46580249   -1.666467 -0.522646232 11.255413 7.939116e-04
chr11:3381329-3382716     -1.609877  0.150272464 11.199302 8.182812e-04
chr12:73096715-73098352   -1.571274  0.279644314 10.903543 9.598051e-04
chr6:3150807-3151528      -1.711999  5.507119017 10.847639 9.892179e-04

# note: only a few of these peaks reach a significant FDR threshold
# not surprising given the small sample size


# create a results table with all the information and stats
# write to disk

a = counts2
colnames(a) = paste("counts",colnames(a),sep=".")
b = peaks2
colnames(b) = paste("MACS.logp",colnames(b),sep=".")
res = data.frame( a , b , p.edgeR )

write.csv( res , file="Smad3.MACS.edgeR.results.2016-7-27_FINAL.csv" )

x = res
loc = matrix( unlist(strsplit(rownames(x),split=":")) , ncol = 2 , byrow = T )
chr = loc[,1]
start = matrix( unlist(strsplit(loc[,2],split="-")) , ncol = 2 , byrow = T )[,1]
end = matrix( unlist(strsplit(loc[,2],split="-")) , ncol = 2 , byrow = T )[,2]
bed = data.frame( chr , start , end )
write.table( bed , sep="\t" , row.names=F , col.names=F , quote=F ,
   file="Smad3_peaks.merge-d200.reproducible_fdr01+10reads_FINAL.bed" )









