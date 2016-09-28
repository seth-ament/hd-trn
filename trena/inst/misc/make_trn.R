# generate a TRN model from the striatum allelic series data
options( stringsAsFactors=F )

# load the expression data
setwd("/proj/price1/CHDI/users/sament/allelic_series_trn")
load("combined_striatum_fpkm_and_metadata.RData")

# load the TFBS model
library( trena )
data( "tfbs.mmu.mgi_symbol" )

# TRN reconstruction
rownames( fpkm ) = anno[,1]
low.expr = which( rowSums( fpkm == 0 ) > ( ncol( fpkm ) / 2 ) )
expr = as.matrix(fpkm[ -low.expr , ])
anno = anno[ -low.expr , ]
testTRN = estLambda( expr , tfbs.mmu.mgi_symbol , n = 100 )
trn0 = fitTRN( expr , tfbs.mmu.mgi_symbol , lambda = 0.1 )
save( expr , anno , trn0 , meta ,  file = "trn_2015-6-4.RData" )

# plot summary statistics
pdf("summary_plots_2015-6-4.pdf")
summary_plots( trn0 )
dev.off()
