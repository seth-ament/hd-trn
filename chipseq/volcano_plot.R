# create a volcano plot with Smad3 edgeR fold changes and p-values

p.edgeR = read.csv("Smad3.MACS.edgeR.results.2016-6-29.csv")

logp = -log10( p.edgeR$PValue )

logfc = p.edgeR$logFC
logfc.jitter = jitter( logfc )

png("volcano.png")
par( bty = "l" , lwd = 2 , mar = c(5,5,2,1) )
plot( 
   x = logfc.jitter , y = logp , 
   xlim = c(-3,3) , ylim = c(0,7) ,
   ylab = "" , xlab = "" , cex.axis = 2 )
# abline( v = 0 )
abline( h = 2 , lty = 1 , lwd = 2 )
abline( v = 0 , lty = 1 , lwd = 2 )
#abline( v = 1 , lty = 2 )
text( x = -1.7 , y = 7 , adj = 0.5 , cex = 3 , xpd = NA ,
   length( which( logfc < 0 & logp > 2 )) )
text( x = 1.7 , y = 7 , adj = 0.5 , cex = 3 , xpd = NA ,
   length( which( logfc > 0 & logp > 2 )) )
mtext( side = 1 , cex = 2 , line = 3.5 , "log2( Fold Change )" )
mtext( side = 2 , cex = 2 , line = 3.5 , "-log10( p-value )" )
dev.off()





