
library( irr )

options(stringsAsFactors=F)
smad3.enrich = read.csv("GO_enrichments.top_Smad3_targets.2016-8-1.csv")



sig = smad3.enrich$set[ smad3.enrich$q < 0.01 ]
setGenes = strsplit( smad3.enrich$setGenes , split = ";" )
names(setGenes) = smad3.enrich$set
setGenes = setGenes[ sig ]

 
#load("/proj/price1/sament/resources/goterms_mmu_2015-11-4.RData")
#n = length(genesets)
#size = sapply( 1:n , function(x) length( genesets[[x]] ) )
#sets = genesets[ size <= 500 & size >= 10 ]



genes = read.table("/proj/price1/sament/resources/refGenomes/mm9/refGene.txt.gz" )
genes = unique( genes[,13] )
genes = as.character( genes )

n = length(sig)

k = matrix( NA , n , n )
for( i in 1:n ) {
  for( j in 1:n ) {
     a = genes %in% setGenes[[ i ]]
     b = genes %in% setGenes[[ j ]]
     k[i,j] = kappa2( cbind( a , b ))$value
  }
  cat( i , "\n" )
}
rownames(k) = colnames(k) = sig

setClust = hclust( as.dist( 1-k ) , method = "average" )

pdf("go_term_clusters.pdf")
plot( setClust )
dev.off()


clust = cutree( setClust , h = 0.7 )
names(clust) = sig
#sort( clust )

rownames(smad3.enrich) = smad3.enrich$set
sig.enrich = smad3.enrich[ sig , ]
sig.enrich$clust = clust
sig.enrich = sig.enrich[ order( sig.enrich$p ) , ]
sig.enrich = sig.enrich[ duplicated(sig.enrich$clust) == F , ]
sig.enrich = sig.enrich[ is.na( sig.enrich[,1] ) == F , ]
write.csv( sig.enrich , file="GO_enrichment.topterms.nonoverlapping.2016-08-13.csv" )



sig.enrich[ , c("set","p") ]


actin_filament-based_process                  4.208374e-11
mRNA_processing                               4.159394e-09
actin_binding                                 2.062277e-08
neuromuscular_process_controlling_balance     1.224662e-07
histone_modification                          1.692694e-07
brain_development                             1.267140e-06
chromatin_binding                             2.713612e-06
actin_filament-based_movement                 2.822851e-06
regulation_of_cell_projection_organization    4.455988e-06
lamellipodium                                 5.454795e-06
protein_serine/threonine_kinase_activity      9.895540e-06
protein_deacetylation                         1.135973e-05
centrosome                                    1.368243e-05
purine_ribonucleotide_catabolic_process       1.622276e-05
kinase_binding                                2.313849e-05
phosphoric_ester_hydrolase_activity           2.514136e-05
neuronal_cell_body                            2.532075e-05
protein_kinase_binding                        2.764772e-05
cellular_protein_catabolic_process            3.023238e-05
transcriptional_repressor_complex             3.718618e-05
respiratory_system_development                5.473966e-05
kinesin_binding                               1.116028e-04
endocytosis                                   1.121766e-04
negative_regulation_of_ERBB_signaling_pathway 1.361494e-04




