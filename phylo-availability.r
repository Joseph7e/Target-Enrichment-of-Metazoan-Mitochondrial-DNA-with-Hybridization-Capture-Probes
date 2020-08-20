install.packages("wesanderson")
install.packages("ape")
library(ape)
library(wesanderson)

tre = read.tree( "C:/Users/joseph7e/Documents/mito-enrichment-data/dunn_tree.tre" )
x = read.table( "C:/Users/joseph7e/Documents/mito-enrichment-data/sampling-mitos.csv", header=TRUE, sep=",", stringsAsFactors=FALSE )
rownames(x) <- x$clade

x = x[ match(tre$tip.label, x$clade), ]

pal <- wes_palette("Zissou1", 5, type = "discrete")

pdf( "Figure1.pdf",width=7.5,height=8 )
plot( tre, font=2, edge.width=2, cex=0.8, x.lim=c(-10,80), label.offset=1, node.depth=1 )
# plot( tre, font=2, edge.width=2, cex=0.8, x.lim=c(-10,40), label.offset=1, node.depth=2 )

internal_color = "darkgray"
nodelabels( tre$node.label, frame="none", cex=0.8, font=2, adj=c(1.1, -.7), col=internal_color )
nodelabels(which(nchar(tre$node.label)>0)+length(tre$tip.label), pch=20, col=internal_color, cex = 1, frame="none", text="")

# Drawing on http://treethinkers.org/tutorials/dataset-visualization/
x_offset = 55 # x position of start of bars
one_size = 0.2 # how big the bar should be for a value of 1
scale = 0.7 # scales overall length of bars

abline(v=( x_offset ), col="gray")
#abline(v=( x_offset +(log(1) + one_size)/scale), col="gray")
abline(v=( x_offset +(log(10) + one_size)/scale), col="gray")
abline(v=( x_offset +(log(100) + one_size)/scale), col="gray")
abline(v=( x_offset +(log(1000) + one_size)/scale), col="gray")
abline(v=( x_offset +(log(10000) + one_size)/scale), col="gray")
abline(v=( x_offset +(log(100000) + one_size)/scale), col="gray")
abline(v=( x_offset +(log(1000000) + one_size)/scale), col="gray")
# abline(v=( x_offset +(log(10000000) + one_size)/scale), col="gray")

axis_pos = 34.4
text_cex = 0.4
text(( x_offset ), axis_pos, "0", srt=270, cex=text_cex, adj=c(1, NA))
#text(( x_offset +(log(1) + one_size/scale)), axis_pos, "1", srt=270, cex=text_cex, adj=c(1, NA))
text(( x_offset +(log(10) + one_size)/scale), axis_pos, "10", srt=270, cex=text_cex, adj=c(1, NA))
text(( x_offset +(log(100) + one_size)/scale), axis_pos, "100", srt=270, cex=text_cex, adj=c(1, NA))
text(( x_offset +(log(1000) + one_size)/scale), axis_pos, "1k", srt=270, cex=text_cex, adj=c(1, NA))
text(( x_offset +(log(10000) + one_size)/scale), axis_pos, "10k", srt=270, cex=text_cex, adj=c(1, NA))
text(( x_offset +(log(100000) + one_size)/scale), axis_pos, "100k", srt=270, cex=text_cex, adj=c(1, NA))
text(( x_offset +(log(1000000) + one_size)/scale), axis_pos, "1m", srt=270, cex=text_cex, adj=c(1, NA))


# segments( x_offset, 1:nrow(x), x_offset + pmax( (log(x[,5]) + one_size)/scale, 0 ), 1:nrow(x), lwd=13, col=pal[1], lend=1 ) # species predicted
segments( x_offset, 1:nrow(x), x_offset + pmax( (log(x[,2]) + one_size)/scale, 0 ), 1:nrow(x), lwd=11, col=pal[2], lend=1 ) # species described
segments( x_offset, 1:nrow(x), x_offset + pmax( (log(x[,3]) + one_size)/scale, 0 ), 1:nrow(x), lwd=9, col=pal[4], lend=1 ) # species with molecular data
segments( x_offset, 1:nrow(x), x_offset + pmax( (log(x[,4]) + one_size)/scale, 0 ), 1:nrow(x), lwd=7, col=pal[5], lend=1 ) # species with sequenced genomes

# legend( "topleft", c("Species predicted","Species described", "Species with molecular data", "Species with sequenced genomes"), col=c(pal[1],pal[2], pal[4], pal[5]), pch = 15, cex=0.7 )
#legend( "topleft", c("Species described", "Species with molecular data", "Species with sequenced genomes"), col=c(pal[2], pal[4], pal[5]), pch = 15, cex=0.7 )
legend( "topleft", c("Species described", "Species with sequenced mitochondrial genomes"), col=c(pal[2], pal[5]), pch = 15, cex=0.7 )


dev.off()
getwd()

genomes = sum(x$genomes)
vertarthro = ( x["Arthropoda","genomes"] + x["Craniata","genomes"] ) / genomes

cat("Numbers in manuscript:")
cat("Number of published animal genomes: " , genomes, fill=TRUE)
cat("Percent of published animal genomes that are from vertebrates and arthropods: " , vertarthro*100, fill=TRUE)
cat("Percent of described animal species with published genomes: " , sum(x$genomes)/ sum(x$species)*100, fill=TRUE)
cat("Number of described animal speciess: " , sum(x$species), fill=TRUE)
cat("Number of animal species with molecular data: " , sum(x$ncbi), fill=TRUE)
cat("Percent of described animal species with molecular data: " , sum(x$ncbi)/sum(x$species)*100, fill=TRUE)
cat("Percent of projected animal species with molecular data: " , sum(x$ncbi)/sum(x$projected, na.rm=TRUE)*100, fill=TRUE)
