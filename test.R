library(gplots)
library(RColorBrewer)
pdf("test.pdf", width=10, height=10)
x <- as.matrix(read.table("test.mat", header=T, row.names=1))
fc = nrow(x) / 2
heatmap.2(x, symm=TRUE, trace="none", symbreaks=TRUE, density.info="none",
	  col = rev(colorRampPalette(brewer.pal(11, "RdBu"))(100)),
	  lmat=rbind(c(4, 4, 4), c(2,1,0), c(0,3,0)), lwid=c(1, 6, 1), lhei=c(2, 10, 1))
heatmap.2(x[1:fc,1:fc], symm=TRUE, col = cm.colors(256), trace="none", symbreaks=TRUE,
             lmat=rbind(c(4, 4, 4), c(2,1,0), c(0,3,0)), lwid=c(1, 6, 1), lhei=c(2, 10, 1))
heatmap.2(x[(fc+1):nrow(x),1:fc], symm=TRUE, col = cm.colors(256), trace="none", symbreaks=TRUE,
             lmat=rbind(c(4, 4, 4), c(2,1,0), c(0,3,0)), lwid=c(1, 6, 1), lhei=c(2, 10, 1))
heatmap.2(x[(fc+1):nrow(x), (fc+1):nrow(x)], symm=TRUE, col = cm.colors(256), trace="none", symbreaks=TRUE,
             lmat=rbind(c(4, 4, 4), c(2,1,0), c(0,3,0)), lwid=c(1, 6, 1), lhei=c(2, 10, 1))
dev.off()
