## Collection of Heatmap from ggplot2, gplots
## Color
myColor <- function(heatmap, common, palette){
  library(RColorBrewer)
  redgreen(75)
  colorRampPalette(brewer.pal(9,"Blues"))(100)
  display.brewer.all()
  cr <- colorRampPalette(colors = c("#2927FF","#FFFFFF","#DF5C5C"), bias=3)
}
# Color function to generate green-red heat maps
my.colorFct <- function(n = 50, low.col = 0.45, high.col=1, saturation = 1) { 
	if (n < 2) stop("n must be greater than 2")
	n1 <- n%/%2
	n2 <- n - n1
	c(hsv(low.col, saturation, seq(1,0,length=n1)), hsv(high.col, saturation, seq(0,1,length=n2))) 
}

## common image heatmap, gplots and heatmap.2 with key and without key
myHeatmap <- function(type, data, colors, key, ..){
 require(reshape)
 library(gplots)
 library(ggplot2)
 library(scales)
 heatmap.2(exprs(esetSel), col=redgreen(75), scale="row", ColSideColors=patientcolors,
           key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)
 mtscaled <- as.matrix(scale(mtcars))
 ## density.info = "none"
 ## bias for key color distribution
 ## las for xlab, ylab font angle,
 ## subtitle
 # create heatmap and don't reorder columns
 # use breaks to ajust key values distribution
 heatmap(mtscaled, Colv=F, scale='none')
# cluster rows
 hc.rows <- hclust(dist(mtscaled))
 plot(hc.rows)
# transpose the matrix and cluster columns
 hc.cols <- hclust(dist(t(mtscaled)))
# for all
 heatmap(mtscaled, Colv=as.dendrogram(hc.cols), scale='none')
# draw heatmap for first cluster
 heatmap(mtscaled[cutree(hc.rows,k=2)==1,], Colv=as.dendrogram(hc.cols), scale='none')
# draw heatmap for second cluster
 heatmap(mtscaled[cutree(hc.rows,k=2)==2,], Colv=as.dendrogram(hc.cols), scale='none')
palette <- colorRampPalette(c('#f0f3ff','#0033BB'))(256)
heatmap(mtscaled, Colv=F, scale='none', col=palette)

nba <- read.csv("http://datasets.flowingdata.com/ppg2008.csv", sep=",")
nba <- nba[order(nba$PTS),]
row.names(nba) <- nba$Name
nba <- nba[,2:20]
nba_matrix <- data.matrix(nba)
nba_heatmap <- heatmap(nba_matrix, Rowv=NA, Colv=NA, col = cm.colors(256), scale="column", margins=c(5,10))
nba_heatmap <- heatmap(nba_matrix, Rowv=NA, Colv=NA, col = heat.colors(256), scale="column", margins=c(5,10))
## illustrator and inkscape

## data <- read.csv("http://dl.dropbox.com/u/2505196/fruit.txt", head=TRUE, sep=",")
data.m = melt(data)
data.m <- ddply(data.m, .(variable), transform, rescale = rescale(value))
p <- ggplot(data.m, aes(variable, people)) + geom_tile(aes(fill = rescale), 
                                                   colour =   "white") 
     p + scale_fill_gradient(low = "white", high = "steelblue")

 ## df <- data.frame(expand.grid(x = 1:4, y = 1:4), v = runif(16, -10, 10))
## ggplot(df, aes(x, y, fill = v, label = sprintf("%.1f", v))) +
      ##   geom_tile() + geom_text() +
      ##   scale_fill_gradient2(low = "blue", high = "red")

## data_matrix <- cbind(c(rnorm(30,-2.5,sd= 0.85),rnorm(30,25,sd= 8),rnorm(30,6,sd= 3)),
##             c(rnorm(30,-2.5,sd= 0.85),rnorm(30,25,sd= 8),rnorm(30,6,sd= 3)),
##             c(rnorm(30,-2.5,sd= 0.85),rnorm(30,25,sd= 8),rnorm(30,6,sd= 3)),
##             c(rnorm(30,-2.5,sd= 0.85),rnorm(30,25,sd= 8),rnorm(30,6,sd= 3)))
## breaks = c(5/6* -5, 5/6* -4, 5/6* -3, 5/6* -2, 5/6* -1, 0 ,5/6* 1, 50/6*1, 50/6*2, 50/6*3, 50/6*4, 50/6*5)
## hv <- heatmap.2(data_matrix, 
##     scale="none",
##     Rowv=NA,
##     Colv=NA,
##     col = rev(brewer.pal(11,"RdBu")),
##     margins=c(5,5),
##     cexRow=0.5, cexCol=1.0, 
##     breaks=breaks,
##     ylab= "Mutations",
##     main = "heatmap",
##     key=TRUE,keysize=1.5, trace="none")

# hclust

## heatmap.2 parameters suggestions

## breaks <- seq(mn, mx, (mx-mn)/(n))
## heatmap.2(em_t,col = cr(n),key=TRUE,Colv=NA, trace="none",
##           breaks=breaks, notecex=0.7, keysize=0.7, density.info="histogram", cexRow=1.00,cexCol=1.00, margins=c(17.0,17.0),lmat=rbind( c(0, 3, 4), c(2,1,1 ) ), lwid=c(1.5, 4, 2 ))

      ## heatmap.2(as.matrix(your data),col =
      ## colorRampPalette(c("white","green","green4","violet","purple"))(100))
      ##       hmcol<-rev(brewer.pal(11,"RdBu"))
      ## heatmap(data,col=hmcol)
      ## You can also use something like

      ## hmcols<-colorRampPalette(c("red","white","blue"))(256)
      ## breaks = breaks, ylab="xxx", main="heatmap", key=T,keysize=1.5, Rowv=NA, Colv=NA,cexRow=0.5, cexCol=1.0
      heatmap.2(as.matrix(exp.match[as.numeric(w.index),]), col=brewer.pal(9,"Blues"), main= paste(gene, "wilcoxon cutoff p.value", cutoff),
                trace='none', notecex=0.2, ColSideColors = patient.category, dendrogram = "both", cexRow=0.4, scale="row")
      # t.test
      t.index <- which(t.p.value < cutoff)
      head(exp.match[which(t.p.value<cutoff), ])
      ## Col colors are patients
      heatmap.2(as.matrix(exp.match[as.numeric(t.index),]), col=brewer.pal(9,"Blues"),main= paste(gene, "t.test cutoff p.value", cutoff),
                trace='none', notecex=0.2, ColSideColors = patient.category, dendrogram = "both", cexRow=0.4, scale="row")
}

## heatmap contour
myContour <- function(data, ..){
d <- structure(list(X = c(-60L, 60L, 7L, -22L, 59L, 29L, -58L, 60L, 
7L, -21L, 61L, 29L, -57L, -22L, 59L, 29L, -56L, 61L, 8L, -20L, 
62L, 30L), Y = c(-18L, -62L, 14L, -60L, 58L, 22L, -18L, -61L, 
14L, -59L, 59L, 22L, -18L, -59L, 60L, 24L, -17L, -60L, 16L, -58L, 
60L, 23L)), .Names = c("X", "Y"), class = "data.frame", row.names = c(NA, 
-22L))
require(MASS)
dens <- kde2d(d$X, d$Y, h=75, n=50)  #overrode default bandwidth
filled.contour(dens)

require(splancs)
data(bodmin)
plot(bodmin$poly,asp=1,type="n")
image(kernel2d(as.points(bodmin),bodmin$poly, h0=2,nx=100,ny=100),
        add=TRUE, col=terrain.colors(20))
pointmap(as.points(bodmin),add=TRUE)
polymap(bodmin$poly,add=TRUE)
bodmin.xy<-coordinates(bodmin[1:2])
apply(bodmin$poly,2,range)
grd1<-GridTopology(cellcentre.offset=c(-5.2,-11.5),cellsize=c(0.2,0.2),cells.dim=c(100,100))
k100<-spkernel2d(bodmin.xy,bodmin$poly,h0=1,grd1)
k150<-spkernel2d(bodmin.xy,bodmin$poly,h0=1.5,grd1)
k200<-spkernel2d(bodmin.xy,bodmin$poly,h0=2,grd1)
k250<-spkernel2d(bodmin.xy,bodmin$poly,h0=2.5,grd1)
if(.sp_lt_0.9()){
    df<-AttributeList(list(k100=k100,k150=k150,k200=k200,k250=k250))
} else{
    df<-data.frame(k100=k100,k150=k150,k200=k200,k250=k250)
}
kernels<-SpatialGridDataFrame(grd1,data=df)
spplot(kernels,checkEmptyRC=FALSE,col.regions=terrain.colors(16),cuts=15)

d <- structure(list(X = c(-60L, 60L, 7L, -22L, 59L, 29L, -58L, 60L, 
7L, -21L, 61L, 29L, -57L, -22L, 59L, 29L, -56L, 61L, 8L, -20L, 
62L, 30L), Y = c(-18L, -62L, 14L, -60L, 58L, 22L, -18L, -61L, 
14L, -59L, 59L, 22L, -18L, -59L, 60L, 24L, -17L, -60L, 16L, -58L, 
60L, 23L)), .Names = c("X", "Y"), class = "data.frame", row.names = c(NA, 
-22L))
require(MASS)
dens <- kde2d(d$X, d$Y, h=75, n=50)  #overrode default bandwidth
filled.contour(dens)

## example
kern.obj <- structure(c(-161.913250909479, 154.013482116162, 31.6474639061300, 
17.7340639366637, -102.170823111156, 17.6809699563749, 90.505728795223, 
143.854796792441, -70.1806511117134, 230.600354761065, 133.500211485414, 
-225.74140063979, 220.599384351733, -55.5956512970632, 128.631103577179, 
-36.9382693513206, 86.1151116370548, -67.9572171234925, 138.313636950703, 
59.4122360493993, -128.418347257186, 28.4313444162254, -253.438542232118, 
-2.62936998134802, 96.6705573949275, 126.350347596454, -76.3053490233138, 
-98.1667749493097, -132.615954657406, -239.003804126569, -32.052834858324, 
152.055005227299, -171.132473363859, -96.0272921226682, -91.4859761718545, 
172.662664785850, 92.3258005260648, -9.33884441249779, -24.4260189034222, 
-171.435971200881, 84.9052731383744, -171.768339197942, -13.5871193263486, 
-51.839925496188, -193.00283491136, 57.1126055897217, -40.890549093622, 
83.600134171797, 6.66515671609591, -261.487889322599, 138.624659821426, 
158.911075756538, 111.598989561161, 62.6150728399137, -155.366548557697, 
95.9501552130317, -32.0820888905296, -85.4929337702259, -178.010310820340, 
100.526315864149, -190.431234842843, 223.959168312304, -10.693030515916, 
-155.820490522984, 87.7527496146106, 293.991051801515, -69.1568338969259, 
77.0440461941863, -137.088789092018, -284.434533670747, -52.9437134391306, 
129.855822783810, 147.208098412254, -144.394565933009, 11.1193096498363, 
-26.6883210946328, 36.3402764034715, -27.5111672678245, 161.017920279498, 
133.961438546933, -139.924061267615, -194.861248844460, -138.902485043792, 
-59.6746738747854, -193.856125217724, 58.9319665388044, -151.870347293954, 
185.500357832384, 77.8198201646078, 217.406148533358, 125.978806993972, 
-96.8970637852723, 85.2079461295587, -71.5845844358825, 90.0263897196243, 
-3.85398693321446, -233.945188963933, -252.371240484100, -152.282817449886, 
-175.448833834566, 74.8285138048232, 218.884530197829, -65.9526397939771, 
113.776709279045, -69.4176647812128, -196.919950610027, 268.779812799767, 
119.294722331688, 272.239590017125, -161.720151454210, -16.8415614869446, 
-13.6117741931230, -96.0124779492617, 157.184316962957, 188.061125110835, 
-214.437550725415, 121.667246008292, 89.747676299885, -4.44232751615345, 
-106.699166027829, -261.718519963324, -42.1719799283892, -78.4863225650042, 
204.811030067503, 265.774235548452, 38.5583057999611, -239.476124290377, 
231.875250348821, 135.243163537234, -42.7497774828225, -59.7301519475877, 
-2.99901310354471, -240.498538082466, -109.713196987286, 172.524304641411, 
113.648047484457, -221.150079695508, 131.948393024504, 62.1528406161815, 
-8.31053741276264, -76.1619768105447, 157.933613704517, -42.225355328992, 
208.729289704934, 10.0781018380076, 98.7709498498589, -74.8700814787298, 
-215.313404565677, -87.6694556325674, -139.495075587183, -28.3679623156786, 
-76.2799751479179, -138.629644783214, -164.171522296965, 16.3864661939442, 
-109.221789333969, -49.0070185158402, -23.0688956100494, 54.3438952881843, 
-145.427243504673, -18.4494345914572, 14.391646720469, -200.727640092373, 
187.278914311901, -75.3078812733293, 4.16369824670255, -191.299003595486, 
169.710802193731, -103.791763912886, 32.9403738956898, -91.6615933645517, 
-222.505887318403, 49.3231621105224, -151.363900210708, -23.9421324804425, 
-207.101033208892, 169.309269497171, -250.131661305204, 11.1456824932247, 
-193.683278560638, -66.6569401044399, -139.672750141472, -115.024601574987, 
-198.41345124878, -205.971520487219, 104.227339709178, 162.442225730047, 
-167.216443363577, -100.033209286630, 152.823372976854, -191.260906308889, 
-234.539421927184, 213.049413822591, 130.761165590957, -234.716210095212, 
6.07512393034995, -49.286244995892, -56.5862323623151, -50.971424812451, 
-168.812829069793), .Dim = c(100L, 2L), .Dimnames = list(NULL, 
    c("x", "y")))

circpol <- structure(c(37.674311717588, 75.1999401385825, 112.428788751435, 
149.213932298913, 185.410196624968, 220.874731610807, 255.467574939044, 
289.052204461029, 321.496076987398, 352.671151375484, 382.454393849214, 
410.728263557213, 437.381176452847, 462.307945665474, 485.410196624968, 
506.596755301209, 525.784008026318, 542.896231479612, 557.865891532951, 
570.633909777092, 581.149896677179, 589.372350437213, 595.268820788687, 
598.816037056963, 600, 598.816037056963, 595.268820788687, 589.372350437213, 
581.149896677179, 570.633909777092, 557.865891532951, 542.896231479612, 
525.784008026318, 506.596755301209, 485.410196624968, 462.307945665474, 
437.381176452847, 410.728263557213, 382.454393849214, 352.671151375484, 
321.496076987398, 289.052204461029, 255.467574939043, 220.874731610807, 
185.410196624968, 149.213932298913, 112.428788751435, 75.1999401385824, 
37.6743117175879, -1.92977144680695e-13, -37.674311717588, -75.1999401385826, 
-112.428788751435, -149.213932298913, -185.410196624969, -220.874731610807, 
-255.467574939044, -289.052204461029, -321.496076987398, -352.671151375484, 
-382.454393849214, -410.728263557213, -437.381176452847, -462.307945665474, 
-485.410196624968, -506.596755301209, -525.784008026318, -542.896231479612, 
-557.865891532951, -570.633909777092, -581.149896677179, -589.372350437213, 
-595.268820788687, -598.816037056963, -600, -598.816037056963, 
-595.268820788687, -589.372350437213, -581.149896677179, -570.633909777092, 
-557.865891532951, -542.896231479612, -525.784008026318, -506.596755301209, 
-485.410196624968, -462.307945665473, -437.381176452847, -410.728263557213, 
-382.454393849214, -352.671151375484, -321.496076987398, -289.052204461029, 
-255.467574939043, -220.874731610807, -185.410196624968, -149.213932298913, 
-112.428788751435, -75.1999401385823, -37.6743117175880, -1.46952762458685e-13, 
37.674311717588, 598.816037056963, 595.268820788687, 589.372350437213, 
581.149896677179, 570.633909777092, 557.865891532951, 542.896231479612, 
525.784008026318, 506.596755301209, 485.410196624968, 462.307945665473, 
437.381176452847, 410.728263557213, 382.454393849214, 352.671151375484, 
321.496076987398, 289.052204461029, 255.467574939044, 220.874731610807, 
185.410196624968, 149.213932298913, 112.428788751435, 75.1999401385825, 
37.674311717588, -9.64885723403475e-14, -37.6743117175880, -75.1999401385826, 
-112.428788751435, -149.213932298913, -185.410196624969, -220.874731610807, 
-255.467574939044, -289.052204461029, -321.496076987398, -352.671151375484, 
-382.454393849214, -410.728263557213, -437.381176452847, -462.307945665474, 
-485.410196624968, -506.596755301209, -525.784008026318, -542.896231479612, 
-557.865891532951, -570.633909777092, -581.149896677179, -589.372350437213, 
-595.268820788687, -598.816037056963, -600, -598.816037056963, 
-595.268820788687, -589.372350437213, -581.149896677179, -570.633909777092, 
-557.865891532951, -542.896231479612, -525.784008026318, -506.596755301209, 
-485.410196624968, -462.307945665473, -437.381176452847, -410.728263557213, 
-382.454393849214, -352.671151375484, -321.496076987398, -289.052204461029, 
-255.467574939043, -220.874731610807, -185.410196624969, -149.213932298913, 
-112.428788751435, -75.1999401385822, -37.6743117175879, -1.10214571844014e-13, 
37.6743117175882, 75.1999401385825, 112.428788751435, 149.213932298913, 
185.410196624968, 220.874731610807, 255.467574939044, 289.052204461029, 
321.496076987398, 352.671151375484, 382.454393849214, 410.728263557213, 
437.381176452847, 462.307945665474, 485.410196624969, 506.596755301209, 
525.784008026318, 542.896231479612, 557.865891532951, 570.633909777092, 
581.149896677179, 589.372350437213, 595.268820788687, 598.816037056963, 
600, 598.816037056963), .Dim = c(101L, 2L), .Dimnames = list(
    NULL, c("x", "y")))

grd <- GridTopology(cellcentre.offset = c(-600, -600), cellsize = c(1, 1), cells.dim = c(1200, 1200))
obj <- kernel2d(pts = kern.obj, poly = circpol, h0 = 100, nx = 600, ny = 600, kernel='quartic')
plot(kern.obj[, "x"], kern.obj[, "y"], xlim = c(-600, 600), ylim = c(-600, 600))
image(obj, add = TRUE, col = terrain.colors(20))
}

require(plots)
require(ggplot2)
require(RColorBrewer)
source("my.colorFct.R")
MH63 <- read.xls("si.xls", header=T, row.names=1, sheet=1)
SY63 <- read.xls("si.xls", header=T, row.names=1, sheet=2)
ZS97 <- read.xls("si.xls", header=T, row.names=1, sheet=3)
my.heat <- function(input)
{
}
MH63.t <- as.matrix(MH63[,-1])
rownames(MH63.t) <- MH63[,1]
MH63.t <- MH63.t[apply(MH63.t > 100, 1, sum)/length(MH63.t[1,])>0.5 & apply(log2(MH63.t), 1, IQR) > 1.5, ]
MH63.ts <- t(scale(t(MH63.t)))

SY63.t <- as.matrix(SY63[,-1])
rownames(SY63.t) <- SY63[,1]
SY63.t <- SY63.t[apply(SY63.t > 100, 1, sum)/length(SY63.t[1,])>0.5 & apply(log2(SY63.t), 1, IQR) > 1.5, ]
SY63.ts <- t(scale(t(SY63.t)))

ZS97.t <- as.matrix(ZS97[,-1])
rownames(ZS97.t) <- ZS97[,1]
ZS97.t <- ZS97.t[apply(ZS97.t > 100, 1, sum)/length(ZS97.t[1,])>0.5 & apply(log2(ZS97.t), 1, IQR) > 1.5, ]
ZS97.ts <- t(scale(t(ZS97.t)))

pdf("heatmap_all.pdf", height = 10, width = 15)
heatmap.2(t(scale(t(SY63.t))), col=redgreen(75), trace = "none", density.info = "none", keysize = 0.8, key=T, Colv = F,
          lmat = rbind(c(4,0,0), c(2,1,0), c(0,3,0)), lwid = c(3,12,1), lhei = c(2,10,1))
title("SY63 heatmap")

heatmap.2(t(scale(t(MH63.t))), col=redgreen(75), trace = "none", density.info = "none", keysize = 0.8, key=T, Colv = F,
          lmat = rbind(c(4,0,0), c(2,1,0), c(0,3,0)), lwid = c(3,12,1), lhei = c(2,10,1))
title("MH63 heatmap")

heatmap.2(t(scale(t(ZS97.t))), col=redgreen(75), trace = "none", density.info = "none", keysize = 0.8, Colv = F,
          lmat = rbind(c(4,0,0), c(2,1,0), c(0,3,0)), lwid = c(3,12,1), lhei = c(2,10,1))
          ## lmat = rbind(c(0,3,0), c(2,1,0), c(0,4,0)), lwid = c(1,6,2), lhei = c(2,6,2))## RowSideColors=mycolhc)
title("ZS97 heatmap")
dev.off()

## layout and par(cfrow, mar) to set the configuration of graphics
pdf("heatmap_han.pdf", height=15, width = 25)
par(mar=c(20,1,3,3))
layout(matrix(c(4,5,6,1,2,3),2,3,byrow=TRUE), c(6,6,6), c(2,12), TRUE)
## layout.show(nf)
x <- 10*(1:nrow(MH63.ts))
y <- 10*(1:ncol(MH63.ts))
image(x, y, MH63.ts, col=rev(redgreen(75)), axes=FALSE, xlab="", ylab="", main="MH63 heatmap")
axis(1, at=seq(min(x), max(x), length=length(rownames(MH63.ts))), labels = rownames(MH63.ts), las=2, ## side = 2,
     outer=F, tick = F, cex.axis=1.2)
axis(2, at=seq(min(y), max(y), length=length(colnames(MH63.ts))), labels = colnames(MH63.ts), las=1, ## side = 2,
     outer=F, tick = F, cex.axis=1.2)
x <- 10*(1:nrow(SY63.ts))
y <- 10*(1:ncol(SY63.ts))
image(x, y, SY63.ts, col=rev(redgreen(75)), axes=FALSE, xlab="", ylab="", main="SY63 heatmap")
axis(1, at=seq(min(x), max(x), length=length(rownames(SY63.ts))), labels = rownames(SY63.ts), las=2, ## side = 2,
     outer=F, tick = F, cex.axis=1.2)
x <- 10*(1:nrow(ZS97.ts))
y <- 10*(1:ncol(ZS97.ts))
image(x, y, ZS97.ts, col=rev(redgreen(75)), axes=FALSE, xlab="", ylab="", main="ZS97 heatmap")
axis(1, at=seq(min(x), max(x), length=length(rownames(ZS97.ts))), labels = rownames(ZS97.ts), las=2, ## side = 2,
     outer=F, tick = F, cex.axis=1.2)
## image(matrix(1:max(MH63.ts)), col=colorRampPalette(brewer.pal())
## par(mar=c(1,1,9,1))
## image(matrix(1:max(SY63.ts)), col=rev(redgreen(75)), ylim=c(0, 0.1), xlim=c(0, max(SY63.ts)))
par(mar=c(3,5,9,9), mai=c(0.7,1,1,3), lheight = 1) ## mai for inches per side
## image(matrix(1:max(ZS97.ts)), col=rev(brewer.pal(11, "RdBu")))
## image(t(matrix(seq(1,max(ZS97.ts), by=0.1))), col=redgreen(75))
z <- matrix(seq(min(ZS97.ts),max(ZS97.ts), length=10))
x <- seq(round(min(ZS97.ts)-1),(max(ZS97.ts)+1), length=10)
y <- 1:ncol(z)
image(x, y, z, col=rev(redgreen(75)), xlab="color key of Rice Heatmap",ylab="", axes=F, cex=3)
axis(1, at=as.integer(x), las=0, ## side = 2,
     outer=F, tick = F, cex.axis=1.2)
dev.off()

     ## x <- 10*(1:nrow(volcano))
     ## y <- 10*(1:ncol(volcano))
     ## image(x, y, volcano, col = terrain.colors(100), axes = FALSE)
     ## contour(x, y, volcano, levels = seq(90, 200, by = 5),
     ##         add = TRUE, col = "peru")
     ## axis(1, at = seq(100, 800, by = 100))
     ## axis(2, at = seq(100, 600, by = 100))

require(grDevices) # for colours
x <- y <- seq(-4*pi, 4*pi, len=27)
r <- sqrt(outer(x^2, y^2, "+"))
image(z = z <- cos(r^2)*exp(-r/6), col=gray((0:32)/32))
box()
title(main = "Maunga Whau Volcano", font.main = 4)

dat <- as.matrix(dat[,-1])
## data filter
dat <- dat[apply(dat > 100, 1, sum)/length(dat[1,])>0.5 & apply(log2(dat), 1, IQR) > 1.5, ]
mydatascale <- t(scale(t(dat))) # Centers and scales data.
hr <- hclust(as.dist(1-cor(t(mydatascale), method="pearson")), method="complete") # Clusters rows by Pearson correlation as distance method.
hc <- hclust(as.dist(1-cor(mydatascale, method="spearman")), method="complete") # Clusters columns by Spearman correlation as distance method.
par(mfrow = c(1,2), las=0)
image(mydatascale, col=my.colorFct())
image(matrix(1:100),col=rev(colorRampPalette(brewer.pal(11, "RdBu"))(50)))
test <- as.vector(t(mydatascale))
## transform data to 2 matched data.frame
test <- data.frame(expand.grid(y=colnames(mydatascale), x=rownames(mydatascale)),v=test)

# beatiful ggplot heatmap
gplot <- function(input,      ## data input, better in matrix structure
name = "ggplot2_heatmap.pdf", ## the name of the export image
)
{
  cin <- t(scale(t(input)))
  pdf("heatmap_han.pdf", height = 15, width = 25)
  par(mfrow = c(1,2))
  ggplot(test, aes(y, x, fill = v, label = sprintf("%.1f", v)))+
  geom_tile() + geom_text() +
  scale_fill_gradient2(low = "blue", high = "red")
  dev.off()
  heatmap(mydatascale, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=my.colorFct())
  heatmap(mydatascale, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=my.colorFct(), scale = "row")
  mycl <- cutree(hr, h=max(hr$height)/1.5); mycolhc <- sample(rainbow(256));
  mycolhc <- mycolhc[as.vector(mycl)]
  heatmap(mydatascale, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=my.colorFct(), scale="row", RowSideColors=mycolhc,
          cexRow = 1.2, cexCol = 1.2, las=1, symm = T)
  heatmap.2(mydatascale, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col= rev(colorRampPalette(brewer.pal(11, "RdBu"))(100)), scale="row", cellnote=round(mydatascale, 1), notecol = "black", trace = "none", density.info = "none", lmat = rbind(c(0,3,0), c(2,1,0), c(0,4,0)), lwid = c(1,6,2), lhei = c(2,6,2))## RowSideColors=mycolhc
}
