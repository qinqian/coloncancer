require(ROCR)
require(pROC)
require(Daim)

perf <- performance( pred, "tpr", "fpr" )

plot(perf, main="ROCR fingerpainting toolkit", colorize=TRUE,
     xlab="Marys axis", ylab="", box.lty=7, box.lwd=5,
     box.col="gold", lwd=17, colorkey.relwidth=0.5, xaxis.cex.axis=2,
     xaxis.col='blue', xaxis.col.axis="blue", yaxis.col='green', yaxis.cex.axis=2,
     yaxis.at=c(0,0.5,0.8,0.85,0.9,1), yaxis.las=1, xaxis.lwd=2, yaxis.lwd=3,
     yaxis.col.axis="orange", cex.lab=2, cex.main=2)