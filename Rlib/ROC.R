require(ROCR)
require(pROC)
require(Daim)

data.frame(ROCR.simple)

rocr <- function(predicts, labels){
  pred <- prediction(predicts, labels)
  perf <- performance(pred, "tpr", "fpr")
  plot(perf)
}