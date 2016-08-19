load("trainData.RData")
load("testData.RData")
load("folds.RData")
load("testData.folds.RData")

BP.other <- list(
  BP=list(data=trainData, folds=folds),
  other=list(data=testData, folds=testData.folds))

save(BP.other, file="BP.other.RData")
