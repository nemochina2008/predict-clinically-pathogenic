load("testData.RData")

N <- nrow(trainData$feature.mat)
n.folds <- 5

set.seed(1)
testData.folds <- sample(rep(1:n.folds, l=N))

save(testData.folds, file="testData.folds.RData")
