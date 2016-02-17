load("trainData.RData")


N <- nrow(trainData$feature.mat)
n.folds <- 5

set.seed(1)
folds <- sample(rep(1:n.folds, l=N))

save(folds, file="folds.RData")
