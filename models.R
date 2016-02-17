works_with_R("3.2.3", glmnet="1.9.5", caret="6.0.41",
             party="1.0.13")

load("trainData.RData")
load("folds.RData")

caret.methods <- c("knn", "rf", "svmRadial", "gbm", "glmnet")
caret.methods <- "rf"
all.methods <- c(caret.methods, "glmnetBinDev", "glmnetAcc")

test.error.list <- list()
for(test.fold in unique(folds)){
  cat(sprintf("%4d / %4d folds\n", test.fold, length(unique(folds))))
  is.test <- folds==test.fold
  is.train <- !is.test
  train.features <- trainData$feature.mat[is.train,]
  train.labels <- trainData$label.vec[is.train]
  train.df <- data.frame(
    label=factor(train.labels), train.features, check.names=FALSE)
  response.tab <- table(train.labels)
  response.dec <- sort(response.tab, decreasing=TRUE)
  one.tab <- response.tab
  one.tab[] <- 1
  weight.list <-
    list(one=one.tab,
         balanced=response.tab)
  sapply(weight.list, sum)
  test.features <- trainData$feature.mat[is.test,]
  test.labels <- trainData$label.vec[is.test]
  test.df <- data.frame(test.features, check.names=FALSE)
  fit.tree <- ctree(label ~ ., train.df)
  fit.forest <- cforest(label ~ ., train.df)
  test.prediction.list <-
    list(major.class=rep(names(response.dec)[1], nrow(test.features)),
         ctree=predict(fit.tree, test.df),
         cforest=predict(fit.forest, test.df, OOB=TRUE))
  if(FALSE){ #how to deal with missing values in glmnet?
    for(weight.name in names(weight.list)){
      weight.tab <- weight.list[[weight.name]]
      ignore <- rowSums(is.na(train.features)) != 0
      not.na.features <- train.features[!ignore,]
      not.na.labels <- train.labels[!ignore]
      weight.vec <- 1/weight.tab[paste(not.na.labels)]
      glmnet.accuracy <-
        cv.glmnet(not.na.features, not.na.labels, weight.vec,
                  family="binomial", nfolds=5, type.measure="class")
      test.prediction.list[[paste0("glmnetAcc.", weight.name)]] <-
        predict(glmnet.accuracy, test.features, type="response", s="lambda.min")
      glmnet.bindev <-
        cv.glmnet(not.na.features, not.na.labels, weight.vec,
                  family="binomial", nfolds=5)
      test.prediction.list[[paste0("glmnetBinDev.", weight.name)]] <-
        predict(glmnet.bindev, test.features, type="response")
    }
  }
  for(model.name in names(test.prediction.list)){
    pred.labels <- test.prediction.list[[model.name]]
    is.negative <- test.labels=="Benign"
    is.positive <- test.labels=="Pathogenic"
    is.fp <- pred.labels=="Pathogenic" & is.negative
    is.fn <- pred.labels=="Benign" & is.positive
    is.incorrect <- is.fp | is.fn
    test.errors <- sum(is.incorrect)
    fn <- sum(is.fn)
    possible.fn <- sum(is.positive)
    fp <- sum(is.fp)
    possible.fp <- sum(is.negative)
    test.error.list[[paste(test.fold, model.name)]] <- data.frame(
      test.fold, model.name,
      test.errors, test.labels=length(test.labels),
      FPR=fp/possible.fp,
      FNR=fn/possible.fn)
  }
}

models <- do.call(rbind, test.error.list)

save(models, file="models.RData")

