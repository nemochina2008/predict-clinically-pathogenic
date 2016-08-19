works_with_R("3.2.3",
             "tdhock/WeightedROC@ef8f35ba7ae85e2995fa66afe13732cebb8b5633",
             party="1.0.13")

load("BP.other.RData")

fold.vec <- sort(unique(BP.other[[1]]$folds))
levs <- c("Benign", "Pathogenic")
is.pathogenic <- levs == "Pathogenic"

BP.roc.list <- list()
BP.pred.list <- list()
for(test.fold in fold.vec){
  fit.list <- list()
  for(train.name in names(BP.other)){
    data.folds.list <- BP.other[[train.name]]
    is.train <- data.folds.list$folds != test.fold
    train.df <- with(data.folds.list$data, data.frame(
      label=factor(label.vec[is.train], levs),
      feature.mat[is.train,],
      check.names=FALSE)
      )
    fit.tree <- ctree(label ~ ., train.df)
    response.tab <- table(train.df$label)
    response.dec <- sort(response.tab, decreasing=TRUE)
    major.prob <- ifelse(names(response.dec)[1]=="Pathogenic", 1, 0)
    ##fit.forest <- cforest(label ~ ., train.df)
    single.features.list <- list()
    for(col.name in colnames(data.folds.list$data$feature.mat)){
      one.feature <- train.df[, col.name]
      na.tab <- table(train.df$label[is.na(one.feature)])
      na.guess <- names(na.tab)[which.max(na.tab)]
      not.na <- train.df[!is.na(one.feature),]
      ord <- order(not.na[, col.name])
      train.ord <- data.table(
        label=not.na$label,
        feature=not.na[[col.name]]
      )[ord,]
      train.counts <- train.ord[, list(
        n.Benign=sum(label=="Benign"),
        n.Pathogenic=sum(label=="Pathogenic")
      ), by=feature]
      thresh.dt.list <- list()
      for(feature.sign in c(1, -1)){
        thresh.counts <- train.counts[order(feature.sign*feature),]
        thresh.counts[, cum.TP := cumsum(n.Pathogenic)]
        thresh.counts[, cum.FP := sum(n.Benign)-cumsum(n.Benign)]
        thresh.counts[, n.incorrect := cum.FP + cum.TP]
        thresh.row <- thresh.counts[which.min(n.incorrect),]
        thresh <- feature.sign*thresh.row$feature
        pred.vec <- ifelse(
          feature.sign*one.feature <= thresh,
          "Benign", "Pathogenic")
        thresh.errors <- sum(pred.vec != train.df$label, na.rm=TRUE)
        stopifnot(thresh.errors == thresh.row$n.incorrect)
        thresh.dt.list[[paste(feature.sign)]] <- data.table(
          feature.sign, thresh, col.name,
          n.incorrect=thresh.row$n.incorrect)
      }
      thresh.dt <- do.call(rbind, thresh.dt.list)
      single.features.list[[col.name]] <- data.table(
        thresh.dt[which.min(n.incorrect),],
        na.guess)
    }#col.name (single feature models)
    single.features <- do.call(rbind, single.features.list)
    for(test.name in names(BP.other)){
      cat(sprintf("fold=%d trainSet=%s testSet=%s\n",
                  test.fold, train.name, test.name))
      test.list <- BP.other[[test.name]]
      is.test <- test.list$folds == test.fold
      test.df <- with(test.list$data, data.frame(
        label=factor(label.vec[is.test]),
        feature.mat[is.test,],
        check.names=FALSE)
        )
      ## prediction function is
      ## ifelse(is.na(feature, na.guess,
      ##  ifelse(feature <= thresh, "Benign", "Pathogenic")
      pred.mat.list <- list(
        major.class=rep(major.prob, nrow(test.df)),
        ctree={
          pred.list <- treeresponse(fit.tree, test.df)
          sapply(pred.list, "[", is.pathogenic)
        })
      pred.thresh.list <- list(
        major.class=0.5,
        ctree=0.5)
      tree.pred <- predict(fit.tree, test.df)
      my.pred <- ifelse(
        pred.mat.list$ctree <= pred.thresh.list$ctree,
        "Benign", "Pathogenic")
      stopifnot(tree.pred == my.pred)
      for(feature.i in 1:nrow(single.features)){
        one.feature <- single.features[feature.i,]
        feature.name <- paste(one.feature$col.name)
        test.feature.vec <- test.df[[feature.name]]
        pred.score <- test.feature.vec*one.feature$feature.sign
        na.score <- ifelse(
          one.feature$na.guess=="Pathogenic",
          one.feature$thresh+1,
          one.feature$thresh-1)
        pred.score[is.na(test.feature.vec)] <- na.score
        stopifnot(!is.na(pred.score))
        pred.mat.list[[feature.name]] <- pred.score
        pred.thresh.list[[feature.name]] <- one.feature$thresh
      }
      label.int <- ifelse(test.df$label == "Pathogenic", 1, -1)
      for(model.name in names(pred.mat.list)){
        pred.score <- pred.mat.list[[model.name]]
        roc <- WeightedROC(pred.score, label.int)
        thresh <- pred.thresh.list[[model.name]]
        label.pred <- ifelse(pred.score <= thresh, -1, 1)
        FP <- sum(label.int == -1 & label.pred == 1)
        TP <- sum(label.int == 1 & label.pred == 1)
        n.positive <- sum(label.int == 1)
        n.negative <- sum(label.int == -1)
        FN <- sum(n.positive)-TP
        FPR <- FP/n.negative
        TPR <- TP/n.positive
        errors <- FP + FN
        n.data <- nrow(test.df)
        pred.dot <- data.frame(
          FP, n.negative,
          TP, FN, n.positive, 
          FPR, TPR,
          errors, n.data,
          error.percent=100*errors/n.data)
        ## ggplot()+
        ##   coord_equal()+
        ##   geom_point(aes(FPR, TPR), data=pred.dot)+
        ##   geom_path(aes(FPR, TPR), data=roc)
        auc <- WeightedAUC(roc)
        BP.roc.list[[paste(model.name, test.name, train.name, test.fold)]] <-
          data.table(model.name, test.name, train.name, test.fold,
                     roc)
        BP.pred.list[[paste(model.name, test.name, train.name, test.fold)]] <-
          data.table(model.name, test.name, train.name, test.fold,
                     auc, pred.dot)
      }#model.name
    }#test.name
  }#train.name
}#test.fold
BP.other.models <- list(
  roc=do.call(rbind, BP.roc.list),
  predictions=do.call(rbind, BP.pred.list))
save(BP.other.models, file="BP.other.models.RData")
