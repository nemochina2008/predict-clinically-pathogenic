works_with_R("3.2.3",
             data.table="1.9.7",
             missMDA="1.10",
             glmnet="1.9.5",
             "tdhock/WeightedROC@ef8f35ba7ae85e2995fa66afe13732cebb8b5633",
             xgboost="0.4.4",
             party="1.0.13")

load("BP.other.RData")

levs <- c("Benign", "Pathogenic")
is.pathogenic <- levs == "Pathogenic"

BP.roc.list <- list()
BP.pred.list <- list()
learned.score.list <- list()
learned.thresh.list <- list()
for(train.name in names(BP.other)){
  data.folds.list <- BP.other[[train.name]]
  is.train <- TRUE
  train.features <- data.folds.list$data$feature.mat[is.train,]
  nb <- estim_ncpPCA(train.features)
  imputed <- MIPCA(train.features, nb$ncp)
  train.imputed <- imputed$res.imputePCA
  train.labels <- factor(data.folds.list$data$label.vec[is.train], levs)
  train.01 <- ifelse(train.labels=="Pathogenic", 1, 0)
  train.df <- data.frame(
    label=train.labels,
    train.features,
    check.names=FALSE)
  response.tab <- table(train.labels)
  other.tab <- response.tab
  names(other.tab) <- rev(names(response.tab))
  weight.list <- list(
    balanced=as.numeric(other.tab[paste(train.df$label)]),
    one=rep(1, nrow(train.df)))
  for(weight.name in names(weight.list)){
    train.weight.vec <- weight.list[[weight.name]]
    fit.tree <- ctree(label ~ ., train.df, weights=train.weight.vec)
    response.dec <- sort(response.tab, decreasing=TRUE)
    major.prob <- ifelse(names(response.dec)[1]=="Pathogenic", 1, 0)
    xg.param <- list(
      objective="binary:logistic",
      eval_metric="auc"
      )
    xg.mat <- xgb.DMatrix(train.features, label=train.01, weight=train.weight.vec, missing=NA)
    fit.xgboost <- xgb.train(
      params=xg.param,
      data=xg.mat,
      nrounds=50
      )
    fit.forest <- cforest(label ~ ., train.df, weights=train.weight.vec)
    fit.glmnet <- cv.glmnet(train.imputed, train.labels, train.weight.vec, family="binomial")
    single.features.list <- list()
    for(col.name in colnames(data.folds.list$data$feature.mat)){
      one.feature <- train.df[, col.name]
      train.and.weights <- data.table(train.df, weight=train.weight.vec)
      na.dt <- train.and.weights[is.na(one.feature),]
      na.guess <- if(nrow(na.dt)){
        na.class.weights <- na.dt[, list(weight=sum(weight)), by=label]
        paste(na.class.weights[which.max(weight), label])
      }else{
        "Benign"
      }
      not.na <- train.and.weights[!is.na(one.feature),]
      ord <- order(not.na[, col.name, with=FALSE])
      train.ord <- data.table(
        label=not.na$label,
        feature=not.na[[col.name]],
        weight=not.na$weight
        )[ord,]
      train.counts <- train.ord[, list(
        n.Benign=.SD[label=="Benign", sum(weight)],
        n.Pathogenic=.SD[label=="Pathogenic", sum(weight)]
        ), by=feature]
      thresh.dt.list <- list()
      for(feature.sign in c(1, -1)){
        thresh.counts <- train.counts[order(feature.sign*feature),]
        ## Low scores are Benign, high scores are Pathogenic.
        thresh.counts[, cum.FN := cumsum(n.Pathogenic)]
        thresh.counts[, cum.FP := sum(n.Benign)-cumsum(n.Benign)]
        thresh.counts[, n.incorrect := cum.FP + cum.FN]
        thresh.row <- thresh.counts[which.min(n.incorrect),]
        thresh <- feature.sign*thresh.row$feature
        pred.vec <- ifelse(
          feature.sign*one.feature <= thresh,
          "Benign", "Pathogenic")
        is.incorrect <- pred.vec != train.df$label
        thresh.errors <- sum(is.incorrect * train.weight.vec, na.rm=TRUE)
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
      cat(sprintf("trainSet=%s weights=%s testSet=%s\n",
                  train.name, weight.name, test.name))
      test.list <- BP.other[[test.name]]
      is.test <- TRUE
      test.features <- test.list$data$feature.mat[is.test,]
      nb <- estim_ncpPCA(test.features)
      test.mipca <- MIPCA(test.features, nb$ncp)
      test.imputed <- test.mipca$res.imputePCA
      test.df <- data.frame(
        label=factor(test.list$data$label.vec[is.test]),
        test.features,
        check.names=FALSE)
      ## prediction function is
      ## ifelse(is.na(feature, na.guess,
      ##  ifelse(feature <= thresh, "Benign", "Pathogenic")
      pred.mat.list <- list(
        major.class=rep(major.prob, nrow(test.df)),
        ctree={
          pred.list <- treeresponse(fit.tree, test.df)
          sapply(pred.list, "[", is.pathogenic)
        },
        xgboost={
          predict(fit.xgboost, test.features, missing=NA)            
        },
        cforest={
          pred.list <- treeresponse(fit.forest, test.df)
          sapply(pred.list, "[", is.pathogenic)
        },
        glmnet={
          pred.mat <- predict(fit.glmnet, test.imputed, type="response")
          as.numeric(pred.mat)
        })
      pred.thresh.list <- list(
        xgboost=0.5,
        major.class=0.5,
        ctree=0.5,
        cforest=0.5,
        glmnet=0.5)
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
        if(is.null(thresh)){
          stop("no threshold for model ", model.name)
        }
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
        BP.roc.list[[paste(model.name, test.name, train.name, weight.name)]] <-
          data.table(model.name, test.name, train.name, weight.name,
                     roc)
        BP.pred.list[[paste(model.name, test.name, train.name, weight.name)]] <-
          data.table(model.name, test.name, train.name, weight.name,
                     auc, pred.dot)
        learned.score.list[[paste(model.name, test.name, train.name, weight.name)]] <-
          data.table(model.name, test.name, train.name, weight.name,
                     row.i=seq_along(pred.score),
                     score=pred.score)
      }#model.name
      learned.thresh.list[[paste(test.name, train.name, weight.name)]] <-
        data.table(test.name, train.name, weight.name,
                   model.name=names(pred.thresh.list),
                   largest.Benign=unlist(pred.thresh.list))
    }#test.name
  }#weight.name
}#train.name

predictions.no.folds <- list(
  roc=do.call(rbind, BP.roc.list),
  error=do.call(rbind, BP.pred.list),
  score=do.call(rbind, learned.score.list),
  threshold=do.call(rbind, learned.thresh.list))

## 2013-2016 is supposed to be the training dataset and 2016 the test one.
combination.dt <- unique(predictions.no.folds$score[, .(train.name, test.name, weight.name)])
setkey(predictions.no.folds$score, train.name, test.name, weight.name)
dir.create("predicted-scores")
for(row.i in 1:nrow(combination.dt)){
  one.combo <- combination.dt[row.i,]
  score.dt <- dcast(predictions.no.folds$score[one.combo], row.i ~ model.name, value.var="score")
  name.vec <- sub(".name", "", names(one.combo))
  value.vec <- unlist(one.combo)
  name.value.vec <- paste0(name.vec, "=", value.vec)
  base <- paste(name.value.vec, collapse="_")
  base.csv <- paste0(base, ".csv")
  out.path <- file.path("predicted-scores", base.csv)
  write.table(
    score.dt,
    out.path,
    sep=",",
    row.names=FALSE,
    col.names=TRUE,
    quote=FALSE)
  fread(out.path)
}

save(predictions.no.folds, file="predictions.no.folds.RData")
