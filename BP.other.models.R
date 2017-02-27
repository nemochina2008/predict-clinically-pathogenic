works_with_R("3.3.2",
             data.table="1.9.7",
             missMDA="1.10",
             glmnet="1.9.5",
             "tdhock/WeightedROC@4dc6e54e248dadd2c60f58b9870fd6b1d7e13eef",
             xgboost="0.4.4",
             party="1.0.13")

not.train.features <- c(
  "MetaSVM_score",
  "MetaLR_score",
  "Eigen-raw",
  "Eigen-PC-raw",
  "GenoCanyon_score",
  "MCAP",
  "REVEL",
  "esp6500siv2_all", "1000g2015aug_all")

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
    all.train.features <- data.folds.list$data$feature.mat[is.train,]
    keep <- ! colnames(all.train.features) %in% not.train.features
    train.features <- all.train.features[, keep]
    n.train.features <- ncol(train.features)
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
      xg.mat <- xgb.DMatrix(
        train.features, label=train.01, weight=train.weight.vec, missing=NA)
      fit.xgboost <- xgb.train(
        params=xg.param,
        data=xg.mat,
        nrounds=50
        )
      fit.forest <- cforest(label ~ ., train.df, weights=train.weight.vec)
      fit.glmnet <- cv.glmnet(
        train.imputed, train.labels, train.weight.vec, family="binomial")
      single.features.list <- list()
      for(col.name in colnames(data.folds.list$data$feature.mat)){
        train.and.weights <- data.table(
          feature=all.train.features[, col.name],
          label=train.df$label,
          weight=train.weight.vec)
        na.dt <- train.and.weights[is.na(feature),]
        na.guess <- if(nrow(na.dt)){
          na.class.weights <- na.dt[, list(weight=sum(weight)), by=label]
          paste(na.class.weights[which.max(weight), label])
        }else{
          "Benign"
        }
        not.na <- train.and.weights[!is.na(feature),]
        ord <- order(not.na$feature)
        train.ord <- data.table(
          label=not.na$label,
          feature=not.na$feature,
          weight=not.na$weight
          )[ord,]
        train.counts <- dcast(
          train.ord,
          feature ~ label,
          fun.aggregate=sum,
          value.var="weight")
        thresh.dt.list <- list()
        for(feature.sign in c(1, -1)){
          thresh.counts <- train.counts[order(feature.sign*feature),]
          ## Low scores are Benign, high scores are Pathogenic.
          thresh.counts[, cum.FN := cumsum(Pathogenic)]
          thresh.counts[, cum.FP := sum(Benign)-cumsum(Benign)]
          thresh.counts[, n.incorrect := cum.FP + cum.FN]
          thresh.row <- thresh.counts[which.min(n.incorrect),]
          thresh <- feature.sign*thresh.row$feature
          score.vec <- feature.sign*not.na$feature
          pred.vec <- ifelse(
            score.vec <= thresh,
            "Benign", "Pathogenic")
          range.vec <- range(score.vec)
          names(range.vec) <- c("Benign", "Pathogenic")
          is.incorrect <- pred.vec != not.na$label
          thresh.errors <- sum(is.incorrect * not.na$weight, na.rm=TRUE)
          stopifnot(thresh.errors == thresh.row$n.incorrect)
          thresh.dt.list[[paste(feature.sign)]] <- data.table(
            feature.sign, thresh, col.name,
            na.score=range.vec[[na.guess]],
            n.incorrect=thresh.row$n.incorrect)
        }
        thresh.dt <- do.call(rbind, thresh.dt.list)
        single.features.list[[col.name]] <- data.table(
          thresh.dt[which.min(n.incorrect),],
          na.count=nrow(na.dt),
          na.guess)
      }#col.name (single feature models)
      single.features <- do.call(rbind, single.features.list)
      for(test.name in names(BP.other)){
        cat(sprintf("fold=%d trainSet=%s weights=%s testSet=%s\n",
                    test.fold, train.name, weight.name, test.name))
        test.list <- BP.other[[test.name]]
        is.test <- test.list$folds == test.fold
        all.test.features <- test.list$data$feature.mat[is.test, ]
        test.features <- all.test.features[, colnames(train.features)]
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
          test.feature.vec <- all.test.features[, feature.name]
          pred.score <- ifelse(
            is.na(test.feature.vec),
            one.feature$na.score,
            test.feature.vec*one.feature$feature.sign)
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
          key <- paste(
            model.name, test.name, train.name, test.fold, weight.name)
          BP.roc.list[[key]] <- data.table(
            model.name, test.name, train.name, test.fold, weight.name,
            roc)
          BP.pred.list[[key]] <- data.table(
            model.name, test.name, train.name, test.fold, weight.name,
            auc, pred.dot)
        }#model.name
      }#test.name
    }#weight.name
  }#train.name
}#test.fold

BP.other.models <- list(
  roc=do.call(rbind, BP.roc.list),
  predictions=do.call(rbind, BP.pred.list))
save(BP.other.models, file="BP.other.models.RData")
