works_with_R("3.3.3",
             data.table="1.10.4",
             missMDA="1.10",
             glmnet="1.9.5",
             "tdhock/WeightedROC@4dc6e54e248dadd2c60f58b9870fd6b1d7e13eef",
             xgboost="0.4.4",
             party="1.2.2")

arg.vec <- "data-2016-11-09/train.txt"
arg.vec <- "data-minusHGMD/2013to2016minus-HGMD.txt"

not.train.features <- c(
  "MetaSVM_score",
  "MetaLR_score",
  "Eigen-raw",
  "Eigen-PC-raw",
  "GenoCanyon_score",
  "VEST3_score",
  "FATHMM_score",
  "M-CAP_score",
  "MutationTaster_score",
  "fathmm-MKL_coding_score",
  "MCAP",
  "REVEL",
  "esp6500siv2_all", "1000g2015aug_all")

arg.vec <- commandArgs(trailingOnly=TRUE)

stopifnot(length(arg.vec)==1)
train.txt <- normalizePath(arg.vec, mustWork=TRUE)
input.dt <- fread(train.txt, na.strings="-")
if(! "clinvar" %in% names(input.dt)){
  stop("no clinvar column found in ", train.txt)
}
not.recognized <- input.dt[! clinvar %in% c("Pathogenic", "Benign"),]
if(nrow(not.recognized)){
  print(not.recognized)
  stop("clinvar should be either Pathogenic or Benign")
}
keep <- sapply(input.dt, is.numeric)
train.features <- as.matrix(
  input.dt[, keep, with=FALSE])
cat(
  "The following", sum(!keep),
  "non-numeric columns will not be used as input features.\n")
print(names(input.dt)[!keep])
cat(
  "Found", nrow(train.features),
  "variants and the following", ncol(train.features),
  "numeric features.\n")
print(colnames(train.features))
train.labels <- factor(input.dt$clinvar, c("Benign", "Pathogenic"))

## This logical vector determines which training features we will use
## for the multivariate models.
not.multivar <- colnames(train.features) %in% not.train.features
is.multivar <- !not.multivar
cat(
  "Not using the following", sum(not.multivar),
  "features for multivariate model training.\n")
print(colnames(train.features)[not.multivar])
feature.name.vec <- colnames(train.features)[is.multivar]

nb <- estim_ncpPCA(train.features)
imputed <- MIPCA(train.features[, is.multivar], nb$ncp)
train.imputed <- imputed$res.imputePCA
train.01 <- ifelse(train.labels=="Pathogenic", 1, 0)
train.df <- data.frame(
  label=train.labels,
  train.features[, is.multivar],
  check.names=FALSE)
response.tab <- table(train.labels)
other.tab <- response.tab
names(other.tab) <- rev(names(response.tab))
weight.list <- list(
  balanced=as.numeric(other.tab[paste(train.df$label)]),
  one=rep(1, nrow(train.df)))
model.list <- list()
univariate.model.list <- list()
for(weight.name in names(weight.list)){
  train.weight.vec <- weight.list[[weight.name]]
  m <- function(model.name){
    paste0(model.name, ".weights=", weight.name)
  }
  model.list[[m("ctree")]] <- ctree(
    label ~ ., train.df, weights=train.weight.vec)
  xg.param <- list(
    objective="binary:logistic",
    eval_metric="auc"
  )
  xg.mat <- xgb.DMatrix(
    train.features[, is.multivar],
    label=train.01, weight=train.weight.vec, missing=NA)
  model.list[[m("xgboost")]] <- xgb.train(
    params=xg.param,
    data=xg.mat,
    nrounds=50
  )
  model.list[[m("cforest")]] <- cforest(
    label ~ ., train.df, weights=train.weight.vec)
  model.list[[m("glmnet")]] <- cv.glmnet(
    train.imputed, train.labels, train.weight.vec, family="binomial")
  single.features.list <- list()
  for(col.name in colnames(train.features)){
    one.feature <- train.features[, col.name]
    train.and.weights <- data.table(
      label=train.labels,
      train.features[, col.name, drop=FALSE], weight=train.weight.vec)
    na.dt <- train.and.weights[is.na(one.feature),]
    na.guess <- if(nrow(na.dt)){
      na.class.weights <- na.dt[, list(weight=sum(weight)), by=label]
      paste(na.class.weights[which.max(weight), label])
    }else{
      "Benign"
    }
    not.na <- train.and.weights[!is.na(one.feature),]
    not.na.feature <- not.na[[col.name]]
    ord <- order(not.na.feature)
    train.ord <- data.table(
      label=not.na$label,
      feature=not.na.feature,
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
      score.vec <- feature.sign*not.na.feature
      pred.vec <- ifelse(
        score.vec <= thresh,
        "Benign", "Pathogenic")
      range.vec <- range(score.vec)
      names(range.vec) <- c("Benign", "Pathogenic")
      is.incorrect <- pred.vec != not.na$label
      thresh.errors <- sum(is.incorrect * not.na$weight)
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
  univariate.model.list[[weight.name]] <- single.features
}#weight.name

model.RData <- paste0(train.txt, ".RData")
cat("Saving models to ", model.RData, "\n", sep="")
save(model.list, univariate.model.list, feature.name.vec, file=model.RData)
