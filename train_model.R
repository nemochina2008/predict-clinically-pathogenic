works_with_R("3.2.3",
             data.table="1.9.7",
             missMDA="1.10",
             glmnet="1.9.5",
             "tdhock/WeightedROC@ef8f35ba7ae85e2995fa66afe13732cebb8b5633",
             xgboost="0.4.4",
             party="1.0.13")

arg.vec <- "data-2016-11-09/train.txt"

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
  "Training using", nrow(train.features),
  "variants and the following", ncol(train.features),
  "numeric features.\n")
print(colnames(train.features))
train.labels <- factor(input.dt$clinvar, c("Benign", "Pathogenic"))

nb <- estim_ncpPCA(train.features)
imputed <- MIPCA(train.features, nb$ncp)
train.imputed <- imputed$res.imputePCA
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
model.list <- list()
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
    train.features, label=train.01, weight=train.weight.vec, missing=NA)
  model.list[[m("xgboost")]] <- xgb.train(
    params=xg.param,
    data=xg.mat,
    nrounds=50
  )
  ## model.list[[m("cforest")]] <- cforest(
  ##   label ~ ., train.df, weights=train.weight.vec)
  model.list[[m("glmnet")]] <- cv.glmnet(
    train.imputed, train.labels, train.weight.vec, family="binomial")
}#weight.name

model.RData <- paste0(train.txt, ".RData")
cat("Saving models to ", model.RData, "\n", sep="")
save(model.list, keep, file=model.RData)
