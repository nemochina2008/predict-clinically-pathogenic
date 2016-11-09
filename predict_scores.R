works_with_R("3.2.3",
             data.table="1.9.7",
             missMDA="1.10",
             glmnet="1.9.5",
             "tdhock/WeightedROC@ef8f35ba7ae85e2995fa66afe13732cebb8b5633",
             xgboost="0.4.4",
             party="1.0.13")

arg.vec <- c(
  "data-2016-11-09/train.txt.RData",
  "data-2016-11-09/test1.txt"
  )

arg.vec <- commandArgs(trailingOnly=TRUE)

stopifnot(length(arg.vec)==2)
model.RData <- normalizePath(arg.vec[1], mustWork=TRUE)
test.txt <- normalizePath(arg.vec[2], mustWork=TRUE)
input.dt <- fread(test.txt, na.strings="-")
objs <- load(model.RData)

feature.name.vec <- names(keep)[keep]
feature.in.test <- feature.name.vec %in% names(input.dt)
feature.missing <- feature.name.vec[!feature.in.test]
if(length(feature.missing)){
  stop(model.RData, " needs features listed above to predict, but ",
       test.txt, " does not have ",
       paste(feature.missing, collapse=", "))
}

test.features <- as.matrix(input.dt[, feature.name.vec, with=FALSE])

nb <- estim_ncpPCA(test.features)
test.mipca <- MIPCA(test.features, nb$ncp)
test.imputed <- test.mipca$res.imputePCA
test.df <- data.frame(
  test.features,
  check.names=FALSE)
## prediction function is
## ifelse(is.na(feature, na.guess,
##  ifelse(feature <= thresh, "Benign", "Pathogenic")
is.pathogenic <- 2
pred.fun.list <- list(
  ctree=function(fit){
    pred.list <- treeresponse(fit, test.df)
    sapply(pred.list, "[", is.pathogenic)
  },
  xgboost=function(fit){
    predict(fit, test.features, missing=NA)            
  },
  cforest=function(fit){
    pred.list <- treeresponse(fit, test.df)
    sapply(pred.list, "[", is.pathogenic)
  },
  glmnet=function(fit){
    pred.mat <- predict(fit, test.imputed, type="response")
    as.numeric(pred.mat)
  })
for(model.name in names(model.list)){
  model.prefix <- sub("[.].*", "", model.name)
  pred.fun <- pred.fun.list[[model.prefix]]
  model <- model.list[[model.name]]
  pred.score <- pred.fun(model)
  input.dt[[paste0(model.name, "_score")]] <- pred.score
  input.dt[[paste0(model.name, "_pred")]] <- ifelse(
    pred.score <= 0.5, "Benign", "Pathogenic")
}#model.name

out.txt <- paste0(test.txt, "_predictions.txt")
write.table(
  input.dt,
  out.txt,
  quote=FALSE,
  sep="\t",
  row.names=FALSE,
  col.names=TRUE)
