works_with_R("3.2.3", data.table="1.9.7")
data.txt.vec <- Sys.glob("data/*.txt")
BP.other <- list()
set.seed(1)
n.folds <- 4
ignore.features <- c(
  "integrated_confidence_value",
  "esp6500siv2_all",
  "1000g2015aug_all",
  "SubmitterCategories")
for(data.txt in data.txt.vec){
  data.dt <- fread(data.txt, na.strings="-")
  data.name <- sub(".txt$", "", basename(data.txt))
  keep <- sapply(data.dt, is.numeric)
  keep[names(data.dt) %in% ignore.features] <- FALSE
  BP.other[[data.name]] <- list(
    data=list(
      feature.mat=as.matrix(data.dt[, keep, with=FALSE]),
      label.vec=data.dt$clinvar),
    folds=sample(rep(1:n.folds, l=nrow(data.dt))))
}
save(BP.other, file="BP.other.RData")
