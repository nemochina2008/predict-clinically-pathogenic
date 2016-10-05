works_with_R("3.2.3", data.table="1.9.7")
data.txt.vec <- Sys.glob("data/*.txt")
BP.other <- list()
set.seed(1)
n.folds <- 4
for(data.txt in data.txt.vec){
  data.dt <- fread(data.txt, na.strings="-")
  data.name <- sub(".txt$", "", basename(data.txt))
  BP.other[[data.name]] <- list(
    data=list(
      feature.mat=as.matrix(data.dt[, sapply(data.dt, is.numeric), with=FALSE]),
      label.vec=data.dt$clinvar),
    folds=sample(rep(1:n.folds, l=nrow(data.dt))))
}

save(BP.other, file="BP.other.RData")
