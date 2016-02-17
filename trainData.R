works_with_R("3.2.3",
             data.table="1.9.7")

BP <- fread("BP_predictors_selected.txt", na.strings="-")

is.input <- sapply(BP, is.numeric)

trainData <- 
  list(feature.mat=as.matrix(BP[, is.input, with=FALSE]),
       label.vec=BP$clinvar)

save(trainData, file="trainData.RData")
