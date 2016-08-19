works_with_R("3.2.3", data.table="1.9.7")

test.dt <- fread("testdata.txt", na.strings="-")

is.input <- sapply(test.dt, is.numeric)

testData <- 
  list(feature.mat=as.matrix(test.dt[, is.input, with=FALSE]),
       label.vec=test.dt$clinvar)

save(testData, file="testData.RData")
