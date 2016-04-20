works_with_R("3.2.2", 
             data.table="1.9.6",
             party="1.0.25")


load("trainData.RData")
train.features <- trainData$feature.mat
train.labels <- trainData$label.vec
train.df <- data.frame(
  label=factor(train.labels), train.features, check.names=FALSE)

fit.forest <- cforest(label ~ ., train.df)

## load your actual test data
test.dt <- fread("test3.txt")
score.name.vec <- grep("[sS]core", names(test.dt), value=TRUE)
dash.name.vec <- gsub(" ", "-", score.name.vec)
rbind(score.name.vec, dash.name.vec)
for(col.i in seq_along(score.name.vec)){
  old.name <- score.name.vec[[col.i]]
  new.name <- dash.name.vec[[col.i]]
  test.dt[[new.name]] <- as.numeric(test.dt[[old.name]])
}
test.dt$predicted.prob.pathogenic <- predict(fit.forest, test.dt, OOB=TRUE) ##PRedict
