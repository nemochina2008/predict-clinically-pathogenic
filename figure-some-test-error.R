works_with_R("3.2.3", ggplot2="1.0.1", data.table="1.9.7")

load("some.models.RData")

models.dt <- data.table(some.models)
models.dt[, percent.test.error := 100*test.errors/test.labels]
stats.dt <- models.dt[, list(mean=mean(percent.test.error)), by=model.name]
order.dt <- stats.dt[order(mean, decreasing=TRUE),]
models.dt[, model.fac := factor(model.name, order.dt$model.name)]
models.dt[, n.features := factor(ifelse(
  model.name=="major.class", 0,
  ifelse(
    model.name %in% c("ctree", "cforest"), ncol(trainData$feature.mat),
    ifelse(grepl("[(]", model.name), 2, 1))))]

feature.colors <- c(
  "0"="#E41A1C",
  "1"="#377EB8",
  "2"="black",
  "#4DAF4A",
  "#984EA3", "#FF7F00", "#FFFF33", 
  "#A65628", "#F781BF", "#999999")
names(feature.colors)[8] <- ncol(trainData$feature.mat)

gg <- ggplot()+
  theme(legend.position="top")+
  scale_color_manual(values=feature.colors)+
  xlab("percent test error (5-fold cross-validation)")+
  ylab("model")+
  geom_point(aes(percent.test.error, model.fac, color=n.features),
             shape=1,
             data=models.dt)

png("figure-some-test-error.png", 400, 300)
print(gg)
dev.off()
