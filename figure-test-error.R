works_with_R("3.2.3", ggplot2="1.0.1", data.table="1.9.7")

load("models.RData")

models.dt <- data.table(models)
models.dt[, percent.test.error := 100*test.errors/test.labels]
stats.dt <- models.dt[, list(mean=mean(percent.test.error)), by=model.name]
order.dt <- stats.dt[order(mean, decreasing=TRUE),]
models.dt[, model.fac := factor(model.name, order.dt$model.name)]

gg <- ggplot()+
  xlab("percent test error (5-fold cross-validation)")+
  geom_point(aes(percent.test.error, model.fac),
             shape=1,
             data=models.dt)

png("figure-test-error.png")
print(gg)
dev.off()
