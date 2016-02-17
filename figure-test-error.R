works_with_R("3.2.3", ggplot2="1.0.1")

load("models.RData")

gg <- ggplot()+
  xlab("percent test error (5-fold cross-validation)")+
  geom_point(aes(100*test.errors/test.labels, model.name),
             shape=1,
             data=models)

png("figure-test-error.png")
print(gg)
dev.off()
