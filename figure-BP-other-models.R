works_with_R("3.2.3", scales="0.4.0", ggplot2="2.1.0", data.table="1.9.7")

load("BP.other.models.RData")
load("BP.other.RData")

trainData <- BP.other[[1]]$data
getFeatures <- function(model.name){
  factor(ifelse(
    model.name=="major.class", 0,
    ifelse(
      model.name %in% c("ctree", "cforest", "glmnet", "xgboost"), ncol(trainData$feature.mat),
      ifelse(grepl("[(]", model.name), 2, 1))))
}

pred.dt <- BP.other.models$predictions
pred.dt[, n.features := getFeatures(model.name)]

roc.dt <- BP.other.models$roc
roc.dt[, n.features := getFeatures(model.name)]

feature.colors <- c(
  "0"="#E41A1C",
  "1"="#377EB8",
  "2"="black",
  "#4DAF4A",
  "#984EA3", "#FF7F00", "#FFFF33", 
  "#A65628", "#F781BF", "#999999")
names(feature.colors)[8] <- ncol(trainData$feature.mat)

model.means <- pred.dt[, list(
  mean.auc=mean(auc),
  mean.error=mean(error.percent)
), by=.(model.name, n.features)]
ggplot()+
  scale_color_manual(values=feature.colors)+
  geom_point(aes(mean.error, mean.auc, color=n.features),
             data=model.means)

set.means <- pred.dt[, list(
  mean.auc=mean(auc),
  mean.error=mean(error.percent)
), by=.(model.name, n.features, train.name, test.name)]
set.means[, rank := rank(mean.auc), by=.(train.name, test.name)]
model.ranks <- set.means[, list(
  mean.rank=mean(rank)
  ), by=.(model.name, n.features)]
best.auc <- set.means[,.SD[which.max(mean.auc), .(test.name, mean.auc)] ,by=test.name]
best.error <- set.means[,.SD[which.min(mean.error), .(test.name, mean.error)] ,by=test.name]
levs <- model.means[order(mean.error), model.name]
levs <- model.ranks[order(mean.rank), model.name]
pred.dt[, model.fac := factor(model.name, levs)]
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(train.name ~ test.name, labeller=label_both)+
  scale_color_manual(values=feature.colors)+
  geom_vline(aes(xintercept=mean.error),
             color="grey",
             data=best.error)+
  geom_point(aes(error.percent, model.fac, color=n.features, fill=weight.name),
             shape=1,
             data=pred.dt)

gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(train.name ~ test.name, labeller=label_both)+
  scale_color_manual(values=feature.colors)+
  geom_vline(aes(xintercept=mean.auc),
             color="grey",
             data=best.auc)+
  scale_fill_manual(values=c(balanced=NA, one="black"))+
  geom_point(aes(auc, model.fac, color=n.features, fill=weight.name),
             shape=21,
             data=pred.dt)
png("figure-BP-other-models-four.png", width=900, height=750)
print(gg)
dev.off()

## pred.dt[, accuracy.percent := 100-error.percent]
## wide.accuracy <- dcast(
##   pred.dt,
##   test.name + test.fold + model.name + n.features ~ train.name,
##   value.var="accuracy.percent")
## same.other <- data.table(
##   side=c("higher", "lower"),
##   same.auc=c(0.65, 0.85),
##   other.auc=c(0.9, 0.6),
##   same.accuracy=c(65, 85),
##   other.accuracy=c(90, 60))
## gg.acc <- ggplot()+
##   geom_text(aes(
##     same.accuracy, other.accuracy,
##     label=paste(side, "accuracy\nfor model trained\non other data set")),
##     color="grey",
##     data=same.other)+
##   geom_abline(slope=1, intercept=0, color="grey")+
##   theme_bw()+
##   theme(panel.margin=grid::unit(0, "lines"))+
##   facet_grid(. ~ test.name, labeller=label_both)+
##   scale_color_manual(values=feature.colors)+
##   geom_point(aes(
##     ifelse(test.name=="BP", BP, other),
##     ifelse(test.name=="BP", other, BP),
##     color=n.features),
##     shape=1,
##     data=wide.accuracy)+
##   coord_equal()+
##   xlab("test accuracy when training on the same data set")+
##   ylab("test accuracy when training on the other data set")
## wide.auc <- dcast(
##   pred.dt,
##   test.name + test.fold + model.name + n.features ~ train.name,
##   value.var="auc")
## png("figure-BP-other-models-accuracy.png", h=300)
## print(gg.acc)
## dev.off()

## gg.auc <- ggplot()+
##   geom_text(aes(
##     same.auc, other.auc,
##     label=paste(side, "AUC\nfor model trained\non other data set")),
##     color="grey",
##     data=same.other)+
##   geom_abline(slope=1, intercept=0, color="grey")+
##   theme_bw()+
##   theme(panel.margin=grid::unit(0, "lines"))+
##   facet_grid(. ~ test.name, labeller=label_both)+
##   scale_color_manual(values=feature.colors)+
##   geom_point(aes(
##     ifelse(test.name=="BP", BP, other),
##     ifelse(test.name=="BP", other, BP),
##     color=n.features),
##     shape=1,
##     data=wide.auc)+
##   coord_equal()+
##   xlab("test AUC when training on the same data set")+
##   ylab("test AUC when training on the other data set")
## png("figure-BP-other-models-auc.png", h=300)
## print(gg.auc)
## dev.off()

## library(animint)
## viz <- list(
##   accuracy=ggplot()+
##     geom_text(aes(
##       same.accuracy, other.accuracy,
##       label=paste(side, "accuracy\nfor model trained\non other data set")),
##       color="grey",
##       data=same.other)+
##     geom_abline(slope=1, intercept=0, color="grey")+
##     theme_bw()+
##     theme(panel.margin=grid::unit(0, "lines"))+
##     facet_grid(. ~ test.name, labeller=label_both)+
##     scale_color_manual(values=model.colors)+
##     guides(color="none")+
##     geom_point(aes(
##       ifelse(test.name=="BP", BP, other),
##       ifelse(test.name=="BP", other, BP),
##       clickSelects=test.fold,
##       tooltip=paste(model.name, "test fold", test.fold),
##       showSelected=model.name,
##       color=model.name),
##       shape=1,
##       size=3,
##       data=wide.accuracy)+
##     theme_animint(height=250)+
##     xlab("test accuracy, trained on same data")+
##     ylab("test accuracy, trained on other data"),
##   auc=ggplot()+
##     geom_text(aes(
##       same.auc, other.auc,
##       label=paste(side, "AUC\nfor model trained\non other data set")),
##       color="grey",
##       data=same.other)+
##     geom_abline(slope=1, intercept=0, color="grey")+
##     theme_bw()+
##     theme(panel.margin=grid::unit(0, "lines"))+
##     facet_grid(. ~ test.name, labeller=label_both)+
##     scale_color_manual(values=model.colors)+
##     guides(color="none")+
##     geom_point(aes(
##       ifelse(test.name=="BP", BP, other),
##       ifelse(test.name=="BP", other, BP),
##       tooltip=paste(model.name, "test fold", test.fold),
##       clickSelects=test.fold,
##       showSelected=model.name,
##       color=model.name),
##       shape=1,
##       size=3,
##       data=wide.auc)+
##     scale_x_continuous(
##       "test AUC, trained on same data",
##       breaks=c(0.5, 0.75, 1),
##       labels=c("0.5", "0.75", "1"))+
##     theme_animint(height=250)+
##     scale_y_continuous(
##       "test AUC, trained on other data",
##       breaks=c(0.5, 0.75, 1),
##       labels=c("0.5", "0.75", "1")),
##   select=ggplot()+
##     theme_bw()+
##     geom_segment(aes(x, model.fac, xend=xend, yend=model.fac,
##                      clickSelects=model.name),
##                  data=seg.dt,
##                  size=12,
##                  alpha=0.1)+
##     theme(panel.margin=grid::unit(0, "lines"))+
##     facet_grid(. ~ train.name + test.name, labeller=function(df){
##       df$train.name <- paste0("train=", df$train.name)
##       df$test.name <- paste0("test=", df$test.name)
##       df
##     })+
##     guides(color="none")+
##     scale_color_manual(values=model.colors)+
##     xlab("Test accuracy (percent correct labels)")+
##     ylab("model")+
##     geom_point(aes(100-error.percent, model.fac, color=model.fac,
##                    clickSelects=test.fold),
##                shape=1,
##                alpha=0.7,
##                size=4,
##                data=pred.dt)+
##     theme_animint(width=700, height=300),
##   ## roc=ggplot()+
##   ##   scale_color_manual(values=model.colors)+
##   ##   geom_path(aes(FPR, TPR, color=model.name,
##   ##                 showSelected=test.fold,
##   ##                 group=paste(model.name,
##   ##                 key=paste(model.name),
##   ##             data=roc.dt),
##   selector.types=list(model.name="multiple"),
##   first=list(model.name=c("VEST3_score", "major.class", "ctree")))
## animint2dir(viz, "figure-BP-other-models")

