works_with_R("3.2.3", scales="0.4.0", ggplot2="2.1.0", data.table="1.9.7")

load("BP.other.models.RData")
load("BP.other.RData")

trainData <- BP.other[[1]]$data
getFeatures <- function(model.name){
  factor(ifelse(
    model.name=="major.class", 0,
    ifelse(
      model.name %in% c("ctree", "cforest", "glmnet", "xgboost"),
      n.train.features,
      ifelse(grepl("[(]", model.name), 2, 1))))
}

pred.dt <- BP.other.models$predictions
pred.dt[, n.features := getFeatures(model.name)]

roc.dt <- BP.other.models$roc[, {
  FPR.grid <- seq(0, 1, l=50)
  TPR.grid <- approx(FPR, TPR, FPR.grid)$y
  data.table(FPR=FPR.grid, TPR=TPR.grid)
}, by=.(train.name, test.name, test.fold, weight.name, model.name)]
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
set.means$min.auc <- min(pred.dt$auc)
set.means$max.auc <- max(pred.dt$auc)
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
  scale_fill_manual(values=c(balanced=NA, one="black"))+
  scale_color_manual(values=feature.colors)+
  geom_vline(aes(xintercept=100-mean.error),
             color="grey",
             data=best.error)+
  geom_point(aes(100-error.percent, model.fac, color=n.features, fill=weight.name),
             shape=21,
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

gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(train.name ~ test.name, labeller=label_both)+
  scale_color_manual(values=feature.colors)+
  scale_linetype_manual(values=c(balanced="dotted", one="solid"))+
  scale_fill_manual(values=c(balanced="white", one="black"))+
  geom_path(aes(
    FPR, TPR, color=n.features, linetype=weight.name,
    group=paste(test.fold, weight.name, model.name)),
            data=roc.dt)+
  geom_point(aes(
    FPR, TPR, color=n.features, fill=weight.name,
    group=paste(weight.name, model.name)),
             data=pred.dt)+
  coord_equal()

png("figure-BP-other-models-roc.png", width=900, height=750)
print(gg)
dev.off()

library(animint)
h <- 1200
viz <- list(
  auc=ggplot()+
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
             data=pred.dt)+
  theme_animint(width=850, height=h)+
  geom_segment(aes(min.auc, model.name,
                   xend=max.auc, yend=model.name,
                   clickSelects=model.name),
                size=9,
                alpha=0.2,
                color="black",
               data=set.means),
  roc=ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  theme_animint(width=850, height=h)+
  facet_grid(train.name ~ test.name, labeller=label_both)+
  scale_color_manual(values=feature.colors)+
  scale_linetype_manual(values=c(balanced="dotted", one="solid"))+
  scale_fill_manual(values=c(balanced="white", one="black"))+
  geom_path(aes(
    FPR, TPR, color=n.features, linetype=weight.name,
    showSelected=model.name,
    group=paste(test.fold, weight.name, model.name)),
            data=roc.dt)+
  geom_point(aes(
    FPR, TPR, color=n.features, fill=weight.name,
    showSelected=model.name),
             data=pred.dt),
  selector.types=list(model.name="multiple"),
  first=list(model.name=c("REVEL", "major.class", "glmnet", "VEST3_score")))
viz$auc

animint2dir(viz, "figure-BP-other-models")

