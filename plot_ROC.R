works_with_R("3.2.3",
             ggplot2="2.1.0",
             data.table="1.9.7",
             "tdhock/animint@7ae982edef45ba0de31e88fed89132b005e4a760",
             "tdhock/namedCapture@05175927a45c301a18e8c6ebae67ea39a842d264",
             "tdhock/WeightedROC@ef8f35ba7ae85e2995fa66afe13732cebb8b5633")

arg.vec <- "data/2016-minusHGMD.txt_predictedScores.RData"

arg.vec <- commandArgs(trailingOnly=TRUE)

stopifnot(length(arg.vec)==1)
predictions.RData <- normalizePath(arg.vec[1], mustWork=TRUE)
objs <- load(predictions.RData)

input.dt <- predictedScores$scores
if(! "clinvar" %in% names(input.dt)){
  stop("no clinvar column found in ", train.txt)
}
not.recognized <- input.dt[! clinvar %in% c("Pathogenic", "Benign"),]
if(nrow(not.recognized)){
  print(not.recognized)
  stop("clinvar should be either Pathogenic or Benign")
}
keep <- sapply(input.dt, is.numeric)
is.benign <- input.dt$clinvar=="Benign"
label.int <- ifelse(is.benign, -1, 1)
is.pathogenic <- input.dt$clinvar=="Pathogenic"

pattern <- paste0(
  "(?<model>.*?)",
  "[.]weights=",
  "(?<weight>[^_]+)")
score.name.vec <- names(keep)[keep]
match.mat <- str_match_named(score.name.vec, pattern, list(
  weight=function(x)ifelse(is.na(x)|x=="", "unknown", x)))
did.not.match <- is.na(match.mat[, "model"])
match.mat[did.not.match, "model"] <- score.name.vec[did.not.match]

roc.list <- list()
auc.list <- list()
for(score.i in seq_along(score.name.vec)){
  score.name <- score.name.vec[[score.i]]
  pred.score <- input.dt[[score.name]]
  if(any(is.na(pred.score))){
    cat("Skipping", score.name, "since it has missing values\n")
  }else{
    pred.label <- predictedScores$labels[[score.name]]
    tp.fp <- WeightedROC(pred.score, label.int)
    meta <- match.mat[score.i,]
    AUC <- WeightedAUC(tp.fp)
    auc.list[[score.name]] <- data.table(
      score.name, meta, AUC,
      TPR=sum(is.pathogenic & pred.label=="Pathogenic")/sum(is.pathogenic),
      FPR=sum(is.benign & pred.label=="Pathogenic")/sum(is.benign))
    roc.list[[score.name]] <- data.table(
      score.name, meta, tp.fp)
  }
}
roc <- do.call(rbind, roc.list)
auc <- do.call(rbind, auc.list)
auc.ord <- auc[, list(max.auc=max(AUC)), by=model][order(max.auc),]

gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ weight, labeller=label_both)+
  geom_abline(aes(
    slope=slope, intercept=intercept),
    data=data.table(slope=1, intercept=0),
    color="grey")+
  coord_equal()+
  geom_path(aes(
    FPR, TPR,
    group=paste(model, weight),
    color=model),
            data=roc)+
  geom_point(aes(
    FPR, TPR,
    color=model),
             shape=1,
             data=auc)
out.png <- paste0(predictions.RData, ".png")
cat("Writing", out.png, "\n")
png(out.png)
print(gg)
dev.off()

auc[, model.fac := factor(model, auc.ord$model)]
auc.limits <- data.table(
  max.AUC=max(auc$AUC),
  min.AUC=min(auc$AUC))
seg.dt <- data.table(auc.ord, auc.limits)
seg.dt[, model.fac := factor(model, auc.ord$model)]
viz <- list(
  auc=ggplot()+
    xlab("Area under the ROC curve")+
    ylab("model")+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    scale_fill_manual(values=c(one="white", balanced="black"))+
    geom_segment(aes(
      min.AUC, model.fac,
      xend=max.AUC, yend=model.fac,
      clickSelects=model),
                 data=seg.dt,
                 size=20,
                 alpha=0.1)+
    geom_vline(aes(xintercept=max.AUC),
               color="grey",
               data=auc.limits)+
    guides(color="none")+
    geom_point(aes(
      AUC, model.fac, fill=weight, color=model, clickSelects=model),
               data=auc),
  roc=ggplot()+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    geom_abline(aes(
      slope=slope, intercept=intercept),
                data=data.table(slope=1, intercept=0),
                color="grey")+
    coord_equal()+
    geom_path(aes(
      FPR, TPR,
      group=paste(model, weight),
      linetype=weight,
      clickSelects=model,
      showSelected=model,
      color=model),
              size=1,
              data=roc)+
    geom_point(aes(
      FPR, TPR,
      showSelected=weight,
      clickSelects=model,
      showSelected2=model,
      color=model),
               shape=1,
               data=auc)
  )
figure.dir <- paste0(predictions.RData, "-figure")
cat("Writing interactive figure in", figure.dir, "\n")
animint2dir(viz, figure.dir)
