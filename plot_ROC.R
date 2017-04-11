works_with_R("3.3.3",
             ggplot2="2.1.0",
             data.table="1.10.4",
             "tdhock/namedCapture@05175927a45c301a18e8c6ebae67ea39a842d264",
             "tdhock/WeightedROC@4dc6e54e248dadd2c60f58b9870fd6b1d7e13eef")

arg.vec <- "data-minusHGMD/2016-minusHGMD.txt_predictions.RData"
arg.vec <- "data/mouse_test2.txt_predictions.RData"

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
wrong.label <- ifelse(is.pathogenic, "Benign", "Pathogenic")

pattern <- paste0(
  "(?<model>.*?)",
  "(?:",
  "_NA(?<na>keep|guess)",
  ")?",
  "[.]weights=",
  "(?<weight>[^_]+)")
score.name.vec <- names(keep)[keep]
match.mat <- str_match_named(score.name.vec, pattern, list(
  weight=function(x)ifelse(is.na(x)|x=="", "unknown", x)))
did.not.match <- is.na(match.mat[, "model"])
match.mat[did.not.match, "model"] <- score.name.vec[did.not.match]
match.mat[match.mat[,"na"]=="", "na"] <- "guess"

roc.list <- list()
auc.list <- list()
for(score.i in seq_along(score.name.vec)){
  score.name <- score.name.vec[[score.i]]
  pred.score <- input.dt[[score.name]]
  pred.label <- predictedScores$labels[[score.name]]
  meta <- match.mat[score.i,]
  AUC <- if(any(is.na(pred.score))){
    cat(score.name, "has missing values which will be treated as incorrect\n")
    pred.label <- ifelse(is.na(pred.label), wrong.label, pred.label)
    NA
  }else{
    tp.fp <- WeightedROC(pred.score, label.int)
    roc.list[[score.name]] <- data.table(
      score.name, meta, tp.fp)
    WeightedAUC(tp.fp)
  }
  pred.pathogenic <- pred.label=="Pathogenic"
  pred.benign <- pred.label=="Benign"
  TP <- as.double(sum(is.pathogenic & pred.pathogenic))
  FP <- as.double(sum(is.benign & pred.pathogenic))
  TN <- as.double(sum(is.benign & pred.benign))
  FN <- as.double(sum(is.pathogenic & pred.benign))
  n.positive <- sum(is.pathogenic)
  n.negative <- sum(is.benign)
  auc.list[[score.name]] <-
    data.table(
      score.name, meta, AUC,
      TPR=TP/n.positive,
      FPR=FP/n.negative,
      F1=2*TP/(2*TP+FP+FN),
      accuracy=(TP+TN)/(n.positive+n.negative),
      precision=TP/(TP+FP),
      specificity=TN/(TN+FP),
      sensitivity=TP/(TP+FN),
      MCC=(TP*TN-FP*FN)/(
        sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
}
roc <- do.call(rbind, roc.list)
auc <- do.call(rbind, auc.list)
auc.ord <- auc[, list(max.auc=max(AUC)), by=model][order(max.auc),]

auc.not.na <- auc[!is.na(AUC)]
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
  scale_fill_manual(values=c(guess="white", keep="grey"))+
  geom_point(aes(
    FPR, TPR,
    color=model),
             shape=21,
             data=auc.not.na)
out.png <- paste0(predictions.RData, ".png")
cat("Writing", out.png, "\n")
png(out.png)
print(gg)
dev.off()

if(require(animint)){
  auc.not.na[, model.fac := factor(model, auc.ord$model)]
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
                 data=auc.not.na),
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
                 data=auc.not.na)
    )
  figure.dir <- paste0(predictions.RData, "-figure")
  cat("Writing interactive figure in", figure.dir, "\n")
  animint2dir(viz, figure.dir)
}

out.csv <- paste0(predictions.RData, ".csv")
cat("Writing prediction statistics in", out.csv, "\n")
fwrite(auc, out.csv)
