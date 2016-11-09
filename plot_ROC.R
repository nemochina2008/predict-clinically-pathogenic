works_with_R("3.2.3",
             ggplot2="2.1.0",
             data.table="1.9.7",
             "tdhock/namedCapture@05175927a45c301a18e8c6ebae67ea39a842d264",
             "tdhock/WeightedROC@ef8f35ba7ae85e2995fa66afe13732cebb8b5633")

arg.vec <- c(
  "data-2016-11-09/test1.txt_predictions.txt"
  )

arg.vec <- commandArgs(trailingOnly=TRUE)

stopifnot(length(arg.vec)==1)
predictions.txt <- normalizePath(arg.vec[1], mustWork=TRUE)
input.dt <- fread(predictions.txt, na.strings=c("-", "NA"))

if(! "clinvar" %in% names(input.dt)){
  stop("no clinvar column found in ", train.txt)
}
not.recognized <- input.dt[! clinvar %in% c("Pathogenic", "Benign"),]
if(nrow(not.recognized)){
  print(not.recognized)
  stop("clinvar should be either Pathogenic or Benign")
}
keep <- sapply(input.dt, is.numeric)
label.int <- ifelse(input.dt$clinvar=="Benign", -1, 1)

pattern <- paste0(
  "(?<model>.*?)",
  "(?:",
  "[.]weights=",
  "(?<weight>[^_]+)",
  ")?",
  "(?:_score)?$")
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
    tp.fp <- WeightedROC(pred.score, label.int)
    meta <- match.mat[score.i,]
    AUC <- WeightedAUC(tp.fp)
    auc.list[[score.name]] <- data.table(
      score.name, meta, AUC)
    roc.list[[score.name]] <- data.table(
      score.name, meta, tp.fp)
  }
}
roc <- do.call(rbind, roc.list)
auc <- do.call(rbind, auc.list)
print(auc[order(AUC),])

gg <- ggplot()+
  theme_bw()+
  geom_abline(aes(
    slope=slope, intercept=intercept),
    data=data.table(slope=1, intercept=0),
    color="grey")+
  coord_equal()+
  geom_path(aes(
    FPR, TPR,
    group=paste(model, weight),
    color=model,
    linetype=weight),
    data=roc)

png(paste0(predictions.txt, ".png"))
print(gg)
dev.off()
