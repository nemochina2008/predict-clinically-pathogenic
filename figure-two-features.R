works_with_R("3.2.3", data.table="1.9.7", ggplot2="1.0.1", party="1.0.13",
             scales="0.2.5")

load("trainData.RData")
load("folds.RData")

train.dt <- with(trainData, data.table(clinvar=factor(label.vec), feature.mat))
##train.dt[, MutationTaster := ifelse(is.na(`MutationTaster-score`), -0.5, `MutationTaster-score`)]
train.dt[, VEST3 := ifelse(is.na(VEST3_score), -0.5, VEST3_score)]

ggplot()+
  scale_y_continuous(breaks=c(-0.5, 0, 0.5, 1),
                     labels=c("NA", 0, 0.5, 1))+
  geom_point(aes(CADD_phred, VEST3, color=clinvar),
             shape=1,
             data=train.dt)

n.grid <- 50
grid.dt <- data.table(expand.grid(
  `CADD-score`=seq(0, max(train.dt[, `CADD-score`]), l=n.grid),
  `SIFT-score-(P.damaging)`=seq(0, 1, l=n.grid)))
grid.contour.list <- list()
for(fun.name in c("ctree", "cforest")){
  fun <- get(fun.name)
  fit <- fun(clinvar ~ `CADD-score` + `SIFT-score-(P.damaging)`,
             data=train.dt)
  pred.list <- treeresponse(fit, grid.dt)
  is.pathogenic <- levels(train.dt$clinvar)=="Pathogenic"
  grid.dt$prob.Pathogenic <- sapply(pred.list, function(df.or.vec){
    unlist(df.or.vec)[is.pathogenic]
  })
  grid.contour.list[[fun.name]] <- data.table(fun.name, grid.dt)
}
grid.contour <- do.call(rbind, grid.contour.list)

gg <- ggplot()+
  geom_tile(aes(`CADD-score`, `SIFT-score-(P.damaging)`,
                fill=prob.Pathogenic),
            shape=".",
            data=grid.contour)+
  geom_point(aes(`CADD-score`, `SIFT-score-(P.damaging)`),
             size=3,
             color="black",
             data=train.dt)+
  geom_point(aes(`CADD-score`, `SIFT-score-(P.damaging)`, color=clinvar),
             size=2,
             data=train.dt)+
  scale_color_manual(values=c(Benign=muted("blue"), Pathogenic=muted("red")),
                     breaks=c("Pathogenic", "Benign"))+
  scale_fill_gradient2("probability\nPathogenic", midpoint=0.5,
                       low=muted("blue"), high=muted("red"),
                       limits=c(0, 1))+
  ## scale_fill_manual(values=c(Benign="white", Pathogenic="black"))+
  ## scale_color_gradient2("probability\nPathogenic", midpoint=0.5,
  ##                       limits=c(0, 1))+
  ## geom_tile(aes(`CADD-score`, `SIFT-score-(P.damaging)`,
  ##               color=prob.Pathogenic),
  ##           shape=".",
  ##           data=grid.contour)+
  ## geom_point(aes(`CADD-score`, `SIFT-score-(P.damaging)`, fill=clinvar),
  ##            shape=21,
  ##            data=train.dt)+
  ## geom_point(aes(`CADD-score`, `SIFT-score-(P.damaging)`,
  ##                color=prob.Pathogenic),
  ##            shape=".",
  ##            data=grid.contour)+
  ## geom_contour(aes(`CADD-score`, `SIFT-score-(P.damaging)`,
  ##                  z=cforest, color=..level..),
  ##              data=grid.dt)+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ fun.name, labeller=function(var, val){
    paste(val, "model")
  })

png("figure-two-features.png", width=700)
print(gg)
dev.off()
