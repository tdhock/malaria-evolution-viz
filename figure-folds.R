works_with_R("3.2.2",
             data.table="1.9.6",
             "tdhock/ggplot2@a8b06ddb680acdcdbd927773b1011c562134e4d2",
             "tdhock/animint@87a18bb72b475044c1c4809dc6ebb19718bd9fd3")

load("locus.fold.vec.RData")

fold.counts <- 
  data.table(TP=sapply(variants.by.locus, with, sum(Validation=="TP")),
             FN=sapply(variants.by.locus, with, sum(Validation=="FN")),
             fold=locus.fold.vec)
fold.counts[, list(variants=sum(variants)), by=fold]

viz <- list(
  roc=ggplot()+
    coord_equal()+
  theme_bw()+
  theme_animint(width=1000, height=1000)+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(test.fold ~ type)+
  geom_path(aes(FPR, TPR, group=method, tooltip=method, color=method),
            size=5,
            data=all.roc)
)

animint2dir(viz, "figure-folds")
