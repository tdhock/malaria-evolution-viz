works_with_R("3.2.2",
             data.table="1.9.6",
             "tdhock/ggplot2@a8b06ddb680acdcdbd927773b1011c562134e4d2",
             "tdhock/WeightedROC@da53b21f5eccaba513623e43326e5b8061d1c611",
             "tdhock/animint@d0b7dfcbb91b6f488f5ccfabc23fa9e91ca74708")

load("locus.fold.vec.RData")

fold.counts <- 
  data.table(TP=sapply(variants.by.locus, with, sum(Validation=="TP")),
             FP=sapply(variants.by.locus, with, sum(Validation=="FP")),
             fold=locus.fold.vec)
fold.counts[, list(TP=sum(TP), FP=sum(FP)), by=fold]

roc.by.filterVar.fold <- split(all.roc, all.roc[, c("filterVar", "test.fold")])
auc.by.filterVar.fold <- list()
for(filterVar.fold in names(roc.by.filterVar.fold)){
  roc <- roc.by.filterVar.fold[[filterVar.fold]]
  auc.by.filterVar.fold[[filterVar.fold]] <-
    data.table(roc[1,],
               metric.name=c("auc", "min.error"), 
               metric.value=c(WeightedAUC(roc), min(with(roc, FP+FN))))
}
all.auc <- do.call(rbind, auc.by.filterVar.fold)
method.ranks <- all.auc[metric.name=="auc", {
  .(method.mean=mean(metric.value))
}, by=method][order(method.mean, decreasing=TRUE),]
setkey(method.ranks, method)
filterVar.ranks <- all.auc[metric.name=="auc", {
  .(filterVar.mean=mean(metric.value))
}, by=.(filterVar, method)]
setkey(filterVar.ranks, method)
all.ranks <- filterVar.ranks[method.ranks][order(method.mean, filterVar.mean), ]
all.auc[, filterVar.fac := factor(filterVar, paste(all.ranks$filterVar))]

error.fun.list <- list(
  FN=function(df)df$FN,
  FP=function(df)df$FP,
  errors=function(df)with(df, FP+FN)
  )
all.error.list <- list()
for(error.type in names(error.fun.list)){
  error.fun <- error.fun.list[[error.type]]
  all.error.list[[error.type]] <-
    data.table(all.roc, error.type, error.value=error.fun(all.roc))
}
all.error <- do.call(rbind, all.error.list)

fp.fn.colors <- c(FP="skyblue",
                  fp="skyblue",
                  fn="#E41A1C",
                  FN="#E41A1C",
                  tn="white",
                  tp="grey",
                  errors="black")
fp.fn.linetypes <- c(errors="solid",
                     false.positive="solid",
                     false.negative="solid",
                     imprecision="dashed")
fp.fn.sizes <- c(errors=1,
                 false.positive=3,
                 false.negative=3,
                 imprecision=1)/1.2
min.lines <- all.error[error.type=="errors", {
  list(min.errors=min(error.value))
}, by=test.fold]

ggplot()+
  geom_hline(aes(yintercept=min.errors),
             data=min.lines,
             color="grey50")+
  theme_bw()+
  theme_animint(height=500)+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(test.fold ~ filterVar, labeller=function(var, val){
    if(var=="test.fold"){
      paste("test fold", val)
    }else{
      paste(val)
    }
  }, scales="free", space="fixed")+
  scale_color_manual(values=fp.fn.colors)+
  geom_line(aes(threshold, error.value,
                group=error.type, color=error.type),
            data=all.error)

viz <- list(
  auc=ggplot()+
    guides(color="none")+
    theme_bw()+
    theme_animint(height=500)+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(.~metric.name, scales="free", space="fixed")+
    scale_y_discrete("method . weights")+
    geom_point(aes(metric.value, filterVar.fac, color=method,
                   showSelected=method,
                   clickSelects=test.fold),
               size=5,
               data=all.auc),
  roc=ggplot()+
    scale_x_continuous("False positive rate",
                       breaks=c(0, 0.25, 0.5, 0.75, 1),
                       labels=c("0", "0.25", "0.5", "0.75", "1"))+
    coord_equal()+
  theme_bw()+
  theme_animint(width=500, height=500)+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(test.fold ~ type, labeller=function(var, val){
    if(var=="test.fold"){
      paste("test fold", val)
    }else{
      paste(val)
    }
  })+
  geom_path(aes(FPR, TPR, clickSelects=test.fold,
                group=method, tooltip=method, color=method),
            size=5,
            data=all.roc),
  selector.types=list(method="multiple"),
  error=ggplot()+
  geom_hline(aes(yintercept=min.errors,
                 showSelected=test.fold),
             data=min.lines,
             color="grey50")+
  theme_bw()+
  theme_animint(width=1500, height=500)+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(. ~ filterVar, labeller=function(var, val){
    if(var=="test.fold"){
      paste("test fold", val)
    }else{
      paste(val)
    }
  }, scales="free", space="fixed")+
  scale_color_manual(values=fp.fn.colors)+
  geom_line(aes(threshold, error.value,
                showSelected=test.fold,
                showSelected2=method,
                group=error.type, color=error.type),
            data=all.error),
  title="3-fold CV estimates variant calling test error"
)

animint2dir(viz, "figure-folds")
