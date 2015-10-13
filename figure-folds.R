works_with_R("3.2.2",
             data.table="1.9.6",
             "tdhock/ggplot2@a8b06ddb680acdcdbd927773b1011c562134e4d2",
             "tdhock/WeightedROC@3452d61638e16f547f73c1a0f3bf852a3751f29d",
             "tdhock/animint@61b8aed14e64a95fdbe9dfcd44c80254858305fd")

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
  min.i <- with(roc, which.min(FP+FN))
  best <- roc[min.i, ]
  auc.by.filterVar.fold[[filterVar.fold]] <-
    data.table(best,
               metric.name=c("auc", "min.error"), 
               metric.value=c(WeightedAUC(roc), with(best, FP+FN)))
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
add.filterVar <- function(df, levs){
  df$filterVar.fac <- factor(df$filterVar, levs)
  df
}
add.filterVar.fac <- function(df){
  add.filterVar(df, rev(paste(all.ranks$filterVar)))
}
add.filterVar.rev <- function(df){
  add.filterVar(df, paste(all.ranks$filterVar))
}
all.roc$method <- factor(all.roc$method, rev(unique(paste(all.ranks$method))))
all.roc$thresh.type <- "selected"
all.roc$metric.name <- "min.error"

all.auc[, thresh.type := "min error"]

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

tallrect.dt <- all.error[error.type=="errors", {
  vals <- threshold
  only.finite <- vals[is.finite(vals)]
  delta.finite <- diff(only.finite)/2
  vals[vals==Inf] <- max(only.finite)+mean(delta.finite)
  Delta <- diff(vals)/2
  breaks <- c(vals[1] - Delta[1], vals[-1] - Delta, vals[length(vals)] + 
        Delta[length(Delta)])
  stopifnot(length(breaks) == length(vals) + 1)
  data.table(threshold, thresh.type="selected",
             xmin = breaks[-length(breaks)], xmax = breaks[-1])
}, by=.(filterVar, method, test.fold)]
tallrect.dt[filterVar=="glmnet.balanced" & test.fold==3,]  
all.error[filterVar=="glmnet.balanced" & test.fold==3 & error.type=="errors",]

min.lines <- all.error[error.type=="errors", {
  list(min.errors=min(error.value))
}, by=test.fold]
largest.min.error <- all.auc[metric.name=="min.error", max(metric.value)]
all.roc$error.or.Inf <- with(all.roc, {
  ifelse(largest.min.error < FP+FN, Inf, FP+FN)
})

## Save for creating an animint test.
## VariantModels <- list(
##   roc=data.frame(all.roc),
##   auc=data.frame(all.auc),
##   error=data.frame(all.error),
##   ranks=data.frame(all.ranks),
##   thresholds=data.frame(tallrect.dt),
##   minima=data.frame(min.lines))
## save(VariantModels, file="VariantModels.RData")

thresh.colors <- c("min error"="black", selected="white")
method.colors <- 
  c(knn="#8DD3C7", #green
    "#FFFFB3", #yellow
    svmRadial="#BEBADA", #pale violet
    ada="#FB8072", #pink-orange
    gbm="#FB8072", #pink-orange
    glmnet="#80B1D3", #blue
    glmnetBinDev="#80B1D3", #blue
    glmnetAcc="#80B1D3", #blue
    MQ="#FDB462", #orange
    QUAL="#B3DE69", #green
    NegFQ="#FCCDE5", #pink-violet
    DP="#D9D9D9", #grey
    rf="#BC80BD", #purple
    "#CCEBC5", #greenish yellow
    "#FFED6F") #gold
fp.fn.colors <- c(FP="skyblue",
                  fp="skyblue",
                  fn="#E41A1C",
                  FN="#E41A1C",
                  tn="white",
                  tp="grey",
                  errors="black")

viz <- list(
  auc=ggplot()+
    ggtitle("Performance on 3 test folds")+
    theme_bw()+
    theme_animint(height=500)+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(.~metric.name, scales="free", space="fixed")+
    scale_y_discrete("method . weights")+
    scale_x_continuous("")+
    scale_color_manual(values=method.colors, guide="none")+
    scale_fill_manual("threshold", values=thresh.colors, guide="none")+
    geom_point(aes(metric.value, filterVar.fac, color=method,
                   fill=thresh.type,
                   showSelected=method,
                   showSelected2=thresh.type,
                   clickSelects=test.fold),
               size=5,
               pch=21,
               data=add.filterVar.rev(all.auc))+
    geom_point(aes(
      error.or.Inf,
      filterVar.fac, 
      showSelected=test.fold,
      key=filterVar,
      showSelected2=thresh.type,
      showSelected3=method,
      showSelected.variable=paste0(filterVar, "_fold", test.fold),
      showSelected.value=threshold,
      fill=thresh.type, color=method),
               size=4,
               pch=21,
               data=add.filterVar.rev(all.roc)),
  roc=ggplot()+
    ggtitle("ROC curves by weights and test fold")+
    scale_y_continuous("True positive rate")+
    scale_x_continuous("False positive rate",
                       breaks=c(0, 0.25, 0.5, 0.75, 1),
                       labels=c("0", "0.25", "0.5", "0.75", "1"))+
    scale_color_manual(values=method.colors)+
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
                  showSelected=method,
                  group=method, tooltip=method, color=method),
              size=5,
              data=all.roc)+
    scale_fill_manual("threshold", values=thresh.colors)+
    geom_point(aes(FPR, TPR, color=method,
                   showSelected=method,
                   clickSelects=test.fold,
                   fill=thresh.type),
               pch=21,
               size=4,
               data=all.auc)+
    geom_point(aes(
      FPR, TPR, clickSelects=test.fold,
      key=method,
      showSelected.variable=paste0(filterVar, "_fold", test.fold),
      showSelected.value=threshold,
      showSelected=method,
      fill=thresh.type, color=method),
               size=3,
               pch=21,
               data=all.roc),
  error=ggplot()+
    geom_hline(aes(yintercept=min.errors,
                   showSelected=test.fold),
               data=min.lines,
               color="grey50")+
    geom_vline(aes(xintercept=threshold,
                   showSelected2=method,
                   showSelected3=thresh.type,
                   showSelected=test.fold),
               data=add.filterVar.fac(all.auc[metric.name=="min.error",]),
               color="grey50")+
    theme_bw()+
    theme_animint(width=1800, height=500)+
    theme(panel.margin=grid::unit(0, "cm"))+
    theme(axis.text.x=element_text(angle=90))+
    facet_grid(. ~ filterVar.fac, labeller=function(var, val){
      sub("balanced", "b", sub("one", "1", val))
    }, scales="free", space="fixed")+
    scale_color_manual(values=fp.fn.colors)+
    geom_line(aes(threshold, error.value,
                  showSelected=test.fold,
                  showSelected2=method,
                  showSelected3=thresh.type,
                  group=error.type, color=error.type),
              data=add.filterVar.fac(all.error))+
    scale_fill_manual(values=method.colors, guide="none")+
    geom_tallrect(aes(
      xmin=xmin, xmax=xmax,
      showSelected=test.fold,
      showSelected2=method,
      showSelected3=thresh.type,
      clickSelects.variable=paste0(filterVar.fac, "_fold", test.fold),
      clickSelects.value=threshold,
      fill=method),
                  alpha=0.5,
                  color=NA,
                  data=add.filterVar.fac(tallrect.dt)),
  selector.types=list(method="multiple", thresh.type="multiple"),
  title="3-fold CV estimates variant calling test error",
  first=with(all.auc[metric.name=="min.error", ], {
    structure(as.list(threshold), names=paste0(filterVar, "_fold", test.fold))
  }),
  duration=with(all.auc[metric.name=="min.error", ], {
    structure(as.list(rep(2000, length(threshold))),
              names=paste0(filterVar, "_fold", test.fold))
  })
)

viz$error+
  facet_grid(test.fold ~ filterVar.fac, labeller=function(var, val){
    if(var=="test.fold"){
      paste("test fold", val)
    }else{
      paste(val)
    }
  }, scales="free", space="fixed")

animint2dir(viz, "figure-folds")
