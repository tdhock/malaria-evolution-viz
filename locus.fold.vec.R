works_with_R("3.2.2",
             kernlab="0.9.19",
             gbm="2.1",
             ada="2.0.3",
             randomForest="4.6.7",
             glmnet="1.9.5",
             data.table="1.9.6",
             caret="6.0.41",
             "tdhock/WeightedROC@1fdd55a7c8d94577071dd393843303d44100c6c1",
             "tdhock/ggplot2@a8b06ddb680acdcdbd927773b1011c562134e4d2",
             "tdhock/animint@87a18bb72b475044c1c4809dc6ebb19718bd9fd3")

library(parallel)
library(doParallel)
options(mc.cores=8)
registerDoParallel()

load("variants.RData")

variants$NegFQ <- -variants$FQ
variants.by.locus <- split(variants, variants$LOCUS_ID)
stop("compute nearest distance from other variant on the same locus")

names(variants.by.locus)
n.folds <- 3
set.seed(1)
locus.fold.vec <- sample(rep(1:n.folds, l=length(variants.by.locus)))
table(locus.fold.vec)

caret.methods <- c("knn", "rf", "svmRadial", "ada", "gbm", "glmnet")
all.methods <- c(caret.methods, "CVglmnet")
method.df <- expand.grid(method=all.methods, weight=c("one", "balanced"))
all.filterVars <- c(with(method.df, {
  paste0(method, ".", weight)
}), c("MQ", "QUAL", "NegFQ", "DP"))

all.roc.list <- list()
for(test.fold in 1:n.folds){
  cat(sprintf("%4d / %4d folds\n", test.fold, n.folds))
  is.test <- locus.fold.vec == test.fold
  train.validation <- do.call(rbind, variants.by.locus[!is.test])
  not.na <- subset(train.validation, !is.na(DP))
  not.na[, label := factor(Validation, c("FP", "TP"))]
  test.variants <- do.call(rbind, variants.by.locus[is.test])
  response.tab <- table(not.na$Validation)
  get.features <- function(dt){
    dt$is.CODING <- ifelse(dt$Coding=="CODING", 1, 0)
    dt$is.INTERGENIC <- ifelse(dt$Coding=="INTERGENIC", 1, 0)
    dt$is.INDEL <- ifelse(dt$Variant_type=="INDEL", 1, 0)
    dt$is.HDR <- dt[["HDR/LCR"]]=="HDR"
    dt$is.LCR <- dt[["HDR/LCR"]]=="LCR"
    col.names <-
      c("MQ", "QUAL", "FQ", "DP",
        "is.CODING", "is.INTERGENIC", "is.INDEL",
        "is.HDR", "is.LCR")
    as.matrix(dt[, col.names, with=FALSE])
  }
  train.validation.features <- get.features(not.na)
  test.features <- get.features(test.variants)
  test.is.na <- apply(is.na(test.features), 1, any)
  weight.list <-
    list(one=rep(1, nrow(not.na)),
         balanced=1/response.tab[paste(not.na$Validation)])
  sapply(weight.list, sum)
  for(weight.name in names(weight.list)){
    weight.vec <- weight.list[[weight.name]]
    glmnet.fit <-
      cv.glmnet(train.validation.features, not.na$label,
                weight.vec, family="binomial")
    not.na$prob <-
      predict(glmnet.fit, train.validation.features, type="response")
    not.na[order(prob),] #higher prob is more likely to be TP.
    test.variants[[paste0("CVglmnet.", weight.name)]] <-
      predict(glmnet.fit, test.features, type="response")
    for(method.name in caret.methods){
      caret.fit <-
        train(train.validation.features, not.na$label,
              method = method.name,
              weights=as.numeric(weight.vec),
              preProcess = c("center", "scale"),
              tuneLength = 10,
              trControl = trainControl(method = "cv"))
      caret.name <- paste0(method.name, ".", weight.name)
      test.variants[[caret.name]] <- NA
      test.variants[[caret.name]][!test.is.na] <-
        predict(caret.fit, test.features[!test.is.na, ], type="prob")[, "TP"]
    }
  }

  for(filterVar in all.filterVars){
    label.num <- test.variants[, ifelse(Validation=="TP", 1, -1)]
    score <- as.numeric(test.variants[[filterVar]])
    not.na <- data.table(label.num, score)[!is.na(score), ]
    roc <- not.na[, WeightedROC(score, label.num)]
    type <- ifelse(grepl("balanced", filterVar), "balanced",
                   ifelse(grepl("one", filterVar), "one", "univariate"))
    method <- sub("[.].*", "", filterVar)
    all.roc.list[[paste(test.fold, filterVar)]] <-
      data.frame(test.fold, filterVar, method, type, roc)
  }
}
all.roc <- do.call(rbind, all.roc.list)

save(variants.by.locus, locus.fold.vec, all.roc,
     file="locus.fold.vec.RData")
