works_with_R("3.2.2",
             glmnet="1.9.5",
             data.table="1.9.6",
             "tdhock/ggplot2@a8b06ddb680acdcdbd927773b1011c562134e4d2",
             "tdhock/animint@0f0866bdab11c73dba8eba27b7e32775de85b040")

load("variants.RData")

variants.by.locus <- split(variants, variants$LOCUS_ID)
names(variants.by.locus)
n.folds <- 4
set.seed(1)
locus.fold.vec <- sample(rep(1:n.folds, l=length(variants.by.locus)))
table(locus.fold.vec)

for(test.fold in 1:n.folds){
  is.test <- locus.fold.vec == test.fold
  train.validation <- do.call(rbind, variants.by.locus[!is.test])
  not.na <- subset(train.validation, !is.na(DP))
  not.na[, label := factor(Validation, c("FP", "TP"))]
  test.variants <- do.call(rbind, variants.by.locus[is.test])
  response.tab <- table(not.na$Validation)
  get.features <- function(dt){
    as.matrix(dt[, c("MQ", "QUAL", "FQ", "DP"), with=FALSE])
  }
  train.validation.features <- get.features(not.na)
  test.features <- get.features(test.variants)
  weight.list <-
    list(one=rep(1, nrow(not.na)),
         balanced=1/response.tab[paste(not.na$Validation)])
  sapply(weight.list, sum)
  for(weight.name in names(weight.list)){
    weight.vec <- weight.list[[weight.name]]
    fit <-
      cv.glmnet(train.validation.features, not.na$label,
                weight.vec, family="binomial")
    not.na$prob <- predict(fit, train.validation.features, type="response")
    not.na[order(prob),] #higher prob is more likely to be TP.
    test.variants[[weight.name]] <- predict(fit, test.features, type="response")
  }
  for(filterVar in c("MQ", "QUAL", "FQ", "DP", "one", "balanced")){
    stop("TODO: ROC curve and dot on the test set")
  }
}
