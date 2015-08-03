works_with_R("3.2.1",
             data.table="1.9.5")

load("variants.RData")

filterVar <- "DP"
filterVar <- "QUAL"
filterVar <- "FQ"
filterVar <- "MQ"

ptab <- table(variants$POS)
duplicated.pos <- ptab[1 < ptab]
variants[POS %in% as.integer(names(duplicated.pos))]

variants[, `:=`(pos.fac=factor(POS, POS),
                LOCUS_ID=factor(LOCUS_ID, unique(LOCUS_ID)))]
variants$filterVar <- variants[[filterVar]]

filterVar.thresh.vec <-
  seq(min(variants$filterVar, na.rm=TRUE),
      max(variants$filterVar, na.rm=TRUE),
      l=40)

amp.coding.counts <- with(variants, table(LOCUS_ID, Coding))

unique.or.multiple <- function(x){
  u <- unique(paste(x))
  if(length(u) == 1){
    u
  }else{
    "multiple"
  }
}

amplicons <- variants[, .(firstVariant=min(POS),
                          lastVariant=max(POS),
                          annotation=unique.or.multiple(Coding),
                          region.type={
                            if(all(`HDR/LCR` == ".")){
                              "."
                            }else{
                              not.dot <- `HDR/LCR`[`HDR/LCR` != "."]
                              ##unique.or.multiple(not.dot)
                              u <- sort(unique(not.dot))
                              paste(u, collapse="/")
                            }
                          },
                          chrom=CHROM[1]),
                      by=LOCUS_ID][, `:=`(
                        position=as.integer((firstVariant+lastVariant)/2),
                        highly.divergent.regions=
                          ifelse(region.type==".", "none", "some")
                        )]
setkey(amplicons, LOCUS_ID)

HDR.LCR <- variants[["HDR/LCR"]]
getRegions <- function(region.type){
  v.diff <- diff(HDR.LCR==region.type)
  HDR.starts <- which(v.diff==1)+1
  HDR.ends <- which(v.diff==-1)
  LOCUS_ID <- variants$LOCUS_ID[HDR.ends]
  data.table(LOCUS_ID,
             regionStart=variants$POS[HDR.starts],
             regionEnd=variants$POS[HDR.ends],
             region.type)
}
regions <- rbind(getRegions("HDR"), getRegions("LCR"))

error.curves.list <- list()
amplicon.errors.list <- list()
filterVar.labels.list <- list()
filtered.variants.list <- list()
for(filterVar.thresh in filterVar.thresh.vec){
  filterVar.labels.list[[paste(filterVar.thresh)]] <-
    data.table(filterVar.thresh, chrom="PyYM_07_v1", position=2e6)
  
  labeled.variants <- data.table(filterVar.thresh, variants)
  labeled.variants[, call:=ifelse(filterVar.thresh < filterVar,
                       "variant", "filtered")]
  labeled.variants[, error.type:=ifelse(is.na(call), "fn",
                       ifelse(call=="variant",
                              ifelse(Validation=="TP", "tp", "fp"),
                              ifelse(Validation=="TP", "fn", "tn")))]
  filtered.variants.list[[paste(filterVar.thresh)]] <- labeled.variants

  amplicon.fp.fn <- 
    labeled.variants[, .(fp=sum(error.type=="fp"),
                         fn=sum(error.type=="fn")),
                     by=LOCUS_ID]
  amplicon.status <- amplicons[amplicon.fp.fn, ]
  amplicon.errors.list[[paste(filterVar.thresh)]] <-
    data.table(filterVar.thresh, amplicon.status)
  
  fp <- with(amplicon.status, sum(fp))
  fn <- with(amplicon.status, sum(fn))
  errors <- fp + fn
  fn.offset <- ifelse(fp < fn, 20, -20)
  fp.offset <- -fn.offset
  for(metric.name in c("fp", "fn", "errors")){
    metric.value <- get(metric.name)
    offset <-
      ifelse(metric.name=="errors", 10,
             ifelse(metric.name=="fp", fp.offset, fn.offset))
    error.curves.list[[paste(filterVar.thresh, metric.name)]] <-
      data.table(filterVar.thresh, metric.name, metric.value, offset)
  }
}
filtered.variants <- do.call(rbind, filtered.variants.list)
amplicon.errors <- do.call(rbind, amplicon.errors.list)
error.curves <- do.call(rbind, error.curves.list)
filterVar.labels <- do.call(rbind, filterVar.labels.list)

only.errors <- error.curves[metric.name=="errors", ]
best.thresh <-
  only.errors[which.min(metric.value), ]$filterVar.thresh

add.legend.vars <- function(dt){
  LID <- paste(dt$LOCUS_ID)
  legend.vars <- amplicons[LID, .(highly.divergent.regions, annotation)]
  data.table(dt, legend.vars)
}

malaria.dt.list <-
  list(error.variants=add.legend.vars(filtered.variants),
       regions=add.legend.vars(regions),
       amplicons=amplicons,
       chroms=YM_Chr,
       error.amplicons=amplicon.errors,
       filterVar.labels=filterVar.labels,
       error.curves=error.curves,
       filterVar=data.table(filterVar, best.thresh))

malaria <- list()
for(data.name in names(malaria.dt.list)){
  malaria[[data.name]] <- data.frame(malaria.dt.list[[data.name]])
}

save(malaria, file="malaria.RData")

