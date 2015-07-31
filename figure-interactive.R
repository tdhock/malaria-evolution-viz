works_with_R("3.2.1",
             data.table="1.9.5",
             "tdhock/ggplot2@a8b06ddb680acdcdbd927773b1011c562134e4d2",
             "tdhock/animint@f7d47b5a9b019fc3b12c0aa2d822572faf7ce4df")

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
                          unfiltered.tp=sum(Validation=="TP"),
                          unfiltered.fn=sum(Validation=="FN"),
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

normalize <- function(LOCUS_ID, position){
  LID <- paste(LOCUS_ID)
  firstVariant <- amplicons[LID, ]$firstVariant
  lastVariant <- amplicons[LID, ]$lastVariant
  mid <- (firstVariant+lastVariant)/2
  left <- firstVariant
  left <- mid-200
  right <- mid+200
  bases <- right - left
  (position-left)/bases
}

ggplot()+
  geom_point(aes(POS, Validation),
             data=variants)+
  facet_grid(. ~ LOCUS_ID, scales="free")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))

ggplot()+
  geom_point(aes(normalize(LOCUS_ID, POS), Validation),
             data=variants)+
  facet_grid(. ~ LOCUS_ID)+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))

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
ggplot()+
  scale_fill_manual(values=fp.fn.colors)+
  scale_y_discrete(drop=FALSE)+
  geom_text(aes(normalize(LOCUS_ID, firstVariant), LOCUS_ID,
                label=paste0(firstVariant, "_")),
            hjust=1,
            data=amplicons)+
  geom_text(aes(normalize(LOCUS_ID, lastVariant), LOCUS_ID,
                label=paste0("_", lastVariant)),
            hjust=0,
            data=amplicons)+
  geom_segment(aes(normalize(LOCUS_ID, regionStart), LOCUS_ID,
                   xend=normalize(LOCUS_ID, regionEnd), yend=LOCUS_ID,
                   color=region.type),
               data=regions)+
  geom_point(aes(normalize(LOCUS_ID, POS), LOCUS_ID,
                 fill=Validation),
             color="black",
             pch=21,
             data=variants)

chrom2int <- function(chrom){
  only.num <- sub("PyYM_([0-9]{2})_v1", "\\1", chrom)
  factor(as.integer(only.num), 1:14)
}

unfiltered.tp <- sum(variants$Validation == "TP")
unfiltered.fn <- sum(variants$Validation == "FN")

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
  
  variants.above <- variants[filterVar.thresh < filterVar, ]
  amplicon.status <- data.frame(amplicons)
  amplicon.status$filtered.tp <- 0
  amplicon.status$fp <- 0
  rownames(amplicon.status) <- amplicon.status$LOCUS_ID
  locus.fp.tp <-
    variants.above[, .(fp=sum(Validation=="FP"),
                       filtered.tp=sum(Validation=="TP")),
                   by=LOCUS_ID]
  amplicon.status[paste(locus.fp.tp$LOCUS_ID), "filtered.tp"] <-
    locus.fp.tp$filtered.tp
  amplicon.status[paste(locus.fp.tp$LOCUS_ID), "fp"] <-
    locus.fp.tp$fp
  amplicon.status$filtered.fn <-
    with(amplicon.status, unfiltered.tp - filtered.tp)
  amplicon.status$fn <-
    with(amplicon.status, unfiltered.fn + filtered.fn)
  amplicon.status$errors <- with(amplicon.status, fn + fp)
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

viz <-
  list(errorCurves=ggplot()+
         ggtitle(paste("error curves, select",
                       filterVar, "threshold"))+
         xlab(paste(filterVar, "threshold"))+
         ylab("incorrectly called variants")+
         make_tallrect(error.curves, "filterVar.thresh")+
         geom_line(aes(filterVar.thresh, metric.value,
                       group=metric.name,
                       color=metric.name),
                   data=error.curves)+
         scale_color_manual(values=fp.fn.colors)+
         geom_text(aes(filterVar.thresh, metric.value+offset,
                       color=metric.name,
                       label=paste(metric.value, metric.name, " "),
                       showSelected=filterVar.thresh),
                   hjust=1,
                   data=error.curves),

       chroms=ggplot()+
         ggtitle("Sanger sequenced amplicons")+
         theme_animint(width=600)+
         geom_text(aes(chrom2int(chrom), position/1e3,
                       label=sprintf("%s threshold = %.1f",
                         filterVar, filterVar.thresh),
                       showSelected=filterVar.thresh),
                   data=filterVar.labels)+
         geom_text(aes(chrom2int(chrom), position/1e3,
                       label=paste(fp, "fp_"),
                       clickSelects=LOCUS_ID,
                       showSelected3=annotation,
                       showSelected2=highly.divergent.regions,
                       showSelected=filterVar.thresh),
                   hjust=1,
                   color=fp.fn.colors[["fp"]],
                   data=subset(amplicon.errors, fp != 0))+
         geom_text(aes(chrom2int(chrom), position/1e3,
                       label=paste0("_" , fn, " fn"),
                       clickSelects=LOCUS_ID,
                       showSelected3=annotation,
                       showSelected2=highly.divergent.regions,
                       showSelected=filterVar.thresh),
                   color=fp.fn.colors[["fn"]],
                   hjust=0,
                   data=subset(amplicon.errors, fn != 0))+
         geom_segment(aes(chrom2int(chrom), 0, 
                          yend=bases/1e3, xend=chrom2int(chrom)),
                      data=YM_Chr)+
         geom_point(aes(chrom2int(chrom), position/1e3,
                        color=highly.divergent.regions,
                        fill=annotation,
                        clickSelects=LOCUS_ID),
                    size=5,
                    data=amplicons)+
         scale_color_manual(values=c(none="white", some="black"))+
         scale_x_discrete("Malaria parasite yoelii yoelii chromosome",
                          drop=FALSE)+
         ylab("position on chromosome (kilo bases = kb)"),

       variants=ggplot()+
         ggtitle("Variants in each sanger sequenced amplicon")+
         theme_animint(width=1000, height=600)+
         scale_fill_manual(values=fp.fn.colors)+
         scale_y_discrete("amplicon LOCUS_ID", drop=FALSE)+
         scale_x_continuous("relative position on amplicon",
                            limits=c(-0.05, 1.05),
                            breaks=c())+
         geom_text(aes(normalize(LOCUS_ID, firstVariant), LOCUS_ID,
                       showSelected=highly.divergent.regions,
                       showSelected2=annotation,
                       label=paste0(firstVariant, "_")),
                   hjust=1,
                   data=amplicons)+
         geom_text(aes(normalize(LOCUS_ID, lastVariant), LOCUS_ID,
                       showSelected=highly.divergent.regions,
                       showSelected2=annotation,
                       label=paste0("_", lastVariant, " --- ",
                                    lastVariant-firstVariant, " bases")),
                   hjust=0,
                   data=amplicons)+
         geom_segment(aes(normalize(LOCUS_ID, firstVariant), LOCUS_ID,
                          xend=normalize(LOCUS_ID, lastVariant), yend=LOCUS_ID,
                          showSelected=highly.divergent.regions,
                          showSelected2=annotation,
                          clickSelects=LOCUS_ID),
                      size=12,
                      alpha=0.6,
                      data=amplicons)+
         geom_segment(aes(normalize(LOCUS_ID, regionStart), LOCUS_ID,
                          xend=normalize(LOCUS_ID, regionEnd), yend=LOCUS_ID,
                          showSelected=highly.divergent.regions,
                          showSelected2=annotation,
                          color=region.type),
                      size=8,
                      data=add.legend.vars(regions))+
         scale_color_manual(values=c("#E41A1C", #red
                              "#377EB8", #blue
                              "#4DAF4A", #green
                              "#984EA3", #purple
                              "#FF7F00", #orange
                              LCR="#FFFF33", #yellow
                              "#A65628",
                              "#F781BF",
                                     HDR="black"))+
         geom_point(aes(normalize(LOCUS_ID, POS), LOCUS_ID,
                        tooltip=paste(Coding, Variant_type),
                        showSelected=highly.divergent.regions,
                        showSelected2=annotation,
                        showSelected3=filterVar.thresh,
                        fill=error.type),
                    color="black",
                    pch=21,
                    size=4,
                    data=add.legend.vars(filtered.variants)),

       first=list(filterVar.thresh=best.thresh),

       title="Malaria parasite NextGenSeq variant calling errors")

animint2dir(viz, "figure-interactive")

## BUG: metric.name and highly.divergent.regions legend entries do not
## fade to opacity: 0.5 after clicking.

## BUG: HDR color legend shows fill not stroke!

## BUG: LCR region.type does not show up at first (but it does after
## clicking the region.type legend).

##animint2gist(viz)

