works_with_R("3.2.1",
             data.table="1.9.5",
             "tdhock/ggplot2@a8b06ddb680acdcdbd927773b1011c562134e4d2",
             "tdhock/animint@f7d47b5a9b019fc3b12c0aa2d822572faf7ce4df")

load("variants.RData")

unique.or.multiple <- function(x)

HDR.LCR <- variants[["HDR/LCR"]]
getRegions <- function(region.type){
  v.diff <- diff(HDR.LCR==region.type)
  HDR.starts <- which(v.diff==1)+1
  HDR.ends <- which(v.diff==-1)
  data.table(chromStart=variants$POS[HDR.starts],
             chromEnd=variants$POS[HDR.ends],
             LOCUS_ID=variants$LOCUS_ID[HDR.ends],
             region.type)
}
regions <- rbind(getRegions("HDR"), getRegions("LCR"))

amp.coding.counts <- with(variants, table(LOCUS_ID, Coding))

unique.or.multiple <- function(x){
  u <- unique(paste(x))
  if(length(u) == 1){
    u
  }else{
    "multiple"
  }
}

amplicons <- variants[, .(chromStart=min(POS),
                          chromEnd=max(POS),
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
                      by=LOCUS_ID][, 
                          position := as.integer((chromStart+chromEnd)/2),
                        ]
amplicons[, highly.divergent.regions :=
            ifelse(region.type==".", "none", "some")]

chrom2int <- function(chrom){
  only.num <- sub("PyYM_([0-9]{2})_v1", "\\1", chrom)
  factor(as.integer(only.num), 1:14)
}

unfiltered.tp <- sum(variants$Validation == "TP")
unfiltered.fn <- sum(variants$Validation == "FN")

error.curves.list <- list()
amplicon.errors.list <- list()
QUAL.labels.list <- list()
for(QUAL.thresh in seq(3, 225, l=40)){
  QUAL.labels.list[[paste(QUAL.thresh)]] <-
    data.table(QUAL.thresh, chrom="PyYM_07_v1", position=2e6)
  variants.above <- variants[QUAL.thresh < QUAL, ]
  amplicon.status <- data.frame(amplicons)
  amplicon.status$filtered.tp <- 0
  amplicon.status$fp <- 0
  rownames(amplicon.status) <- amplicon.status$LOCUS_ID
  locus.fp.tp <-
    variants.above[, .(fp=sum(Validation=="FP"),
                       filtered.tp=sum(Validation=="TP")),
                   by=LOCUS_ID]
  amplicon.status[locus.fp.tp$LOCUS_ID, "filtered.tp"] <-
    locus.fp.tp$filtered.tp
  amplicon.status[locus.fp.tp$LOCUS_ID, "fp"] <-
    locus.fp.tp$fp
  amplicon.status$filtered.fn <-
    with(amplicon.status, unfiltered.tp - filtered.tp)
  amplicon.status$fn <-
    with(amplicon.status, unfiltered.fn + filtered.fn)
  amplicon.status$errors <- with(amplicon.status, fn + fp)
  amplicon.errors.list[[paste(QUAL.thresh)]] <-
    data.table(QUAL.thresh, amplicon.status)
  
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
    error.curves.list[[paste(QUAL.thresh, metric.name)]] <-
      data.table(QUAL.thresh, metric.name, metric.value, offset)
  }
}
amplicon.errors <- do.call(rbind, amplicon.errors.list)
error.curves <- do.call(rbind, error.curves.list)
QUAL.labels <- do.call(rbind, QUAL.labels.list)

thresh.30 <- with(error.curves, QUAL.thresh[which.min(abs(QUAL.thresh - 30))])

viz <-
  list(errors=ggplot()+
         ggtitle("error curves, select QUAL threshold")+
         xlab("QUAL threshold")+
         ylab("incorrectly called variants")+
         make_tallrect(error.curves, "QUAL.thresh")+
         geom_line(aes(QUAL.thresh, metric.value,
                       group=metric.name,
                       color=metric.name),
                   data=error.curves)+
         geom_text(aes(QUAL.thresh, metric.value+offset,
                       color=metric.name,
                       label=paste(metric.value, metric.name, " "),
                       showSelected=QUAL.thresh),
                   hjust=1,
                   data=error.curves),

       chroms=ggplot()+
         ggtitle("Sanger sequenced amplicons")+
         theme_animint(width=600)+
         geom_text(aes(chrom2int(chrom), position/1e3,
                       label=sprintf("QUAL threshold = %.1f", QUAL.thresh),
                       showSelected=QUAL.thresh),
                   data=QUAL.labels)+
         geom_text(aes(chrom2int(chrom), position/1e3,
                       label=paste(fp, "fp_"),
                       clickSelects=LOCUS_ID,
                       showSelected3=annotation,
                       showSelected2=highly.divergent.regions,
                       showSelected=QUAL.thresh),
                   hjust=1,
                   data=subset(amplicon.errors, fp != 0))+
         geom_text(aes(chrom2int(chrom), position/1e3,
                       label=paste0("_" , fn, " fn"),
                       clickSelects=LOCUS_ID,
                       showSelected3=annotation,
                       showSelected2=highly.divergent.regions,
                       showSelected=QUAL.thresh),
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

       first=list(QUAL.thresh=thresh.30),

       title="Malaria parasite NextGenSeq variant calling errors")

animint2dir(viz, "figure-interactive")

## BUG: metric.name and highly.divergent.regions legend entries do not
## fade to opacity: 0.5 after clicking.

##

##animint2gist(viz)
