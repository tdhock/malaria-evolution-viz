works_with_R("3.2.2", data.table="1.9.6")

original <- fread("N67toYMaug_22Feb13_reflist_2June14_eBAQ_MASTER_NOChr2_fullhomhet_annotated_for_Toby.txt")
names(original)[1] <- "CHROM"
n.positive <- original[Validation %in% c("FN", "TP"), .N]

variant.intervals <- fread("variant_intervals_for_Toby.txt")
variant.intervals[, bases := Stop - Start]
variant.intervals[, LOCUS_ID := factor(LOCUS_ID, LOCUS_ID)]

original[, LOCUS_ID := factor(LOCUS_ID, levels(variant.intervals$LOCUS_ID))]

unique.or.multiple <- function(x){
  u <- unique(paste(x))
  if(length(u) == 1){
    u
  }else{
    "multiple"
  }
}
amplicons <- original[, {
  .(firstVariant=min(POS),
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
    chrom=CHROM[1])
}, by=LOCUS_ID][, `:=`(
     position=as.integer((firstVariant+lastVariant)/2),
     highly.divergent.regions=
       ifelse(region.type==".", "none", "some")
  )]

HDR.LCR <- original[["HDR/LCR"]]
getRegions <- function(region.type){
  v.diff <- diff(HDR.LCR==region.type)
  HDR.starts <- which(v.diff==1)+1
  HDR.ends <- which(v.diff==-1)
  LOCUS_ID <- original$LOCUS_ID[HDR.ends]
  data.table(LOCUS_ID,
             regionStart=original$POS[HDR.starts],
             regionEnd=original$POS[HDR.ends],
             region.type)
}
regions <- rbind(getRegions("HDR"), getRegions("LCR"))

vcf.file.vec <- Sys.glob("Annotated_viz_VCFs/*.vcf")
vcf.list <- list()
count.tab.list <- list()
for(vcf.file in vcf.file.vec){
  one.vcf <- fread(vcf.file)
  if(ncol(one.vcf)==8){
    setnames(one.vcf, c(
      "CHROM", "POS",
      "Ref", "Alt",
      "Quality", "Depth", "MQ",
      "Validation"))
  }else{
    setnames(one.vcf, c(
      "CHROM", "POS",
      "Ref", "Alt",
      "DP", "QUAL", "MQ", "FQ",
      "Variant_type",
      "Validation"))
  }
  count.tab.list[[vcf.file]] <- table(one.vcf$Validation)
  stopifnot(nrow(original[Validation %in% c("FN", "TP"), .N]) == n.positive)
  vcf.list[[vcf.file]] <- one.vcf
}
(count.tab <- do.call(rbind, count.tab.list))

no.locus <- vcf.list[["Annotated_viz_VCFs/GATK_Diploid_annotated_viz.vcf"]]
no.locus[, POS1 := POS ]
setkey(no.locus, CHROM, POS, POS1)
only.limits <- amplicons[, {
  data.table(chrom, LOCUS_ID, firstVariant, lastVariant)
}]
setkey(only.limits, chrom, firstVariant, lastVariant)
variants <- foverlaps(no.locus, only.limits)
variants$Coding <- "unknown"
variants$Variant_type <- "unknown"

## Comment the next line to use one of the new data files defined in
## vcf.list above.
variants <- original

## "#CHROM - Chromosome name from the reference genome, in the format
## of PyYM_##_v1.  The naming convention is organism (genus + species,
## Plasmodium yoelii in this case, Py) then isolate (YM) and the
## version (1).  The hash is just there to easily grep out the header,
## etc.

## POS - Position of the variant in the reference, in the case of
## indels this is the start position.

## LOCUS_ID - a unique identifier for each amplicon (locus), named for
## the chromosome and then sequentially.

## Variant_type --binary variable, SNP or INDEL

## DP - Coverage depth

## QUAL - Phred-scaled quality score

## MQ - Mapping Quality

## FQ - Consensus Quality

## HDR/LCR - Whether the variant position is with in an HDR (highly
## divergent regions) or LCR (low-complexity region), otherwise (.)

## Validation - TP (true positive), FP (false positive), or FN (False
## Negative)

##  Error - Boolean (True for FN and FP, False for TP)

## Coding - Whether the variant is coding, intergenic, or in an
## intron.

## NA in quality scores (DP, QUAL, MQ, FQ) are in validation=FN only
## since even in the unfiltered NGS data set, there are some true
## variants (via Sanger seq) which are not detected at all.

YM_Chr <- fread("YM_Chr.bed.txt")
setnames(YM_Chr, c("chrom", "bases"))

save(YM_Chr, variants, variant.intervals,
     amplicons, regions,
     file="variants.RData")
