works_with_R("3.2.1", data.table="1.9.5")

unfiltered <- fread("N67toYMaug_22Feb13_reflist_2June14_eBAQ_MASTER_NOChr2_fullhomhet_annotated_for_Toby.txt")

variants <- unfiltered

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

save(YM_Chr, variants, file="variants.RData")
