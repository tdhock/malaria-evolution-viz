works_with_R("3.2.1", data.table="1.9.5")

unfiltered <- fread("N67toYMaug_22Feb13_reflist_2June14_eBAQ_MASTER_NOChr2_fullhomhet_annotated_for_Toby.txt")

## NA occurs in False negatives only since even in the unfiltered NGS
## data set, there are some true variants (via Sanger seq) which are
## not detected at all.


