figure-folds/index.html: figure-folds.R locus.fold.vec.RData
	R --no-save < $<
locus.fold.vec.RData: locus.fold.vec.R
	R --no-save < $<
figure-interactive/index.html: figure-interactive.R variants.RData
	R --no-save < $< 
variants.RData: variants.R
	R --no-save < $<
