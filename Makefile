figure-interactive/index.html: figure-interactive.R variants.RData
	R --no-save < $< 
variants.RData: variants.R
	R --no-save < $<
