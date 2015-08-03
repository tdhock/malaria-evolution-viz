figure-interactive/index.html: figure-interactive.R malaria.RData
	R --no-save < $< 
malaria.RData: malaria.R variants.RData
	R --no-save < $<
variants.RData: variants.R
	R --no-save < $<
