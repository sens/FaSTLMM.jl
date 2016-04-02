PANDOC := pandoc -V geometry:margin=1in --variable mainfont=Georgia 
LATEXENGINE := --latex-engine=xelatex

%.html: %.Rmd
	Rscript -e 'library(knitr); library(markdown); rmarkdown::render("$<")'
clean:
	rm *.pdf
