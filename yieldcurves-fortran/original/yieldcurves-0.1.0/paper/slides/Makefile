# Makefile for yieldcurves conference slides
QUARTO      := RSTUDIO_PANDOC="/Applications/quarto/bin/tools" quarto
RSCRIPT     := Rscript

SLIDES_QMD  := slides.qmd
SLIDES_PDF  := slides.pdf
FIGURE_DIR  := figures
FIGURE_SCRIPT := make_slide_figures.R

.PHONY: all figures clean open

all: $(SLIDES_PDF)

$(SLIDES_PDF): $(SLIDES_QMD) header.tex $(FIGURE_DIR)/hero_figure.pdf $(FIGURE_DIR)/qrcode_paper.png
	$(QUARTO) render $(SLIDES_QMD) --to beamer

figures: $(FIGURE_DIR)/hero_figure.pdf $(FIGURE_DIR)/qrcode_paper.png

$(FIGURE_DIR)/hero_figure.pdf $(FIGURE_DIR)/qrcode_paper.png: $(FIGURE_SCRIPT)
	$(RSCRIPT) $(FIGURE_SCRIPT)

open: $(SLIDES_PDF)
	open $(SLIDES_PDF)

clean:
	rm -f slides.tex slides.log slides.aux slides.nav slides.out \
	      slides.snm slides.toc slides.vrb slides.bbl slides.blg
