
BUILDDIR=build

STAT_FILES=$(shell ls ../stats/*)

all: $(BUILDDIR) gen_stats_tables thesis_draft

thesis_draft: gen_stats_tables
	cp draft_thesis_half/* $(BUILDDIR)
	cd $(BUILDDIR); latexmk -pdf main.tex
	cd $(BUILDDIR); make -j8 -f main.makefile
	cd $(BUILDDIR); latex -halt-on-error main.tex
	cd $(BUILDDIR); bibtex main.aux
	cd $(BUILDDIR); pdflatex -halt-on-error main.tex
	cd $(BUILDDIR); pdflatex -halt-on-error main.tex
	cd $(BUILDDIR); make -j8 -f main.makefile
	cd $(BUILDDIR); pdflatex -halt-on-error main.tex
	acroread $(BUILDDIR)/main.pdf

gen_stats_tables: scripts/gen_stats_tables.sh $(STAT_FILES)
	@$<

$(BUILDDIR):
	@mkdir -p $@

clean:
	rm -rf $(BUILDDIR)
