# docs/Makefile — render markdown docs to stakeholder-styled PDFs.
#
# Uses the GML / NOAA report template (Times serif body, navy headings,
# GML+NOAA header band) defined in bin/md2pdf.py. Dependencies are
# auto-installed by uv on first run from the PEP 723 metadata inside the
# script — `brew install uv` if you don't have it yet.
#
# Quick usage:
#   make -C docs                          # list .md docs and their PDF targets
#   make -C docs all                      # rebuild every PDF
#   make -C docs co_inversion_proposal.pdf   # build a single doc
#
# From the repo root the same targets are also reachable as:
#   make -f docs/Makefile -C docs <target>
#
# Header logos: docs/assets/GML-logo.jpg + docs/assets/noaa-logo.svg.

UV          ?= uv
RENDER      := $(UV) run --quiet --script bin/md2pdf.py
MD_SOURCES  := $(wildcard *.md)
PDF_TARGETS := $(MD_SOURCES:.md=.pdf)

.DEFAULT_GOAL := help

.PHONY: help
help:
	@echo "make targets:"
	@echo "  all              build every PDF in docs/"
	@echo "  <name>.pdf       build a single doc (e.g. make co_inversion_proposal.pdf)"
	@echo "  clean            remove all generated PDFs"
	@echo ""
	@echo "PDFs we know how to build:"
	@printf "  %s\n" $(PDF_TARGETS)

.PHONY: all
all: $(PDF_TARGETS)

# Pattern rule: foo.md -> foo.pdf via the shared template renderer.
%.pdf: %.md bin/md2pdf.py assets/GML-logo.jpg assets/noaa-logo.svg
	$(RENDER) $< $@

.PHONY: clean
clean:
	rm -f $(PDF_TARGETS) assets/.noaa-logo-hires.png
