# LaTeX Makefile

# Basic TeX file prefix
PROJ = elev

# R input path for figures and tables
RPATH = analysis/analysis.R

# Output paths for generated figures and tables
IPATH = manuscript/img
TPATH = manuscript/include

# Create paths of output .pdf files by changing suffix from .R to .pdf
# and prefix from `analysis` (RPATH) to `img` (IPATH)
# These files don't exist yet but the list of files in FIGS is needed 
# as a dependency for $(PROJ).pdf  

# Main 
all: manuscript/$(PROJ).pdf 

# Create pdf
manuscript/$(PROJ).pdf: manuscript/$(PROJ).tex $(FIGS) $(TABS)
	latexmk -pdf -quiet -bibtex $(PROJ).tex

# Create figures and tables
$(IPATH)/%.pdf: analysis.R
	Rscript $<

# Remove generated latex files and generated figures and tables
clean:
	latexmk -C 
	rm -f $(IPATH)/*.pdf
	rm -f $(TPATH)/*.tex
