#!/bin/bash

# Build the paper with analyses, images, tables, references etc.

# Run analyses
Rscript analysis/analysis.R

# Fix table snippets

## Remove $ around numbers
sed -i 's/\$//g' manuscript/include/*.tex

## Align numerical columns on decimals
sed -i 's/ccc/lS[table-format=3.2]S[table-format=3.2]/g' manuscript/include/comp_radius.tex

sed -i 's/cccccccc/lS[table-format=3.2]S[table-format=3.2]S[table-format=3.2]S[table-format=3.2]S[table-format=3.2]S[table-format=3.2]S[table-format=3.2]/g' manuscript/include/*_model_comparison.tex

sed -i 's/ccccccccc/lS[table-format=3.2]S[table-format=3.2]S[table-format=3.2]S[table-format=3.2]S[table-format=3.2]S[table-format=3.2]S[table-format=3.2]S[table-format=3.2]/g' manuscript/include/site_char.tex

sed -i 's/ccccc/llr@{\\hspace{0.2\\tabcolsep}}lr@{\\hspace{0.2\\tabcolsep}}lr@{\\hspace{0.2\\tabcolsep}}l/g' manuscript/include/species_elevcode_tally.tex

## Put column headers in curly braces
sed -i '10s/.*/{Site} \& {Trees ha\\textsuperscript{-1}} \& {$C_R$ (m)} \\\\/' manuscript/include/comp_radius.tex

sed -i '10s/.*/{Fixed eff.} \& {AIC} \& {$\Delta{}AIC_r$} \& {$W_i$} \& {$R^2_c$} \& {$R^2_m$} \% {Slope} \& {SE} \\\\/' manuscript/include/*_model_comparison.tex

sed -i '10s/.*/{Site code} \& {Elev.} \& {Precip. (mm year\\textsuperscript{-1})} \& {Annual temp. (\\textdegree{}C)} \& {Slope (\\textdegree)} \& {Soil C (\\%)} \& {Soil N (\\%)} \& {Soil pH} \& {Trees ha\\textsuperscript{-1}} \\\\/' manuscript/include/site_char.tex

sed -i '10s/.*/{Species code} \& {Species} \& \\multicolumn{2}{c}{Bottom} \& \\multicolumn{2}{c}{Middle} \& \\multicolumn{2}{c}{Top} \\\\/' manuscript/include/species_elevcode_tally.tex

## Format column positions by equals sign
sed -i 's/=/ \& =/g' manuscript/include/species_elevcode_tally.tex
sed -i 's/NA /NA \& /g' manuscript/include/species_elevcode_tally.tex

## Put NA in empty cells in site characteristics
sed -i '12,21s/  / NA /g' manuscript/include/site_char.tex

## Adjust column sep 
sed -i 's/extracolsep{5pt}/extracolsep{0pt}/' manuscript/include/site_char.tex

## Format fixed effects in *_model_comparison.tex
sed -i 's/elev\\_scale.comp\\_adult\\_metric\\_log\\_scale.lai\\_scale \&/Elev. + ISI + LAI \&/' manuscript/include/*_model_comparison.tex

sed -i 's/elev\\_scale.comp\\_adult\\_metric\\_log\\_scale \&/Elev. + ISI \&/' manuscript/include/*_model_comparison.tex

sed -i 's/elev\\_scale.lai\\_scale \&/Elev. + LAI \&/' manuscript/include/*_model_comparison.tex

sed -i 's/comp\\_adult\\_metric\\_log\\_scale.lai\\_scale \&/ISI + LAI \&/' manuscript/include/*_model_comparison.tex

sed -i 's/elev\\_scale \&/Elev. \&/' manuscript/include/*_model_comparison.tex

sed -i 's/lai\\_scale \&/LAI \&/' manuscript/include/*_model_comparison.tex

sed -i 's/comp\\_adult\\_metric\\_log\\_scale \&/ISI \&/' manuscript/include/*_model_comparison.tex

## Add captions
sed -i 's/caption{}/caption{Competition radius used for adult competition measurements for each site based on the number of trees per hectare.}/g' manuscript/include/comp_radius.tex

sed -i "s/caption{}/caption{Site environmental characteristics for each 1 ha plot sampled. NA indicates that no data was available. Adapted from \\citet{Whitaker2014}.}/g" manuscript/include/site_char.tex

sed -i 's/caption{}/caption{The sites at which species seedlings were sampled, with the number of seedlings successfully sampled per site.}/g' manuscript/include/species_elevcode_tally.tex

# Run latex
latexmk -pdf -quiet -bibtex -cd manuscript/elev.tex
