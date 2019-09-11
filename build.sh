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

sed -i 's/ccccccccc/lrcS[table-format=3.2]S[table-format=3.2]S[table-format=3.2]S[table-format=3.2]S[table-format=3.2]S[table-format=3.2]/g' manuscript/include/site_char.tex

sed -i 's/ccccc/llr@{\\hspace{0.2\\tabcolsep}}lr@{\\hspace{0.2\\tabcolsep}}lr@{\\hspace{0.2\\tabcolsep}}l/g' manuscript/include/species_elevcode_tally.tex

sed -i 's/cccccc/lcS[table-format=3.2]S[table-format=3.2]S[table-format=3.2]S[table-format=3.2]/g' manuscript/include/best_mod_multi_output.tex

## Put column headers in curly braces
sed -i '10s/.*/{Site} \& {Trees ha\\textsuperscript{-1}} \& {$C_R$ (m)} \\\\/' manuscript/include/comp_radius.tex

sed -i '10s/.*/{Fixed eff.} \& {AIC} \& {$\Delta{}AIC_r$} \& {$W_i$} \& {$R^2_c$} \& {$R^2_m$} \& {Slope} \& {SE} \\\\/' manuscript/include/*_model_comparison.tex

sed -i '10s/.*/{Site} \& {Elev. (m (m))} \& {Precip. (mm y\\textsuperscript{-1})} \& { Mean temp. (\\textdegree{}C)} \& {Soil C (\\%)} \& {Soil N (\\%)} \& {Soil pH} \& {Trees ha\\textsuperscript{-1}} \\\\/' manuscript/include/site_char.tex

sed -i '10s/.*/{Species code} \& {Species} \& \\multicolumn{2}{c}{Bottom} \& \\multicolumn{2}{c}{Middle} \& \\multicolumn{2}{c}{Top} \\\\/' manuscript/include/species_elevcode_tally.tex

sed -i '10s/.*/{Response} \& {Fixed effects} \& {$\\Delta{}AIC_r$} \& {$R^2_c$} \& {$R^2_m$} \\\\/' manuscript/include/best_mod_multi_output.tex

## Format column positions by equals sig (m)n
sed -i 's/=/ \& =/g' manuscript/include/species_elevcode_tally.tex
sed -i 's/NA /NA \& /g' manuscript/include/species_elevcode_tally.tex

## Put NA in empty cells in site characteristics
sed -i '12,21s/  / NA /g' manuscript/include/site_char.tex

## Adjust column sep 
sed -i 's/extracolsep{5pt}/extracolsep{-2pt}/' manuscript/include/site_char.tex

## Format fixed effects in *_model_comparison.tex
sed -i 's/elev\\_scale.comp\\_adult\\_metric\\_log\\_scale.lai\\_scale \&/Elev. + ISI + LAI \&/' manuscript/include/*_model_comparison.tex

sed -i 's/elev\\_scale.comp\\_adult\\_metric\\_log\\_scale \&/Elev. + ISI \&/' manuscript/include/*_model_comparison.tex

sed -i 's/elev\\_scale.lai\\_scale \&/Elev. + LAI \&/' manuscript/include/*_model_comparison.tex

sed -i 's/comp\\_adult\\_metric\\_log\\_scale.lai\\_scale \&/ISI + LAI \&/' manuscript/include/*_model_comparison.tex

sed -i 's/elev\\_scale \&/Elev. \&/' manuscript/include/*_model_comparison.tex

sed -i 's/lai\\_scale \&/LAI \&/' manuscript/include/*_model_comparison.tex

sed -i 's/comp\\_adult\\_metric\\_log\\_scale \&/ISI \&/' manuscript/include/*_model_comparison.tex

## Format fixed effects in best_model_multi_output.tex
sed -i 's/d\\_fvfm/\$F_v\/F_m\$/' manuscript/include/best_mod_multi_output.tex

sed -i 's/leaf\\_chl/Chl\\textsubscript{$\\alpha$}/' manuscript/include/best_mod_multi_output.tex

sed -i 's/leaf\\_height\\_ratio/Leaf:height ratio/'  manuscript/include/best_mod_multi_output.tex

sed -i 's/leaf\\_area\\_cm2\\_log/log(Leaf area)/'  manuscript/include/best_mod_multi_output.tex

sed -i 's/stem\\_vol\\_cm3/Stem vol./'  manuscript/include/best_mod_multi_output.tex

sed -i 's/leaf\\_thick\\_mean\\_mm/Leaf thickness/'  manuscript/include/best_mod_multi_output.tex

sed -i 's/elev\\_scale.comp\\_adult\\_metric\\_log\\_scale.lai\\_scale/Elev. + ISI + LAI/'  manuscript/include/best_mod_multi_output.tex

sed -i 's/lai\\_scale/LAI/'  manuscript/include/best_mod_multi_output.tex

## Add decimal places to some integers
#site_char.tex

## Make species names italics
sed -i -r '12,15s/\w+\s\w+/\\textit{&}/g' manuscript/include/species_elevcode_tally.tex
sed -i -r '16s/Myrcia/\\textit{&}/g' manuscript/include/species_elevcode_tally.tex
sed -i -r '17,18s/\w+\s\w+/\\textit{&}/g' manuscript/include/species_elevcode_tally.tex

## Add captions
sed -i 's/caption{}/caption{Competition radius ($C_R$) used for adult competition measurements for each site based on the number of trees per hectare.}/g' manuscript/include/comp_radius.tex

sed -i 's/caption{}/caption{Site environmental characteristics for each 1 ha plot sampled. NA indicates that no data was available. Adapted from \\citet{Whitaker2014}.}/g' manuscript/include/site_char.tex

sed -i 's/caption{}/caption{The sites at which tree seedlings were sampled for each species, with the number of seedlings successfully sampled per site.}/g' manuscript/include/species_elevcode_tally.tex


# Run latex
latexmk -pdf -quiet -bibtex -cd manuscript/elev.tex
