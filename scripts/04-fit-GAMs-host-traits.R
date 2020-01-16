# This code was originally written by Olival et al. (2017)
# and was adapted (lightly) by Liam Shaw 2019 (liam.philip.shaw at gmail dot com)
# for this project.

# See: https://zenodo.org/record/807517 for the original code repository this code was sourced from

# This file is based on: scripts/04-fit-models.R from that repository

# I am grateful to Olival et al. for making their original code available under an MIT License, which also applies here. 
# https://opensource.org/licenses/MIT
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

library(mgcv)
library(dplyr)
library(stringi)
library(parallel)
library(purrr)
library(ggplot2)
library(viridis)
library(knitr)
library(svglite)
source('scripts/Olival-functions.R')

olival.db <- readRDS("data/Olival-et-al-database.rds")
hosts <- olival.db$hosts
olival.hosts <- gsub("_", "", as.character(hosts$hHostNameFinal))

# This study
this.study.hosts <- unique(pathogen_vs_host_db$HostSpeciesPHB)
# Overlap
hosts.to.use <- hosts$hHostNameFinal[which(olival.hosts %in% this.study.hosts)]
hosts <- hosts[which(hosts$hHostNameFinal %in% hosts.to.use),]
hosts$hHostNamePHB <- gsub("_", "", as.character(hosts$hHostNameFinal))

# Now can add in our info for total bacterial / viral richness
this.study.hosts.to.use <- gsub("_", "", as.character(hosts$hHostNameFinal))
# Find total bacterial and viral richness for each species
boxplot.df.bacteria.overlap <- boxplot.df.bacteria[which(boxplot.df.bacteria$HostSpeciesPHB %in% this.study.hosts.to.use),]
hosts.bacteria <- hosts[which(hosts$hHostNamePHB %in% boxplot.df.bacteria.overlap$HostSpeciesPHB),]
hosts.bacteria$TotBacteriaPerHost <- boxplot.df.bacteria.overlap$nPathogens

####---- All Bacteria GAM - All Associations ----

data_set = hosts.bacteria %>%
  filter(hMarOTerr == "Terrestrial",
         hWildDomFAO == "wild",
         !is.na(PdHoSa.cbCst_order))


#  Create dummy variables for orders to use as random effects
# But note that these become non-conformable
# Could just leave them out?
data_set$hOrder <- droplevels(data_set$hOrder)
dummys = as.data.frame(with(data_set, model.matrix(~hOrder))[,-1])
# The problem (from pairs plot) is hOrderSCANDENTIA and hOrderPERAMELEMORPHIA which both have zero
#dummys = dummys[, colSums(dummys) > 1]
data_set = cbind(data_set, dummys)
dummy_terms = paste0("s(", names(dummys), ", bs = 're')")
names(dummy_terms) <- names(dummys)

## Create data.frame of all possible models
terms = list(
  s1 = "s(LnAreaHost, bs='cs', k = 7)",
  s2 = "s(hMassGramsPVR, bs='cs', k = 7)",
  s3 = c(
    "s(S100, bs='cs', k = 7)",
    "s(S80, bs='cs', k = 7)",
    "s(S50, bs='cs', k = 7)",
    "s(S40, bs='cs', k = 7)",
    "s(S20, bs='cs', k = 7)",
    "s(S, bs='cs', k = 7)"
  ),
  bias = c("s(hDiseaseZACitesLn, bs='cs', k=7)"),
  stringsAsFactors = FALSE)

# The dummy terms might be causing problems - so leave it out
terms = c(dummy_terms, terms)

all_bacteria = fit_all_gams_poisson(data_set,
                            outcome_variable = "TotBacteriaPerHost",
                            terms)

saveRDS(all_bacteria, "intermediates/all_bacteria_models.rds")

# Viruses
boxplot.df.virus.overlap <- boxplot.df.virus[which(boxplot.df.virus$HostSpeciesPHB %in% this.study.hosts.to.use),]
hosts.virus <- hosts[which(hosts$hHostNamePHB %in% boxplot.df.virus.overlap$HostSpeciesPHB),]
hosts.virus$TotVirusPerHostThisStudy <- boxplot.df.virus.overlap$nPathogens

####---- All Viruses GAM - All Associations ----

data_set = hosts.virus %>%
  filter(hMarOTerr == "Terrestrial",
         hWildDomFAO == "wild",
         !is.na(PdHoSa.cbCst_order))

#  Create dummy variables for orders to use as random effects
dummys = as.data.frame(with(data_set, model.matrix(~hOrder))[,-1])
#dummys = dummys[, colSums(dummys) > 1]
data_set = cbind(data_set, dummys)
dummy_terms = paste0("s(", names(dummys), ", bs = 're')")
names(dummy_terms) <- names(dummys)

## Create data.frame of all possible models
terms = list(
  s1 = "s(LnAreaHost, bs='cs', k = 7)",
  s2 = "s(hMassGramsPVR, bs='cs', k = 7)",
  s3 = c(
    "s(S100, bs='cs', k = 7)",
    "s(S80, bs='cs', k = 7)",
    "s(S50, bs='cs', k = 7)",
    "s(S40, bs='cs', k = 7)",
    "s(S20, bs='cs', k = 7)",
    "s(S, bs='cs', k = 7)"
  ),
  bias = c("s(hDiseaseZACitesLn, bs='cs', k=7)"),
  stringsAsFactors = FALSE)

terms = c(dummy_terms, terms)

all_viruses = fit_all_gams_poisson(data_set,
                           outcome_variable = "TotVirusPerHostThisStudy",
                           terms)

saveRDS(all_viruses, "intermediates/all_viruses_models.rds")
