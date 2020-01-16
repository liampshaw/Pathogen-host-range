# This code was originally written by Olival et al. (2017)
# and was adapted (lightly) by Liam Shaw 2019 (liam.philip.shaw at gmail dot com)
# for this project.

# See: https://zenodo.org/record/807517 for the original code repository this code was sourced from

# This file is based on: scripts/11-make-ExtendedTable01-models.R
# from that repository

# I am grateful to Olival et al. for making their original code available under an MIT License, which also applies here. 
# https://opensource.org/licenses/MIT
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

library(mgcv)
library(stringi)
library(dplyr)
library(purrr)
library(officer)
source('scripts/Olival-functions.R')

top_models <- lapply(c("intermediates/all_bacteria_models.rds", "intermediates//all_viruses_models.rds", "intermediates/zoonotic-GAM-fits-bacteria.rds", "intermediates/zoonotic-GAM-fits-virus.rds"), function(mods) {
  readRDS(mods)$model[[1]]
})

model_names = c("Bacterial richness per species",
                "Viral richness per species",
                "Zoonotic potential (bacteria)",
                "Zoonotic potential (viruses)")

model_tables = map2(top_models, model_names, function(modd, model_name) {
  summ = summary(modd)
  summ$p.table
  summ$s.table
  rel_dev = get_relative_contribs(modd)
  
  bind_rows(data_frame(Term = stri_extract_first_regex(rownames(summ$p.table), "(?<=\\()[^\\)]+(?=\\))"),
                       Value = round(summ$p.table[,1], 2),
                       `Z statistic` = round(summ$p.table[,3], 2),
                       `Chi-sq statistic` = NA,
                       `P-value` = ifelse(summ$p.table[,4] > 0.001, as.character(round(summ$p.table[,4], digits=3)), "<0.001"),
                       `Effective Degrees of Freedom` = NA,
                       `Percent Dev. Explained` = as.character(NA),
                       model = model_name),
            data_frame(Term = stri_extract_first_regex(rownames(summ$s.table), "(?<=s\\()[^\\)]+(?=\\))"),
                       Value = NA,
                       `Z statistic` = NA,
                       `Chi-sq statistic` = round(summ$s.table[,3], 2),
                       `P-value` = ifelse(summ$s.table[,4] > 0.001, as.character(round(summ$s.table[,4], digits=3)), "<0.001"),
                       `Effective Degrees of Freedom` = round(summ$s.table[,1], 2),
                       `Percent Dev. Explained`= paste0(round(100*rel_dev$rel_deviance_explained, 2), "%"),
                       model=model_name))
  
})

model_rows = map_int(model_tables, nrow)
model_tables2 = model_tables %>%
  map(~ rbind(.[1,], .)) %>%
  map(function(x) {
    x$model = c(x$model[1], rep(NA, nrow(x) -1))
    return(x)
  }) %>%
  bind_rows %>%
  mutate_each(funs(as.character), -Term, -model) #%>%
  #arrange(model, Term !="Intercept") %>%
  dplyr::select(8, 1:7)

names(model_tables2)[1] <- ""

print(kable(model_tables2))
write.csv(model_tables2, file='figures/Table-3-summary-best-fit-GAMs.csv')


