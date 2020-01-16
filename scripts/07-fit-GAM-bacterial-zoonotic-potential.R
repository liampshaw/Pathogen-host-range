# This code was originally written by Olival et al. (2017)
# and was adapted (lightly) by Liam Shaw 2019 (liam.philip.shaw at gmail dot com)
# for this project.

# See: https://zenodo.org/record/807517 for the original code repository this code was sourced from

# This file is based on: scripts/04-fit-models.R (fitting) and
#                        scripts/09-make-Figure04-viral-traits (plotting) 
# from that repository

# I am grateful to Olival et al. for making their original code available under an MIT License, which also applies here. 
# https://opensource.org/licenses/MIT
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# Make ideal data frames for predicting zoonotic potential
source('scripts/Olival-functions.R')

# Read in viral zoonotic potential dataset 
# N.B. Generated from pathogen_vs_host_db association dataset, with additional
# info for each species. Scripts not included as part of this repository.
bacteria.df <- read.csv('data/zoonotic-potential-dataset-bacteria.csv', header=T, stringsAsFactors = T)

PHB_centers <- names(bacteria.df)[stri_detect_regex(names(bacteria.df), "PHB")]

# Produce dataset for fitting GAMs
# PHBs: numeric
# Records: numeric
# Motility: binary
# Spore: binary
# Vector borne: binary
# Gram: dummy (see below)
# Oxygen: dummy (see below)
# Could also include family, but it removes other variables and each 
# has very small contribution...
abundant <- table(bacteria.df$family)>10
abundant.families <- names(abundant)[which(abundant==TRUE)]
bacteria.df$fA <- sapply(bacteria.df$family, function(x) 
  ifelse(x %in% abundant.families, x, "Other"))
data_set <- bacteria.df[c("PHB.max", "PHB.mean", "PHB.median",
                          "PubMed.records.log", "SRA.records.log", "Nucleotide.records.log", 
                          "motility",
                          "spore",
                          "vector.borne",
                          "gram",
                          "oxygen",
                          "zoonotic", 
                          "fA")]
# Exclude some (noted in text)
data_set <- data_set[which(data_set$motility %in% c("Motile", "Non-motile")),]
data_set <- data_set[which(data_set$spore %in% c("No", "Yes")),]
data_set <- data_set[which(data_set$gram %in% c("Gram-negative", "Gram-positive", "Gram-variable", "Lack cell well")),]
# Drop levels
data_set$motility <- droplevels(data_set$motility)
data_set$spore <- droplevels(data_set$spore)
data_set$vector.borne <- as.factor(data_set$vector.borne)
data_set$gram <- droplevels(data_set$gram)
data_set$oxygen <- droplevels(data_set$oxygen)

# To stop values breaking the fit (as interpreted in formulas - parsed later)
data_set$oxygen <- as.factor(gsub(" ", "", data_set$oxygen))
data_set$gram <- as.factor(gsub("-", "", data_set$gram))

# Model terms 
terms = list(
  PHB     = paste0("s(", PHB_centers, ", bs='tp', k=7)"),
  bias   = c("s(PubMed.records.log, bs='tp', k=7)", "s(SRA.records.log, bs='tp', k=7)", "s(Nucleotide.records.log, bs='tp', k=7"),
  motility = c("s(motility, bs='re')"),
  spore = "s(spore, bs='re')", 
  vector = c("s(vector.borne, bs='re')"), 
  stringsAsFactors = FALSE)


# Add dummy terms for oxygen and Gram
dummys = as.data.frame(cbind(with(data_set, model.matrix(~oxygen))[,-1],
                             with(data_set, model.matrix(~gram))[,-1]))

# Combine with dataset
original_data_set = data_set
data_set = cbind(original_data_set, dummys)
dummy_terms = paste0("s(", names(dummys), ", bs = 're')")
names(dummy_terms) <- names(dummys)



# Add to other model terms
original_terms <- terms
terms <- c(original_terms, dummy_terms)

# Motility and vector borne should become numeric: this otherwise
# causes errors below when you try to fit to new data to get the partial
# effects of the variables
data_set$motility <- as.numeric(data_set$motility)-1
data_set$vector.borne <- as.numeric(data_set$vector.borne)-1


set.seed(0)

gam.fit <- fit_all_gams(data_set,
                        outcome_variable = "zoonotic",
                        terms)
saveRDS(gam.fit, file='intermediates/zoonotic-GAM-fits-bacteria.rds')

# Can also run with family (not run)
# Family
# Maybe family too?
#dummys_family = as.data.frame(cbind(with(data_set, model.matrix(~oxygen))[,-1],
#                                    with(data_set, model.matrix(~gram))[,-1], 
#                                    with(data_set, model.matrix(~fA))[,-1]))
# With bacterial families (not used)
#data_set_family = cbind(original_data_set, dummys_family)
#dummy_terms_family = paste0("s(", names(dummys_family), ", bs = 're')")
#dummy_terms_family = paste0("s(", names(dummys_family), ", bs = 're')")
#names(dummy_terms_family) <- names(dummys_family)
#terms_family <- c(original_terms, dummy_terms_family)
# With family (not run)
#data_set_family$motility <- as.numeric(data_set_family$motility)-1
#data_set_family$vector.borne <- as.numeric(data_set_family$vector.borne)-1
#gam.fit.family <- fit_all_gams(data_set_family,
#                               outcome_variable = "zoonotic",
#                               model_family = binomial,
#                               terms_family)
#saveRDS(gam.fit.family, file='bacterial-GAM-fits-family.rds')

# Now analyse (from 09-make-Figure04-viral-traits.R)
partials_theme = theme(text = element_text(family="Helvetica", size=9),
                       panel.border=element_blank(),
                       panel.background=element_blank(),
                       #   axis.title.y = element_blank(),
                       panel.grid = element_blank(),
                       axis.ticks.x = element_line(size=0.3),
                       axis.ticks.y = element_blank(),
                       axis.text = element_text(color="black"),
                       axis.title.x = element_text(lineheight=1.2),
                       legend.position="none")
bgam = gam.fit$model[[1]]
de_bgam =  get_relative_contribs(bgam) %>%
  mutate(dev_explained = rel_deviance_explained * summary(bgam)$dev.expl) %>%
  mutate(dev_explained = paste0(stri_trim_both(formatC(dev_explained*100, format = "fg", digits=2)), "%"))
de_bgam # Really interesting! PHB explains 17% of deviation, vs. 13% for publication bias


binary_vars = c("vector.borne", "motility", 
                data.frame(de_bgam)$term[grep("gram", data.frame(de_bgam)$term)],
                data.frame(de_bgam)$term[grep("oxygen", data.frame(de_bgam)$term)])
preds <- predict(bgam, type="iterms", se.fit=TRUE)
intercept <- attributes(preds)$constant

preds = preds %>% map(as.data.frame) %>%
  map(~setNames(., stri_replace_first_regex(names(.) ,"s\\(([^\\)]+)\\)", "$1")))
model_data = bgam$model
gterms = attributes(terms(bgam))
ilfun <- bgam$family$linkinv
lfun <- bgam$family$linkfun

binary_terms = which(stri_detect_regex(names(model_data), paste0("(", paste0(binary_vars, collapse="|"), ")")))

smooth_data = model_data[, -c(binary_terms, gterms$response, gterms$offset)]
smooth_ranges = data.frame(map(smooth_data, ~seq(min(.), max(.), length.out = 100)))

binary_data = model_data[, binary_terms, drop=FALSE]
binary_ranges = setNames(as.data.frame(diag(length(binary_terms))), names(model_data)[binary_terms])

offset_name = stri_replace_first_regex(names(model_data)[gterms$offset], "offset\\(([^\\)]+)\\)", "$1")


smooth_ranges = cbind(smooth_ranges,
                      setNames(as.data.frame(lapply(c(names(binary_data), offset_name), function(x) rep(0, nrow(smooth_ranges)))), c(names(binary_data), offset_name)))

binary_ranges = cbind(binary_ranges,
                      setNames(as.data.frame(lapply(c(names(smooth_data), offset_name), function(x) rep(as.factor(0), nrow(binary_ranges)))), c(names(smooth_data), offset_name)))

# Prediction is complicated by the factor levels somehow not matching between original and new data
# "factor levels 0 not in original fit"
# Which has now been fixed
smooth_preds <- predict(bgam, newdata=smooth_ranges, type="iterms", se.fit=TRUE)  %>% map(as.data.frame) %>%
  map(~setNames(., stri_replace_first_regex(names(.) ,"s\\(([^\\)]+)\\)", "$1")))

# Hoever, still broken here: these are factors, need to be converted to numeric
binary_ranges$PubMed.records.log <- 0
binary_ranges$PHB.median <- 0

binary_preds <- predict(bgam, newdata=binary_ranges, type="iterms", se.fit=TRUE) %>% map(as.data.frame) %>%
  map(~setNames(., stri_replace_first_regex(names(.) ,"s\\(([^\\)]+)\\)", "$1")))

partials <- as.data.frame(lapply(1:ncol(preds$fit), function(cl) {
  x = bgam$residuals - rowSums(preds$fit[,-cl]) - intercept
  # (rowSums(preds$fit[,-cl]))
  # x = (lfun(model_data[[gterms$response]]) - model_data[[gterms$offset]])
  #  x = model_data[[gterms$response]]
  x
}))
names(partials) <- names(preds$fit)

smooth_titles = c("PHB median (log)", "PubMed records (log)")
names(smooth_titles) = names(smooth_data)
smooth_plots = map(names(smooth_data), function(smooth_term) {
  pl =  ggplot() +
    #  geom_hline(yintercept = (intercept), size=0.5, col="red") +
    geom_hline(yintercept = 0, size=0.1, col="grey50") +
    geom_point(mapping = aes(x=model_data[[smooth_term]], y = (partials[[smooth_term]])),
               shape=21, fill="black", col="black", alpha=0.25, size=1.25, stroke=0.3) +
    geom_ribbon(mapping = aes(x = smooth_ranges[[smooth_term]],
                              ymin = (smooth_preds$fit[[smooth_term]] - 2 * smooth_preds$se.fit[[smooth_term]]),
                              ymax = (smooth_preds$fit[[smooth_term]] + 2 * smooth_preds$se.fit[[smooth_term]])),
                alpha = 0.75, fill=ifelse(smooth_term=="vGenomeAveLengthLn", "grey", viridis(5)[4])) +
    geom_line(mapping = aes(x = smooth_ranges[[smooth_term]], y = (smooth_preds$fit[[smooth_term]])), size=0.3) +
    #    annotate("label", x = max(model_data[[smooth_term]]), y = -7.5, label = paste0("DE = ", de_bgam$dev_explained[de_bgam$term == smooth_term]), hjust = 1, size=1.5,  label.size=0, fill="#FFFFFF8C") +
    #  geom_rug(mapping = aes(x =model_data[[smooth_term]]), alpha=0.3) +
    xlab(smooth_titles[smooth_term]) +
    scale_y_continuous(limits=c(-12,12), oob=scales::rescale_none) +
    theme_bw() + partials_theme+
    ylab("Strength of effect")
  return(pl)
  
})


# Binary traits
bin_data = binary_preds %>% map(function(x) {
  x = x[, stri_detect_regex(names(x), paste0("(", paste0(binary_vars, collapse="|"), ")")), drop=FALSE]
  n = names(x)
  x = rowSums(x)
  data_frame(response=x, variable=n)
})



bin_data$fit$se = bin_data$se.fit$response
bin_data = bin_data$fit
bin_data$response = bin_data$response
bin_data$labels = bin_data$variable
#bin_data$labels = stri_replace_first_regex(bin_data$labels, "hHuntedIUCN", "Hunted")
bin_data$signif = summary(bgam)$s.table[stri_detect_regex(rownames(summary(bgam)$s.table), paste0("(", paste0(binary_vars, collapse="|"), ")")), "p-value"] < 0.05
bin_data = bin_data %>%
  arrange(desc(signif), response) %>%
  mutate(no = 1:nrow(bin_data))

bin_partials = lapply(binary_terms, function(x) {
  vals = partials[as.logical(model_data[[x]]), x-1]
  variable = names(model_data)[x]
  variable_rep = rep(variable, length(vals))
  no_rep = rep(bin_data$no[bin_data$variable == variable], length(vals))
  data_frame(variable=variable_rep, partial=vals, no=no_rep)
}) %>% bind_rows

bin_partials %<>%
  filter(variable %in% bin_data$variable[bin_data$signif])
bin_data %<>%
  filter(signif) %>%
  left_join(de_bgam, by=c('variable'='term'))

bin_data$labels <- c("Facultatively\nanaerobic", "Microaerophilic", "Vector-borne")
bin_plot = ggplot() +
  geom_hline(yintercept = 0, size=0.1, col="grey50") +
  geom_point(data=bin_partials, mapping=aes(x=no, y=(partial)), position=position_jitter(width=0.35), shape=21, fill="black", col="black", alpha=0.25, size=1.25, stroke=0.3) +
  geom_rect(data = bin_data, mapping=aes(xmin = no - 0.45, xmax  = no + 0.45, ymin=(response-2*se), ymax=(response+2*se)), fill = viridis(5)[4], alpha = 0.75) +
  geom_segment(data = bin_data, mapping=aes(x=no - 0.45, xend = no + 0.45, y=(response), yend=(response)), col="black", size=0.3) +
  scale_x_continuous(breaks = bin_data$no, labels = bin_data$labels) +
  scale_y_continuous(limits=c(-12,12), oob=scales::rescale_none, name="") +
  partials_theme + theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(),
                         axis.text.x = element_text(color="black", lineheight = 1.2,angle=0,
                                                    vjust=0.5, margin=margin(t=0), family="Helvetica", size=4),
                         axis.text.y = element_text(family="Helvetica", size=5.6),
                         panel.border = element_blank(), axis.line.x=element_blank(), axis.line.y=element_blank(),
                         panel.background = element_rect(fill = "transparent", colour = NA),
                         plot.background = element_rect(fill = "transparent", colour = NA))

# COMBINE THE PLOTS
smooth_plots[[2]] <- smooth_plots[[2]] + ylab("")

allplots_bacterial = cowplot::plot_grid(smooth_plots[[1]], smooth_plots[[2]], bin_plot, nrow=1, rel_widths = c(5.3, 4.35), labels = c("(d)", "(e)", "(f)"), label_size=7)
