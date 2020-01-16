# This code was originally written by Olival et al. (2017)
# and was adapted (lightly) by Liam Shaw 2019 (liam.philip.shaw at gmail dot com)
# for this project.

# See: https://zenodo.org/record/807517 for the original code repository this code was sourced from

# This file is based on: scripts/06-make-Figure02-all-gams.R from that repository
# The difference here is we compare effects for variables for viral richness per host and bacterial richness per host (rather than zoonotic proportion as the second row of plots as in Olival et al.)


# I am grateful to Olival et al. for making their original code available under an MIT License, which also applies here. 
# https://opensource.org/licenses/MIT
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

library(ggplot2)
library(tidyr)
library(purrr)
library(stringi)
library(cowplot)
library(viridis)
library(dplyr)
library(svglite)
library(mgcv)
library(magrittr)
set.seed(0)
source('scripts/Olival_fit_gam.R')

SHOW_DEV_EXPL = FALSE

partials_theme = theme(text = element_text(family="Helvetica", size=9),
                       panel.border=element_blank(),
                       panel.background=element_blank(),
                       axis.title.y = element_blank(),
                       panel.grid = element_blank(),
                       axis.ticks.x = element_line(size=0.3),
                       axis.ticks.y = element_blank(),
                       axis.text = element_text(color="black"),
                       axis.title.x = element_text(lineheight = 1.2),
                       legend.position="none"
                       #plot.margin=margin(l=0)
)

blankPlot <- ggplot()+geom_blank(aes(1,1)) +
  cowplot::theme_nothing()

bgam = readRDS('intermediates/all_viruses_models.rds')$model[[1]]

de_bgam =  get_relative_contribs(bgam) %>%
  mutate(dev_explained = rel_deviance_explained * summary(bgam)$dev.expl) %>%
  mutate(dev_explained = paste0(stri_trim_both(formatC(dev_explained*100, format = "fg", digits=2)), "%"))

# All Viruses Plot
binary_vars = c("hOrder")

preds <- predict(bgam, type="iterms", se.fit=TRUE)
intercept <- attributes(preds)$constant

preds = preds %>% map(as.data.frame) %>%
  map(~setNames(., stri_replace_first_regex(names(.) ,"s\\(([^\\)]+)\\)", "$1")))
model_data_vir = bgam$model
gterms = attributes(terms(bgam))
ilfun <- bgam$family$linkinv
lfun <- bgam$family$linkfun

binary_terms = which(stri_detect_regex(names(model_data_vir), paste0("(", paste0(binary_vars, collapse="|"), ")")))

smooth_data_vir = model_data_vir[, -c(binary_terms, gterms$response, gterms$offset)]
smooth_ranges = data.frame(map(smooth_data_vir, ~seq(min(.), max(.), length.out = 100)))

binary_data = model_data_vir[, binary_terms]
binary_ranges = setNames(as.data.frame(diag(length(binary_terms))), names(model_data_vir)[binary_terms])

offset_name = stri_replace_first_regex(names(model_data_vir)[gterms$offset], "offset\\(([^\\)]+)\\)", "$1")

smooth_ranges = cbind(smooth_ranges,
                      setNames(as.data.frame(lapply(c(names(binary_data), offset_name), function(x) rep(0, nrow(smooth_ranges)))), c(names(binary_data), offset_name)))
binary_ranges = cbind(binary_ranges,
                      setNames(as.data.frame(lapply(c(names(smooth_data_vir), offset_name), function(x) rep(0, nrow(binary_ranges)))), c(names(smooth_data_vir), offset_name)))

smooth_preds <- predict(bgam, newdata=smooth_ranges, type="iterms", se.fit=TRUE)  %>% map(as.data.frame) %>%
  map(~setNames(., stri_replace_first_regex(names(.) ,"s\\(([^\\)]+)\\)", "$1")))
binary_preds <- predict(bgam, newdata=binary_ranges, type="iterms", se.fit=TRUE) %>% map(as.data.frame) %>%
  map(~setNames(., stri_replace_first_regex(names(.) ,"s\\(([^\\)]+)\\)", "$1")))

partials <- as.data.frame(lapply(1:ncol(preds$fit), function(cl) {
  x = lfun(model_data_vir[[gterms$response]]) - rowSums(preds$fit[,-cl]) - intercept
  # (rowSums(preds$fit[,-cl]))
  # x = (lfun(model_data_vir[[gterms$response]]) - model_data_vir[[gterms$offset]])
  #  x = model_data_vir[[gterms$response]]
  x
}))
names(partials) <- names(preds$fit)

smooth_titles = list("Disease\ncitations (log)", "PVR body mass", bquote('Range (log' ~ km^{2} ~ ')'), "Mammal\nsympatry")
names(smooth_titles) = names(smooth_data_vir)
smooth_plots_vir = map(names(smooth_data_vir), function(smooth_term_vir) {
  pl =  ggplot() +
    geom_hline(yintercept = 0, size=0.1, col="grey50") +
    geom_point(mapping = aes(x=model_data_vir[[smooth_term_vir]], y = (partials[[smooth_term_vir]])),
               shape=21, fill="red", col="red", alpha=0.25, size=1, stroke=0.1) +
    geom_ribbon(mapping = aes(x = smooth_ranges[[smooth_term_vir]],
                              ymin = (smooth_preds$fit[[smooth_term_vir]] - 2 * smooth_preds$se.fit[[smooth_term_vir]]),
                              ymax = (smooth_preds$fit[[smooth_term_vir]] + 2 * smooth_preds$se.fit[[smooth_term_vir]])),
                alpha = 0.5, fill=viridis(5)[4]) +
    geom_line(mapping = aes(x = smooth_ranges[[smooth_term_vir]], y = (smooth_preds$fit[[smooth_term_vir]])), size=0.3) +
    xlab(smooth_titles[[smooth_term_vir]]) +
    scale_y_continuous(limits=c(-2.2,2.2), oob=scales::rescale_none) +
    theme_bw() + partials_theme
  
  #if (SHOW_DEV_EXPL) pl <- pl + annotate("label", x = max(model_data_vir[[smooth_term_vir]]), y = -2, label = paste0("DE = ", de_bgam$dev_explained[de_bgam$term == smooth_term_vir]), hjust = 1, size=1.5,  label.size=0, fill="#FFFFFF8C")
  #  geom_rug(mapping = aes(x =model_data_vir[[smooth_term]]), alpha=0.3) +
  
  return(pl)
  
})

smooth_plots_vir[[1]] = smooth_plots_vir[[1]] + ylab("Strength of effect on\nviral richness") +
  theme(axis.title.y=element_text(angle=90, lineheight=1.2, margin=margin(r=3)))
smooth_plots_vir[[2]] = smooth_plots_vir[[2]] + theme(axis.title.y = element_blank())
smooth_plots_vir[[3]] = smooth_plots_vir[[3]] + theme(axis.title.y = element_blank())
smooth_plots_vir[[4]] = smooth_plots_vir[[4]] + theme(axis.title.y = element_blank())
bin_vir_data = binary_preds %>% map(function(x) {
  x = x[, stri_detect_regex(names(x), paste0("(", paste0(binary_vars, collapse="|"), ")"))]
  n = names(x)
  x = rowSums(x)
  data_frame(response=x, variable=n)
})



bin_vir_data$fit$se = bin_vir_data$se.fit$response
bin_vir_data = bin_vir_data$fit
bin_vir_data$response = bin_vir_data$response
bin_vir_data$labels = stri_replace_first_regex(bin_vir_data$variable, "hOrder", "")
#bin_vir_data$labels = stri_replace_first_regex(bin_vir_data$labels, "hHuntedIUCN", "Hunted")
bin_vir_data$signif = summary(bgam)$s.table[stri_detect_regex(rownames(summary(bgam)$s.table), paste0("(", paste0(binary_vars, collapse="|"), ")")), "p-value"] < 0.05
bin_vir_data = bin_vir_data %>%
  arrange(desc(signif), response) %>%
  mutate(no = 1:nrow(bin_vir_data))

bin_vir_partials = lapply(binary_terms, function(x) {
  vals = partials[as.logical(model_data_vir[[x]]), x-1]
  variable = names(model_data_vir)[x]
  variable_rep = rep(variable, length(vals))
  no_rep = rep(bin_vir_data$no[bin_vir_data$variable == variable], length(vals))
  data_frame(variable=variable_rep, partial=vals, no=no_rep)
}) %>% bind_rows

bin_vir_data = bin_vir_partials %>%
  group_by(variable) %>%
  summarize(minval = min(partial)) %>%
  inner_join(bin_vir_data, by="variable") %>%
  mutate(minval = pmin(minval, response - 2*se)) %>%
  left_join(de_bgam, by=c('variable' = 'term'))

bin_plot_vir = ggplot() +
  geom_hline(yintercept = 0, size=0.1, col="grey50") +
  geom_point(data=bin_vir_partials, mapping=aes(x=no, y=(partial)), position=position_jitter(width=0.25),
             shape=21, fill="red", col="red", alpha=0.25, size=1, stroke=0.1) +
  geom_rect(data = bin_vir_data, mapping=aes(xmin = no - 0.35, xmax  = no + 0.35, ymin=(response-2*se), ymax=(response+2*se), fill=signif), alpha = 0.5) +
  geom_segment(data = bin_vir_data, mapping=aes(x=no - 0.35, xend = no + 0.35, y=(response), yend=(response)), col="black", size=0.3) +
  scale_fill_manual(breaks = c(TRUE, FALSE), values=c(viridis(5)[4], "grey")) +
  scale_x_continuous(breaks = bin_vir_data$no, labels = stri_trans_totitle(bin_vir_data$labels)) +
  scale_y_continuous(limits=c(-2,2.2), breaks=seq(-2,2, by=1), oob=scales::rescale_none, name="") +
  theme_bw() + partials_theme +
  theme(axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.text.x=element_text(colour="black", angle=90, hjust=1), legend.position="none",
        axis.title.y=element_blank())

#if (SHOW_DEV_EXPL) bin_plot_vir = bin_plot_vir + geom_label(data = bin_vir_data, mapping=aes(x=no, y = response + 2*se + 0.4, label=dev_explained), color="black", family="Helvetica", size=1.5, label.size=0, fill="#FFFFFF8C") +


vir_plots <- plot_grid(plot_grid(plotlist = smooth_plots_vir, nrow=1, align="h", rel_widths = c(1.22,1,1,1),
                                 labels=c("(a)", "(b)", "(c)", "(d)"), label_size=7, hjust=0),
                       bin_plot_vir, nrow=1, rel_widths = c(4.22,1.3), labels=c("", "(e)"), label_size=7, hjust=0)


#---
# For bacteria


bgam = readRDS('intermediates/all_bacteria_models.rds')$model[[1]]

de_bgam =  get_relative_contribs(bgam) %>%
  mutate(dev_explained = rel_deviance_explained * summary(bgam)$dev.expl) %>%
  mutate(dev_explained = paste0(stri_trim_both(formatC(dev_explained*100, format = "fg", digits=2)), "%"))

# All bacuses Plot
binary_vars = c("hOrder")

preds <- predict(bgam, type="iterms", se.fit=TRUE)
intercept <- attributes(preds)$constant

preds = preds %>% map(as.data.frame) %>%
  map(~setNames(., stri_replace_first_regex(names(.) ,"s\\(([^\\)]+)\\)", "$1")))
model_data_bac = bgam$model
gterms = attributes(terms(bgam))
ilfun <- bgam$family$linkinv
lfun <- bgam$family$linkfun

binary_terms = which(stri_detect_regex(names(model_data_bac), paste0("(", paste0(binary_vars, collapse="|"), ")")))

smooth_data_bac = model_data_bac[, -c(binary_terms, gterms$response, gterms$offset)]
smooth_ranges = data.frame(map(smooth_data_bac, ~seq(min(.), max(.), length.out = 100)))

binary_data = model_data_bac[, binary_terms]
binary_ranges = setNames(as.data.frame(diag(length(binary_terms))), names(model_data_bac)[binary_terms])

offset_name = stri_replace_first_regex(names(model_data_bac)[gterms$offset], "offset\\(([^\\)]+)\\)", "$1")

smooth_ranges = cbind(smooth_ranges,
                      setNames(as.data.frame(lapply(c(names(binary_data), offset_name), function(x) rep(0, nrow(smooth_ranges)))), c(names(binary_data), offset_name)))
binary_ranges = cbind(binary_ranges,
                      setNames(as.data.frame(lapply(c(names(smooth_data_bac), offset_name), function(x) rep(0, nrow(binary_ranges)))), c(names(smooth_data_bac), offset_name)))

smooth_preds <- predict(bgam, newdata=smooth_ranges, type="iterms", se.fit=TRUE)  %>% map(as.data.frame) %>%
  map(~setNames(., stri_replace_first_regex(names(.) ,"s\\(([^\\)]+)\\)", "$1")))
binary_preds <- predict(bgam, newdata=binary_ranges, type="iterms", se.fit=TRUE) %>% map(as.data.frame) %>%
  map(~setNames(., stri_replace_first_regex(names(.) ,"s\\(([^\\)]+)\\)", "$1")))

partials <- as.data.frame(lapply(1:ncol(preds$fit), function(cl) {
  x = lfun(model_data_bac[[gterms$response]]) - rowSums(preds$fit[,-cl]) - intercept
  # (rowSums(preds$fit[,-cl]))
  # x = (lfun(model_data_bac[[gterms$response]]) - model_data_bac[[gterms$offset]])
  #  x = model_data_bac[[gterms$response]]
  x
}))
names(partials) <- names(preds$fit)

smooth_titles = list("Disease\ncitations (log)", bquote('Range (log' ~ km^{2} ~ ')'), "Mammal\nsympatry")
names(smooth_titles) = names(smooth_data_bac)
smooth_plots_bac = map(names(smooth_data_bac), function(smooth_term_bac) {
  pl =  ggplot() +
    geom_hline(yintercept = 0, size=0.1, col="grey50") +
    geom_point(mapping = aes(x=model_data_bac[[smooth_term_bac]], y = (partials[[smooth_term_bac]])),
               shape=21, fill="black", col="black", alpha=0.25, size=1, stroke=0.1) +
    geom_ribbon(mapping = aes(x = smooth_ranges[[smooth_term_bac]],
                              ymin = (smooth_preds$fit[[smooth_term_bac]] - 2 * smooth_preds$se.fit[[smooth_term_bac]]),
                              ymax = (smooth_preds$fit[[smooth_term_bac]] + 2 * smooth_preds$se.fit[[smooth_term_bac]])),
                alpha = 0.5, fill=viridis(5)[4]) +
    geom_line(mapping = aes(x = smooth_ranges[[smooth_term_bac]], y = (smooth_preds$fit[[smooth_term_bac]])), size=0.3) +
    xlab(smooth_titles[[smooth_term_bac]]) +
    scale_y_continuous(limits=c(-3,2.2), oob=scales::rescale_none) +
    theme_bw() + partials_theme
  
  #if (SHOW_DEV_EXPL) pl <- pl + annotate("label", x = max(model_data_bac[[smooth_term_bac]]), y = -2, label = paste0("DE = ", de_bgam$dev_explained[de_bgam$term == smooth_term_bac]), hjust = 1, size=1.5,  label.size=0, fill="#FFFFFF8C")
  #  geom_rug(mapping = aes(x =model_data_bac[[smooth_term]]), alpha=0.3) +
  
  return(pl)
  
})

smooth_plots_bac[[1]] = smooth_plots_bac[[1]] + ylab("Strength of effect on\nbacterial richness") +
  theme(axis.title.y=element_text(angle=90, lineheight=1.2, margin=margin(r=3)))
smooth_plots_bac[[2]] = smooth_plots_bac[[2]] + theme(axis.title.y = element_blank())
smooth_plots_bac[[3]] = smooth_plots_bac[[3]] + theme(axis.title.y = element_blank())
# Add blank 4th plot
smooth_plots_bac[[4]] = blankPlot
bin_bac_data = binary_preds %>% map(function(x) {
  x = x[, stri_detect_regex(names(x), paste0("(", paste0(binary_vars, collapse="|"), ")"))]
  n = names(x)
  x = rowSums(x)
  data_frame(response=x, variable=n)
})



bin_bac_data$fit$se = bin_bac_data$se.fit$response
bin_bac_data = bin_bac_data$fit
bin_bac_data$response = bin_bac_data$response
bin_bac_data$labels = stri_replace_first_regex(bin_bac_data$variable, "hOrder", "")
#bin_bac_data$labels = stri_replace_first_regex(bin_bac_data$labels, "hHuntedIUCN", "Hunted")
bin_bac_data$signif = summary(bgam)$s.table[stri_detect_regex(rownames(summary(bgam)$s.table), paste0("(", paste0(binary_vars, collapse="|"), ")")), "p-value"] < 0.05
bin_bac_data = bin_bac_data %>%
  arrange(desc(signif), response) %>%
  mutate(no = 1:nrow(bin_bac_data))

bin_bac_partials = lapply(binary_terms, function(x) {
  vals = partials[as.logical(model_data_bac[[x]]), x-1]
  variable = names(model_data_bac)[x]
  variable_rep = rep(variable, length(vals))
  no_rep = rep(bin_bac_data$no[bin_bac_data$variable == variable], length(vals))
  data_frame(variable=variable_rep, partial=vals, no=no_rep)
}) %>% bind_rows

bin_bac_data = bin_bac_partials %>%
  group_by(variable) %>%
  summarize(minval = min(partial)) %>%
  inner_join(bin_bac_data, by="variable") %>%
  mutate(minval = pmin(minval, response - 2*se)) %>%
  left_join(de_bgam, by=c('variable' = 'term'))

bin_plot_bac = ggplot() +
  geom_hline(yintercept = 0, size=0.1, col="grey50") +
  geom_point(data=bin_bac_partials, mapping=aes(x=no, y=(partial)), position=position_jitter(width=0.25),
             shape=21, fill="black", col="black", alpha=0.25, size=1, stroke=0.1) +
  geom_rect(data = bin_bac_data, mapping=aes(xmin = no - 0.35, xmax  = no + 0.35, ymin=(response-2*se), ymax=(response+2*se), fill=signif), alpha = 0.5) +
  geom_segment(data = bin_bac_data, mapping=aes(x=no - 0.35, xend = no + 0.35, y=(response), yend=(response)), col="black", size=0.3) +
  scale_fill_manual(breaks = c(TRUE, FALSE), values=c(viridis(5)[4], "grey")) +
  scale_x_continuous(breaks = bin_bac_data$no, labels = stri_trans_totitle(bin_bac_data$labels)) +
  scale_y_continuous(limits=c(-3,2.2), breaks=seq(-2,2, by=1), oob=scales::rescale_none, name="") +
  theme_bw() + partials_theme +
  theme(axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.text.x=element_text(colour="black", angle=90, hjust=1), legend.position="none",
        axis.title.y=element_blank())

bac_plots <- plot_grid(plot_grid(plotlist = smooth_plots_bac, nrow=1, align="h", rel_widths = c(1.22,1,1,1),
                                 labels=c("(f)", "(g)", "(h)"), label_size=7, hjust=0),
                       bin_plot_bac, nrow=1, rel_widths = c(4.22,1.3), labels=c("", "(i)"), label_size=7, hjust=0)

# COMBINE THE PLOTS
allplots = cowplot::plot_grid(vir_plots, bac_plots, nrow=2, rel_widths = c(5.3, 4.35))

# Save pdf
pdf(file='figures/Figure-4-GAMs-richness-by-species.pdf', width =7, height=5, pointsize=7)
allplots
dev.off()

# Save png
png(file="figures/Figure-4-GAMs-richness-by-species.png", width = 600, height=500, pointsize = 7, res=300)
allplots
dev.off()

