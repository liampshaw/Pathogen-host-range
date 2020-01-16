# This code was originally written by Olival et al. (2017)
# and was adapted (lightly) by Liam Shaw 2019 (liam.philip.shaw at gmail dot com)
# for this project.

# See: https://zenodo.org/record/807517 for the original code repository this code was sourced from
# I am grateful to Olival et al. for making their code available under an MIT License.
# https://opensource.org/licenses/MIT
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

fit_all_gams <- function(data_set, outcome_variable, terms) {

  fit_gam = function(frm) {
    try(gam(formula=as.formula(frm), family="binomial", data_set, select=TRUE), silent=TRUE)
  }

  terms_grid = do.call(expand.grid, terms)

  #Create model forumulas from the grid
  formulas = apply(as.matrix(terms_grid), 1, function(row) paste(row, collapse = " + ")) %>%
    stri_replace_all_regex("\\s[\\+\\s]+\\s", " + ") %>%
    {paste(outcome_variable, "~", .)} %>%
    rearrange_formula %>%
    unique

  models = data_frame(formula = formulas)

  n_cores = detectCores()
  n_cores_use = round(nrow(models) / ceiling(nrow(models) / (n_cores - 1)))
  options(mc.cores = n_cores_use)
  message("Using ", n_cores_use, " cores to fit ", nrow(models), " models")

  models_vec = mclapply(models$formula, fit_gam)

  models = models %>%
    mutate(model = models_vec)
  #print(models)

  # Calculate models
  models = models %>%
    filter(map_lgl(model, ~ !("try-error" %in% class(.) | is.null(.)))) %>%
    mutate(aic = map_dbl(model, MuMIn::AICc),
           daic = aic - min(aic),
           weight = exp(-daic/2)) %>%
    arrange(aic)
  #print("Calculated models")
  #print(models)

  # Remove unused terms from models and reduce to unique ones
  models_reduced = models %>%
    dplyr::select(model) %>%
    mutate(formula = map_chr(model, ~rearrange_formula(rm_low_edf(.)))) %>%
    distinct(formula, .keep_all = TRUE)
  #print("Remove unused terms")
  #print(models)

  n_cores_use = round(nrow(models_reduced) / ceiling(nrow(models_reduced) / (n_cores - 1)))
  options(mc.cores = n_cores_use)
  message("Using ", n_cores_use, " cores to fit ", nrow(models_reduced), " reduced models")

  # Reduce the remaining models
  models_reduced = models_reduced %>%
    mutate(model = mclapply(model, reduce_model))
  #print("Reduced remaining models")
  #print(models_reduced)

  models_reduced = models_reduced %>%
    filter(map_lgl(model, ~ !("try-error" %in% class(.) | is.null(.)))) %>%
    mutate(aic = map_dbl(model, MuMIn::AICc)) %>%
    mutate(formula = map_chr(model, ~rearrange_formula(rm_low_edf(.)))) %>%
    distinct(formula, .keep_all=TRUE) %>%
    arrange(aic) %>%
    mutate(daic = aic - min(aic),
           weight = exp(-daic/2),
           terms = shortform(map(model, ~ rearrange_formula(.$formula))),
           relweight = ifelse(daic > 2, 0, weight/sum(weight[daic < 2])),
           relweight_all = weight/sum(weight),
           cumweight = cumsum(relweight_all))
  #print("Final reduction")
  #print(models_reduced)

  return(models_reduced)

}

fit_all_gams_poisson <- function(data_set, outcome_variable, terms) {
  
  fit_gam = function(frm) {
    try(gam(formula=as.formula(frm), family="poisson", data_set, select=TRUE), silent=TRUE)
  }
  
  terms_grid = do.call(expand.grid, terms)
  
  #Create model forumulas from the grid
  formulas = apply(as.matrix(terms_grid), 1, function(row) paste(row, collapse = " + ")) %>%
    stri_replace_all_regex("\\s[\\+\\s]+\\s", " + ") %>%
    {paste(outcome_variable, "~", .)} %>%
    rearrange_formula %>%
    unique
  
  models = data_frame(formula = formulas)
  
  n_cores = detectCores()
  n_cores_use = round(nrow(models) / ceiling(nrow(models) / (n_cores - 1)))
  options(mc.cores = n_cores_use)
  message("Using ", n_cores_use, " cores to fit ", nrow(models), " models")
  
  models_vec = mclapply(models$formula, fit_gam)
  
  models = models %>%
    mutate(model = models_vec)
  #print(models)
  
  # Calculate models
  models = models %>%
    filter(map_lgl(model, ~ !("try-error" %in% class(.) | is.null(.)))) %>%
    mutate(aic = map_dbl(model, MuMIn::AICc),
           daic = aic - min(aic),
           weight = exp(-daic/2)) %>%
    arrange(aic)
  #print("Calculated models")
  #print(models)
  
  # Remove unused terms from models and reduce to unique ones
  models_reduced = models %>%
    dplyr::select(model) %>%
    mutate(formula = map_chr(model, ~rearrange_formula(rm_low_edf(.)))) %>%
    distinct(formula, .keep_all = TRUE)
  #print("Remove unused terms")
  #print(models)
  
  n_cores_use = round(nrow(models_reduced) / ceiling(nrow(models_reduced) / (n_cores - 1)))
  options(mc.cores = n_cores_use)
  message("Using ", n_cores_use, " cores to fit ", nrow(models_reduced), " reduced models")
  
  # Reduce the remaining models
  models_reduced = models_reduced %>%
    mutate(model = mclapply(model, reduce_model))
  #print("Reduced remaining models")
  #print(models_reduced)
  
  models_reduced = models_reduced %>%
    filter(map_lgl(model, ~ !("try-error" %in% class(.) | is.null(.)))) %>%
    mutate(aic = map_dbl(model, MuMIn::AICc)) %>%
    mutate(formula = map_chr(model, ~rearrange_formula(rm_low_edf(.)))) %>%
    distinct(formula, .keep_all=TRUE) %>%
    arrange(aic) %>%
    mutate(daic = aic - min(aic),
           weight = exp(-daic/2),
           terms = shortform(map(model, ~ rearrange_formula(.$formula))),
           relweight = ifelse(daic > 2, 0, weight/sum(weight[daic < 2])),
           relweight_all = weight/sum(weight),
           cumweight = cumsum(relweight_all))
  #print("Final reduction")
  #print(models_reduced)
  
  return(models_reduced)
  
}


# Returns a model formula from a GAM with the low_edf terms removed
rm_low_edf <- function(mod, edf_cutoff = 0.001) {
  fr = as.character(formula(mod))
  lhs = fr[2]
  rhs = fr[3]
  edfs = pen.edf(mod)
  low_edfs = edfs[edfs < edf_cutoff]
  vars_to_remove = stri_replace_all_fixed(unique(stri_extract_first_regex(names(low_edfs), "(?<=s\\()[^\\)]+(?=\\))")), ",",", ")
  vars_regex = paste0("(", paste(vars_to_remove, collapse="|"), ")")
  new_rhs = stri_replace_all_regex(rhs, paste0("\\s*s\\(", vars_regex, "\\,[^\\)]+\\)\\s*\\+?"), "")
  new_rhs= stri_replace_all_fixed(new_rhs, "+, k = 7) ", "")
  new_rhs= stri_replace_all_fixed(new_rhs, "+ +s", "+ s")
  new_formula = paste(lhs, "~", new_rhs)
  new_formula = stri_replace_all_regex(new_formula, "[\\s\\n]+", " ")
  new_formula = stri_replace_all_regex(new_formula, "[+\\s]*$", "")
  return(new_formula)
}

rm_terms <- function(mod, terms) {
  fr = as.character(formula(mod))
  lhs = fr[2]
  rhs = fr[3]
  vars_regex = paste0("(", paste(terms, collapse="|"), ")")
  new_rhs = stri_replace_all_regex(rhs, paste0("\\s*s\\(", vars_regex, "\\,[^\\)]+\\)\\s*\\+?"), "")
  new_rhs= stri_replace_all_fixed(new_rhs, "+, k = 7) ", "")
  new_rhs= stri_replace_all_fixed(new_rhs, "+ +s", "+ s")
  new_formula = paste(lhs, "~", new_rhs)
  new_formula = stri_replace_all_regex(new_formula, "[\\s\\n]+", " ")
  new_formula = stri_replace_all_regex(new_formula, "[+\\s]*$", "")
  return(new_formula)
}

#' Alphabetizes the right-hand side of a formula so as to compare formulas across models
rearrange_formula = function(formula) {
  if(class(formula) == "formula") {
    formula = as.character(formula)
    formula = paste(formula[2], "~", formula[3], collapse=" ")
  }
  lhs = stri_extract_first_regex(formula, "^[^\\s~]+")
  rhs = stri_replace_first_regex(formula, "[^~]+~\\s+", "")
  terms = stri_split_regex(rhs, "[\\s\\n]+\\+[\\s\\n]+")
  terms = lapply(terms, sort)
  new_formula = mapply(function(lhs, terms) {paste(lhs, "~", paste(terms, collapse = " + "))}, lhs, terms)
  new_formula = stri_replace_all_regex(new_formula, "[\\s\\n]+", " ")
  new_formula = stri_replace_all_fixed(new_formula, "+ +s", "+ s")
  new_formula = stri_replace_all_fixed(new_formula, "+ +", "+")
  names(new_formula) <- NULL
  return(stri_trim(new_formula))
}

# Re-fits a gam model, dropping terms that have been selected out
reduce_model <- function(mod, edf_cutoff = 0.001, recursive=TRUE) {
  low_edf_vars = any(pen.edf(mod) < edf_cutoff)
  if(recursive) {
    while(low_edf_vars) {
      mod = update(mod, formula = as.formula(rm_low_edf(mod, edf_cutoff)))
      low_edf_vars = any(pen.edf(mod) < edf_cutoff)
    }
  } else {
    if(low_edf_vars) {
      mod = update(mod, formula = as.formula(rm_low_edf(mod, edf_cutoff)))
    }
  }
  return(mod)
}

# Makes a reduced version of the RHS of a formula
shortform = function(formula) {
  rhs = stri_replace_first_regex(formula, "[^~]+~\\s+", "")
  rhs = stri_replace_all_regex(rhs, "s\\(([^\\,]+)\\,[^\\)]+\\)", "s($1)")
  rhs= stri_replace_all_fixed(rhs, "+ +s", "+ s")
  stri_replace_all_fixed(rhs, "(1 + | + 1)", "")
}

check_vals <- function(x) {
  w = which(is.nan(x) | is.na(x) | is.infinite(x))
  z = x[w]
  names(z) = w
  return(z)
}


get_relative_contribs <- function(gam_model) {

  terms <- attributes(gam_model$terms)
  offset <- attributes(gam_model$terms)$offset

  formula_chr <- as.character(formula(gam_model))
  model_data <- gam_model$model
  offset_name = stri_replace_first_regex(names(model_data)[terms$offset], "offset\\(([^\\)]+)\\)", "$1")
  names(model_data) <- stri_replace_all_regex(names(model_data), "offset\\(([^\\)]+)\\)", "$1")
  gam_model = gam(formula(gam_model), data=model_data, family=gam_model$family, select=FALSE)

  y = model_data[,terms$response]
  preds = predict(gam_model, type="iterms")
  intercept = attributes(preds)$constant
  smooth_pars  <- gam_model$sp
  rhs <- formula_chr[3]
  lhs <- formula_chr[2]
  devs <- map_dbl(terms$term.labels, function(term) {
    # sub_preds <- rowSums(preds[, !stri_detect_fixed(colnames(preds), term), drop=FALSE]) + intercept
    # if(length(offset_name) ==1 ) sub_preds = sub_preds + model_data[, offset]
    # sub_preds <- gam_model$family$linkinv(sub_preds)
    # res <- gam_model$family$dev.resids(y, sub_preds, 1)
    # s <- attr(res, "sign")
    # if (is.null(s))
    #   s <- sign(y - sub_preds)
    # res <- sqrt(pmax(res, 0, na.rm=TRUE)) * s
    # dev <- sum(res^2)
    term = c(terms$term.labels[terms$term.labels != term], offset_name)
    term_regex= paste0("(", paste(term, collapse="|"), ")")
    new_formula <- stri_extract_all_regex(rhs, paste0("(s|offset)\\(", term_regex, "[^\\)]*\\)"))[[1]] %>%
      paste(collapse = " + ") %>%
      {paste(lhs, "~", .)} %>%
      as.formula()
    smooth_par <- smooth_pars[stri_detect_regex(names(smooth_pars), term_regex)]
    smooth_par <- smooth_par[!stri_detect_regex(names(smooth_par), "\\)[23456789]+$")]
    new_model <- gam(formula = new_formula, sp=smooth_par, data=model_data, family=gam_model$family, select=FALSE)
    dev = deviance(new_model)
    return(dev)
  })
  null_formula = paste(lhs, "~ 1")
  if(length(offset_name) ==1 ) null_formula = paste0(null_formula, "+ offset(", offset_name, ")")
  null_model = gam(formula = as.formula(null_formula), data=model_data, family=gam_model$family, select=FALSE)
  null_model_dev = deviance(null_model)
  orig_dev<- deviance(gam_model)
  dev_expl = (devs - orig_dev)
  dev_expl = dev_expl/sum(dev_expl)
  data_frame(term = terms$term.labels, rel_deviance_explained = dev_expl)

}



