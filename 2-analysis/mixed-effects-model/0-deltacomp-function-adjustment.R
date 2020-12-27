################################################################################################################################################
## ADJUSTED AND/OR ADDITIONAL FUNCTIONS ON DELTACOMP PACKAGE VERSION 0.2.0
## TO ACCOMMODATE LINEAR MIXED-EFFECTS MODEL
################################################################################################################################################

is_lmer_mod <- function(x) {
  return(inherits(x, 'lmerMod'))
}

fit_lm <- function(y_str, X, random_effect=NULL) {
  ## fit linear mixed effects model if random_effect is specified
  if (!is.null(random_effect)) {
    covars <- colnames(X)
    covars <- covars[covars != y_str]
    re_colnames_tmp <- sapply(random_effect, function(i) strsplit(i, '|', fixed=TRUE))
    re_colnames <- as.vector(sapply(re_colnames_tmp, function(i) gsub('\\)', '', i[2])))
    fixed_effects <- covars[covars != re_colnames]
    fe_formula <- paste(fixed_effects, collapse = '+')
    lme_formula <- as.formula(
      paste(paste0(y_str, " ~ ", fe_formula), random_effect, sep = ' + ')
    )
    
    options(warn = 2) # make warnings errors to be caught
    lme_X <- try( {
      lmer(lme_formula, data = X) # as.data.frame(X))
    }, silent = TRUE )
    options(warn = 0) # make warnings warnings again
    
    # check sucessful fit
    if (is_lmer_mod(lme_X)) {
      cat("---\nSummary of the linear model:\n---\n")
      print(summary(lme_X))
      return(lme_X)
    } else {
      warning("### Model fitting unsuccessful, halting function as calculations will be unreliable ###")
      stop(attr(lme_X, "condition")) # return error msg as string
    }
  }
  else {
    lm_formula <- as.formula(paste(y_str, "~ ."))
    
    options(warn = 2) # make warnings errors to be caught
    lm_X <- try( {
      lm(lm_formula, data = X) # as.data.frame(X))
    }, silent = TRUE )
    options(warn = 0) # make warnings warnings again
    
    # check sucessful fit
    if (is_lm_mod(lm_X)) {
      cat("---\nSummary of the linear model:\n---\n")
      print(summary(lm_X))
      return(lm_X)
    } else {
      warning("### Model fitting unsuccessful, halting function as calculations will be unreliable ###")
      stop(attr(lm_X, "condition")) # return error msg as string
    }
  }
}

append_ilr_coords <- function(dataf, comps, psi) {
  n <- nrow(dataf)
  n_c <- length(comps)
  ilr_comps <- compositions::ilr(dataf[, comps], V = psi)
  ilr_comps <- as.data.frame(ilr_comps[1:n, ])
  ilr_names <- paste0("ilr", 1:(n_c - 1))
  colnames(ilr_comps) <- ilr_names
  ## 2020-12-03: Comment this line
  ## intermittent warning message: 
  ## Warning message:
  ##   In cbind(ilr_comps, dataf) :
  ##   number of rows of result is not a multiple of vector length (arg 1)
  ## causing incorrect dataframe dimension and error in the further commands
  ## dataf <- cbind(ilr_comps, dataf)
  dataf <- dplyr::bind_cols(ilr_comps, dataf)
  return(dataf)
}

predict_delta_comps <- function(dataf, # data.frame of data
                                y, # character name of outcome in dataf
                                comps, # character vector of names of compositions in dataf
                                covars = NULL, # character vector of names of covariates (non-comp variables) in dataf
                                deltas = c(0, 10, 20) / (24 * 60), # changes in compositions to be computed pairwise
                                comparisons = c("prop-realloc", "one-v-one")[1],
                                alpha = 0.05,
                                random_effect = NULL, ## character vector of the specified random effects
                                return_model = FALSE ## if TRUE, only returns the fitted model
) {
  comparisons <- get_comp_type(comparisons)
  n <- nrow(dataf)
  n_comp <- length(comps)
  n_delta <- length(deltas)
  n_covar <- ifelse(is.null(covars), 0, length(covars))
  if (any(abs(deltas) > 1)) {
    stop("deltas must be specified as positive and negative proportions of a composition. i.e., values in (-1, 1).")
  }
  dataf <- rm_na_data_rows(dataf, c(y, comps, covars))
  ## 2020-12-03: comment the following code, unused in this code and sometimes give intermittent error:
  ## Error in get_avg_covs(dataf, covars) : 
  ## Covariate misspecification in data: please have all covariates specified as either factors or numeric variables
  # m_cov <- NULL
  # if (n_covar > 0) {
  #   m_cov <- get_avg_covs(dataf, covars)
  # }
  dataf <- standardise_comps(dataf, comps)
  mean_comps <- compositions::mean.acomp(compositions::acomp(dataf[, 
                                                                   comps]), robust = FALSE)
  if (!all.equal(1, sum(mean_comps), tolerance = 1e-05)) 
    stop("Calculated mean composition does not sum to 1")
  sbp <- create_seq_bin_part(n_comp)
  psi <- compositions::gsi.buildilrBase(sbp)
  dataf <- append_ilr_coords(dataf, comps, psi)
  ilr_names <- paste0("ilr", 1:(n_comp - 1))
  X <- dataf[, colnames(dataf) %in% c(y, ilr_names, covars)]
  lm_X <- fit_lm(y, X, random_effect) ## lm_X <- fit_lm(y, X)
  if (return_model) {
    return(lm_X)
  }
  lm_quants <- extract_lm_quantities(lm_X, alpha = alpha)
  delta_mat <- get_delta_mat(deltas, comparisons, comps, mean_comps)
  n_preds <- nrow(delta_mat)
  poss_comps <- get_all_comparison_mat(deltas, comparisons, 
                                       comps, mean_comps)
  m_comps <- matrix(rep(mean_comps, n_preds), nrow = n_preds, 
                    byrow = TRUE)
  m_delta <- m_comps + delta_mat
  m_delta_less_0 <- rowSums(m_delta < 0)
  if (any(m_delta_less_0 > 0)) {
    warning(paste("By using the supplied deltas, there are NEGATIVE compositional", 
                  "values so these predictions are non-sensical.", 
                  "Consider using smaller delta values or proportional changes in compositions."))
  }
  if (!all.equal(rep(1, n_preds), rowSums(m_delta), tolerance = 1e-05)) 
    stop("Calculated mean composition does not sum to 1")
  ilr_means <- compositions::ilr(m_comps, V = psi)
  ilr_delta <- compositions::ilr(m_delta, V = psi)
  attr(ilr_delta, "class") <- NULL
  attr(ilr_means, "class") <- NULL
  x0_star <- get_x0_star(lm_quants$dmX, n_preds, ilr_names, 
                         ilr_delta, ilr_means)
  if (is_lmer_mod(lm_X)) {
    y0_star <- x0_star %*% colMeans(lm_quants$beta_hat[[1]])
  }
  else {
    y0_star <- x0_star %*% lm_quants$beta_hat 
  }
  se_y0_star <- get_se_y0_star(x0_star, lm_quants$s_e, lm_quants$XtX_inv)
  realloc_nms <- get_realloc_nms(comps, comparisons, poss_comps)
  delta_list <- get_pred_deltas(delta_mat, realloc_nms)
  preds <- cbind(as.data.frame(realloc_nms, stringsAsFactors = FALSE), 
                 delta_list, alpha, get_pred_bounds(y0_star, lm_quants$crit_val, 
                                                    se_y0_star, bound = 0), get_pred_bounds(y0_star, 
                                                                                            lm_quants$crit_val, se_y0_star, bound = -1), get_pred_bounds(y0_star, 
                                                                                                                                                         lm_quants$crit_val, se_y0_star, bound = 1))
  colnames(preds) <- c("comp+", "comp-", "delta", "alpha", 
                       "delta_pred", "ci_lo", "ci_up")
  preds$sig <- ifelse(preds$ci_lo <= 0 & preds$ci_up >= 0, 
                      "", "*")
  ret_obj <- preds
  class(ret_obj) <- c(class(preds), "deltacomp_obj")
  attr(ret_obj, "dataf") <- dataf
  attr(ret_obj, "y") <- y
  attr(ret_obj, "comps") <- comps
  attr(ret_obj, "covars") <- covars
  attr(ret_obj, "deltas") <- deltas
  attr(ret_obj, "comparisons") <- comparisons
  attr(ret_obj, "alpha") <- alpha
  attr(ret_obj, "irl_basis") <- psi
  return(ret_obj)
}