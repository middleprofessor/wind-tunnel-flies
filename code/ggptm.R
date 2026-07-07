## ----setup, include=FALSE--------------------------------------------------------------------------------------
# wrangling
library(data.table)
library(stringr)

# plot function
library(ggplot2)
library(ggsci)
library(ggpubr)
library(ggforce)
library(insight)
library(cowplot)



## ----palettes--------------------------------------------------------------------------------------------------
# get some palettes
pal_okabe_ito <- c(
  "#E69F00",
  "#56B4E9",
  "#009E73",
  "#F0E442",
  "#0072B2",
  "#D55E00",
  "#CC79A7"
)
pal_okabe_ito_blue <- pal_okabe_ito[c(5,6,1,2,3,7,4)]
pal_npg <- pal_npg("nrc")(10)
pal_aaas <- pal_aaas("default")(10)
pal_jco <- pal_jco("default")(10)
pal_frontiers <- pal_frontiers("default")(7)


## --------------------------------------------------------------------------------------------------------------
remove_parentheses <- function(x){
  if(substr(x, 1, 1) == "("){
    x <- substr(x, 2, nchar(x))
  }
  if(substr(x, nchar(x), nchar(x)) == ")"){
    x <- substr(x, 1, nchar(x)-1)
  }
  return(x)
}

pretty_pvalues <- function(p){
  p[p >= 0.1] <- round(p[p > 0.1], 2)
  p[p < 0.1 & p >= 0.01] <- round(p[p < 0.1 & p >= 0.01], 3)
  p[p < 0.01 & p >= 0.001] <- round(p[p < 0.01 & p >= 0.001], 4)
  p_string <- paste0("p=",as.character(p))
  p_string[p < 0.001] <- "p<0.001"
  return(p_string)
}

sci_to_10 <- function(n) {
  # https://stackoverflow.com/questions/29785555/in-r-using-scientific-notation-10-rather-than-e
  output <- format(n, scientific = TRUE)
  output <- sub("1e+", "10^", output) #Replace e with 10^
  output <- sub("\\+0?", "", output) #Remove + symbol and leading zeros on expoent, if > 1
  output <- sub("-0?", "-", output) #Leaves - symbol but removes leading zeros on expoent, if < 1
  return(output)
}


## ----ggcheck_the_qq, warning = FALSE---------------------------------------------------------------------------
ggcheck_the_qq = function(m1,
                   line = "robust",
                   n_boot = 200){
  n <- nobs(m1)
  m1_res <- residuals(m1)
  #sigma_m1_res <- sigma(m1)

  normal_qq <- ppoints(n) %>%
    qnorm()
  sample_qq <- m1_res[order(m1_res)]
  
  # mean + sd
  parametric_slope <- sd(sample_qq)
  parametric_intercept <- mean(sample_qq)
  
  # quartiles
  m1_quartiles <- quantile(m1_res, c(.25, .75))
  qnorm_quartiles <- qnorm( c(.25, .75))
  m1_diff <- m1_quartiles[2] - m1_quartiles[1]
  qnorm_diff <- qnorm_quartiles[2] - qnorm_quartiles[1] # = 1.349
  quartile_slope <- m1_diff/qnorm_diff
  quartile_intercept <- median(m1_quartiles) # median of quartiles not quantiles
  
  # robust uses MASS:rlm (default arguments?)
  qq_observed <- data.table(normal_qq = normal_qq,
                            sample_qq = sample_qq)
  m2 <- rlm(sample_qq ~ normal_qq, data = qq_observed)
  robust_intercept <- coef(m2)[1]
  robust_slope <- coef(m2)[2]
  
  # re-sample ribbon
  set.seed(1)
  resample_qq_model <- numeric(n_boot*n)
  Y <- simulate(m1, n_boot)
  fd <- model.frame(m1) %>%
    data.table
  inc <- 1:n
  for(sim_i in 1:n_boot){
    # parametric bound
    fd[, (1) := Y[,sim_i]]
    m1_class <- class(m1)[1]
    if(m1_class == "lm"){
      ff <- lm(formula(m1), data = fd) 
    }
    if(m1_class == "lmerModLmerTest" | m1_class == "lmerMod"){
      ff <- lmer(formula(m1), data = fd)
    }
    y_res <- residuals(ff)
    resample_qq <- y_res[order(y_res)]
    resample_qq_model[inc] <- resample_qq
    inc <- inc + n
    
    # robust bound
    qq_resampled <- data.table(normal_qq = normal_qq,
                              resample_qq = resample_qq)
    m2_resample <- rlm(resample_qq ~ normal_qq, data = qq_resampled)
    
  }

  qq_sim <- data.table(normal_qq = normal_qq,
                       resample_qq_model = resample_qq_model)
  
  qq_ci_model <- qq_sim[, .(median = median(resample_qq_model),
                      lower = quantile(resample_qq_model, 0.025),
                      upper = quantile(resample_qq_model, 0.975)),
                  by = normal_qq]
  m2_boot <- rlm(median ~ normal_qq, data = qq_ci_model)
  robust_intercept_boot <- coef(m2_boot)[1]
  robust_slope_boot <- coef(m2_boot)[2]
 
  ggplot(data = qq_observed,
         aes(x = normal_qq, y = sample_qq)) +
    
    # ribbon
    geom_ribbon(data = qq_ci_model,
                aes(ymin = lower,
                    ymax = upper,
                    y = median,
                    fill = "band"),
                fill = "gray",
                alpha = 0.6) +
    # draw points
    geom_point() +
    
   # robust
    geom_abline(aes(intercept = robust_intercept,
                    slope = robust_slope,
                    color = "robust"),
                show.legend = FALSE,
                size = 0.75) +
    # robust_boot
    # geom_abline(aes(intercept = robust_intercept_boot,
    #                 slope = robust_slope_boot,
    #                 color = "robust boot"),
    #             show.legend = TRUE,
    #             size = 0.75) +
    xlab("Normal Quantiles") +
    ylab("Sample Quantiles") +
    
    scale_color_manual(values = pal_okabe_ito[c(1:2,5:6)]) +
    theme_minimal_grid() +
    NULL
  
}


## --------------------------------------------------------------------------------------------------------------

ggcheck_the_glm_qq = function(m1,
                   n_sim = 250,
                   se = FALSE,
                   normal = FALSE){
  
  simulationOutput <- simulateResiduals(fittedModel = m1, n = n_sim)
  observed = simulationOutput$scaledResiduals %>%
    sort()
  
  m1_data <- insight::get_data(m1)
  n_points <- nrow(m1_data)
  
  q <- n_points + 1
  x <- seq(1/q, 1 - 1/q, by = 1/q)
  theoretical <- qunif(x)
  
  if(normal == TRUE){
    observed <- qnorm(observed)
    theoretical <- qnorm(theoretical)
  }
  
  gg <- ggscatter(data = data.frame(
    Theoretical = theoretical,
    Observed = observed),
    x="Theoretical", 
    y="Observed",
    title = "Quantile-Residual Uniform-QQ Plot"
  ) +
    
    geom_abline(slope = 1, intercept = 0) +
    
    NULL

  if(se == TRUE){
    qr <- matrix(as.numeric(NA), nrow = n_points, ncol = n_sim)
    fake_counts <- simulationOutput$simulatedResponse
    m1_form <- find_formula(m1)$conditional
    m1_y <- find_response(m1)
    m1_model_name <- model_name(m1)
    y_col <- which(names(m1_data) == m1_y)
    for(j in 1:n_sim){
      m1_data[, 1] <- fake_counts[, j]
      if(m1_model_name == "glm"){
        m1_fake <- glm(m1_form, family = "poisson", data = m1_data)
      }
      if(m1_model_name == "negbin"){
        m1_fake <- glm.nb(m1_form, data = m1_data)
      }
      simulated_output_fake <- simulateResiduals(m1_fake, n_sim)
      qr[,j] <- sort(simulated_output_fake$scaledResiduals)
    }
    qq_ci_model <- data.table(
      Theoretical = theoretical,
      median = apply(qr, 1, median),
      lower = apply(qr, 1, quantile, 0.025),
      upper = apply(qr, 1, quantile, 0.975)
    )
    gg <- gg +
      # ribbon
      geom_ribbon(data = qq_ci_model,
                  aes(ymin = lower,
                      ymax = upper,
                      y = median,
                      fill = "band"),
                  fill = "gray",
                  alpha = 0.6)
  }
  
  gg
  return(gg)
}



## ----ggcheck_the_spreadlevel-----------------------------------------------------------------------------------
ggcheck_the_spreadlevel <- function(m1,
                   n_boot = 200){
  n <- nobs(m1)
  m1_res <- residuals(m1)
  m1_scaled <- m1_res/sd(m1_res)
  m1_root <- sqrt(abs(m1_scaled))
  m1_fitted <- fitted(m1)
  
  m2 <- lm(m1_root ~ m1_fitted)
  m2_intercept <- coef(m2)[1]
  m2_slope <- coef(m2)[2]

  plot_data <- data.table(
    m1_res = sqrt(abs(m1_scaled)),
    fitted = m1_fitted
  )
  
    ggplot(data = plot_data,
         aes(x = fitted, y = m1_res)) +
    
    # ribbon
    # geom_ribbon(data = qq_ci_model,
    #             aes(ymin = lower,
    #                 ymax = upper,
    #                 y = median,
    #                 fill = "band"),
    #             fill = "gray",
    #             alpha = 0.6) +
    # # draw points
      geom_point() +
    
      geom_smooth(method = lm) +
   # robust 
    # geom_abline(aes(intercept = robust_intercept,
    #                 slope = robust_slope,
    #                 color = "robust"),
    #             show.legend = TRUE,
    #             size = 0.75) +
    # robust_boot
    # geom_abline(aes(intercept = robust_intercept_boot,
    #                 slope = robust_slope_boot,
    #                 color = "robust boot"),
    #             show.legend = TRUE,
    #             size = 0.75) +
    xlab("Fitted") +
    ylab("root abs-scaled-residual") +
    
    
    scale_color_manual(values = pal_okabe_ito[c(1:2,5:6)]) +
    theme_minimal_grid() +
    NULL
}


## ----ggcheck_the_model-----------------------------------------------------------------------------------------
ggcheck_the_model <- function(m1){
  gg1 <- ggcheck_the_qq(m1)
  gg2 <- ggcheck_the_spreadlevel(m1)
  cowplot::plot_grid(gg1, gg2, nrow = 1)
}


## --------------------------------------------------------------------------------------------------------------

create_model_data <- function(
    data,
    response_label = "response_label",
    factor1_label = "factor1_label",
    factor2_label = "factor2_label",
    block_label = "block_label",
    nest_label = "nest_label",
    cov_label = "cov_label",
    two_factors = FALSE,
    include_block = FALSE,
    include_nest = FALSE,
    include_cov = FALSE
){
  y_cols <- c(response_label, factor1_label)
  if(two_factors == TRUE){y_cols <- c(y_cols, factor2_label)}
  if(include_block == TRUE){y_cols <- c(y_cols, block_label)}
  if(include_nest == TRUE){y_cols <- c(y_cols, nest_label)}
  if(include_cov == TRUE){y_cols <- c(y_cols, cov_label)}
  model_data <- data[, y_cols] |>
    data.table()
  
  # make generic columns for easier handling of data
  model_data[, data_set := ifelse(include_nest == TRUE,
                                  "tech_reps",
                                  "exp_reps")]
  model_data[, y := get(response_label)]
  model_data[, factor_1 := get(factor1_label)]
  model_data[, factor_2 := ifelse(two_factors == TRUE,
                                  get(factor2_label),
                                  NA)]
  model_data[, plot_factor := ifelse(two_factors == TRUE,
                                     paste(factor_1, factor_2,
                                           sep = "\n"),
                                     factor_1)]
  model_data[, block_id := ifelse(include_block == TRUE,
                                  get(block_label),
                                  NA)]
  model_data[, nest_id := ifelse(include_nest == TRUE,
                                 get(nest_label),
                                 NA)]
  if(include_cov){
    model_data[, cov_id := get(cov_label)]
  }else{
    model_data[, cov_id := as.numeric(NA)]
  }
  
  # reorder factor levels for 2-factor plots
  levels_1 <- levels(model_data$factor_1) |>
    as.character()
  levels_2 <- levels(model_data$factor_2) |>
    as.character()
  levels_table <- expand.grid(levels_1, levels_2, stringsAsFactors = FALSE) |>
    data.table()
  levels_table[, plot_levels := paste(Var1, Var2, sep = "\n")]
  model_data[, plot_factor := factor(plot_factor,
                                     levels = levels_table[, plot_levels])]
  return(model_data)
}



## --------------------------------------------------------------------------------------------------------------

create_plot_data <- function(m1, ptm){
  gg_data <- get_data(m1) |>
    data.table()

  # create generic columns
  gg_data[, y := get(ptm$response_label)]
  gg_data[, factor_1 := get(ptm$factor1_label) |>
            factor()]
  if(ptm$two_factors == TRUE){
    gg_data[, factor_2 := get(ptm$factor2_label) |>
              factor()]
  }
  gg_data[, plot_factor := factor_1]
  if(ptm$two_factors == TRUE){
    gg_data[, plot_factor := paste(factor_1, factor_2, sep = "\n")]
    # reorder factor levels for 2-factor plots
    levels_1 <- levels(gg_data$factor_1) |>
      as.character()
    levels_2 <- levels(gg_data$factor_2) |>
      as.character()
    levels_table <- expand.grid(levels_1, levels_2, stringsAsFactors = FALSE) |>
      data.table()
    levels_table[, plot_levels := paste(Var1, Var2, sep = "\n")]
    gg_data[, plot_factor := factor(plot_factor,
                                    levels = levels_table[, plot_levels])]
  }
  gg_data[, plot_factor_id := as.integer(plot_factor) |>
            as.character()]
  return(gg_data)
}



## --------------------------------------------------------------------------------------------------------------
create_emm_data <- function(m1_emm, ptm){

  if(is.data.frame(m1_emm) == TRUE){
    gg_emm <- data.table(m1_emm)
  }else{
    gg_emm <- summary(m1_emm) |>
      data.table()
  }
  # create plot_factor column - this will be horizontal axis
  gg_emm[, plot_factor := get(ptm$factor1_label)]
  gg_emm[, factor_1 := get(ptm$factor1_label) |>
            factor()]
  if(ptm$two_factors == TRUE){
    gg_emm[, plot_factor := paste(get(ptm$factor1_label), get(ptm$factor2_label),
                                  sep = "\n")]
    gg_emm[, plot_factor := factor(plot_factor, levels = plot_factor)]
  }
  gg_emm[, plot_factor_id := as.integer(plot_factor) |>
           as.character()]
  # create generic columns for mean and CIs
  if("emmean" %in% names(gg_emm)){ # linear models
    gg_emm[, mean := emmean]
  }
  if("response" %in% names(gg_emm)){ # generalized linear models
    gg_emm[, mean := response]
  }
  if("rate" %in% names(gg_emm)){ # generalized linear models
    gg_emm[, mean := rate]
  }
  if("lower.CL" %in% names(gg_emm)){
    gg_emm[, lo := lower.CL]
    gg_emm[, hi := upper.CL]
  }
  if("asymp.LCL" %in% names(gg_emm)){
    gg_emm[, lo := asymp.LCL]
    gg_emm[, hi := asymp.UCL]
  }
  
  return(gg_emm)
}


## ----combine-contrasts-----------------------------------------------------------------------------------------
combine_contrasts <- function(m1_pairs){
  part_1 <- m1_pairs[[1]]
  part_2 <- m1_pairs[[2]]
  part_1_c <- cbind(
    group1 = part_1[, 2],
    group2 = ".",
    part_1[, -2]
  )
  part_2_c <- cbind(
    group1 = ".",
    group2 = part_2[, 2],
    part_2[, -2]
  )
  m1_pairs_c <- rbind(part_1_c, part_2_c)
  colnames(m1_pairs_c)[1:2] <- c(names(part_1)[2], names(part_2)[2])
  return(m1_pairs_c)
}


## ----create-pairs-data-----------------------------------------------------------------------------------------
create_pairs_data <- function(m1_pairs,
                              hide_pairs, # the rows to hide
                              ptm){
  if(is.data.frame(m1_pairs) == TRUE){
    gg_pairs <- data.table(m1_pairs)
  }else{
    gg_pairs <- summary(m1_pairs) %>%
      data.table()
  }
  if(!any(is.na(hide_pairs))){
    gg_pairs <- gg_pairs[-hide_pairs,]    
  }
  
  # is the contrast a difference or ratio?
  contrast_is <- "difference"
  if("ratio" %in% names(gg_pairs)){
    contrast_is <- "ratio"
  }
  
  # create group1 and group2 columns (the two groups of contrast)
  if(contrast_is == "difference"){
    groups <- unlist(str_split(gg_pairs$contrast, " - "))
  }
  if(contrast_is == "ratio"){
    groups <- unlist(str_split(gg_pairs$contrast, " / "))
  }
  groups <- lapply(groups, remove_parentheses) |>
    unlist()
  i_seq <- 1:length(groups)
  gg_pairs[, group1_label := groups[i_seq%%2 != 0]]
  gg_pairs[, group2_label := groups[i_seq%%2 == 0]]
  if(ptm$simple == TRUE){ # works w both pairwise and revpairwise
    simple_group_1 <- names(gg_pairs)[1]
    simple_group_2 <- names(gg_pairs)[2]
    gg_pairs[get(simple_group_1) != ".", group1_label := paste(group1_label, get(simple_group_1))]
    gg_pairs[get(simple_group_1) != ".", group2_label := paste(group2_label, get(simple_group_1))]
    gg_pairs[get(simple_group_2) != ".", group1_label := paste(get(simple_group_2), group1_label)]
    gg_pairs[get(simple_group_2) != ".", group2_label := paste(get(simple_group_2), group2_label)]
  }
  if(ptm$two_factors == TRUE){
    gg_pairs[, group1_label := str_replace(group1_label, " ", "\n")]
    gg_pairs[, group2_label := str_replace(group2_label, " ", "\n")]
  }
  gg_pairs[, group1 := match(group1_label, ptm$plot_factor_levels)]
  gg_pairs[, group2 := match(group2_label, ptm$plot_factor_levels)]
  gg_pairs[, p.print := pretty_pvalues(p.value)]
  gg_pairs[, p.print := format_p(p.value, whitespace = FALSE)]
  
  return(gg_pairs)
}



## --------------------------------------------------------------------------------------------------------------
create_nest_data <- function(m1, gg_data, ptm){
  gg_nest_data <- gg_data[, .(y = mean(get(ptm$response_label), na.rm = TRUE)),
                          by = c(ptm$factor1_label, ptm$factor1_label, "factor_1",
                                 "plot_factor", "plot_factor_id", ptm$nest_id)]
  # for prediction to get y, we need factor nest and block labels
  # convert to original column labels
  setnames(gg_nest_data, old = "y", new = ptm$response_label)
  y_hat <- predict(m1, gg_nest_data, type = "response")
  gg_nest_data[, nest_mean := y_hat]
  return(gg_nest_data)
}


## --------------------------------------------------------------------------------------------------------------
# need to find maximum y-value from experimental reps, technical reps, or CIs
add_y_pos <- function(gg_pairs, gg_data, gg_emm, gg_nest, ptm){
  if(ptm$nested == FALSE | (ptm$nested == TRUE & ptm$show_nest == TRUE)){
    max_data <- max(gg_data[, y], na.rm = TRUE)
    min_data <- min(gg_data[, y], na.rm = TRUE)
  }else{
    max_data <- NA
    min_data <- NA
  }
  if(ptm$nested == TRUE){
    max_nest <- max(gg_nest[, nest_mean], na.rm = TRUE)
    min_nest <- min(gg_nest[, nest_mean], na.rm = TRUE)
  }else{
    max_nest <- NA
    min_nest <- NA
  }
  max_ci <- max(gg_emm[, hi], na.rm = TRUE)
  min_ci <- min(gg_emm[, lo], na.rm = TRUE)
  max_y <- max(c(max_data, max_nest, max_ci), na.rm = TRUE)
  min_y <- min(c(min_data, min_nest, min_ci), na.rm = TRUE)
  p_increment <- 0.06*(max_y - min_y)
  gg_pairs[, y_pos := max_y + .I * p_increment]
  return(gg_pairs)
}




## --------------------------------------------------------------------------------------------------------------
get_ptm_parameters <- function(m1, m1_pairs){
  gg_data <- get_data(m1) |>
    data.table()
  
  ptm <- list()
  ptm$response_label <- find_response(m1)
  predictors <- find_predictors(m1)
  predictor_names <- names(predictors)
  if(predictor_names[1] == "conditional"){
    predictors_fixed <- predictors$conditional
  }
  if(predictor_names[1] == "fixed"){
    predictors_fixed <- predictors$fixed
  }
  ptm$factor1_label <- predictors_fixed[1]
  ptm$factor2_label <- predictors_fixed[2]

  # is factor2 a covariate?
  ptm$covariate <- NA
  ptm$offset <- FALSE
  if(!is.na(ptm$factor2_label)){
    if(is.numeric(gg_data[, get(ptm$factor2_label)])){
      ptm$covariate <- ptm$factor2_label
      ptm$factor2_label <- NA
      
      # is covariate an offset?
      if(any(str_detect(as.character(m1$formula), "offset"))){
        ptm$offset <- TRUE
      }
    }
  }
  ptm$two_factors <- ifelse(is.na(ptm$factor2_label), FALSE, TRUE)
  random <- find_random(m1)$random
  ptm$random <- ifelse(is.null(random), NA, random)
  
  factor_list <- c(ptm$factor1_label, ptm$factor2_label, ptm$random) |>
    na.omit()
  
  # nesting or blocked?
  if(is.na(ptm$random)){
    ptm$nested <- FALSE
    ptm$nest_id <- NA
    ptm$blocked <- FALSE
    ptm$block_id <- NA
  }else{
    counts <- gg_data[!is.na(get(ptm$response_label)), .(N = .N),
                      by = factor_list]
    ptm$nested <- ifelse(any(counts$N > 1), TRUE, FALSE) # could be block block if...
    if(ptm$nest == TRUE){
      ptm$nest_id <- ptm$random
    }else{
      ptm$blocked <- TRUE
      ptm$block_id <- ptm$random
    }
  }
  
  # simple effects?
  if(is.data.frame(m1_pairs) == TRUE){
    gg_pairs <- data.table(m1_pairs)
  }else{
    gg_pairs <- summary(m1_pairs) %>%
      data.table()
  }
  ptm$simple <- ifelse(names(gg_pairs)[1] != "contrast", TRUE, FALSE)
  
  
  return(ptm)
}


## --------------------------------------------------------------------------------------------------------------
ggptm <- function(m1,
                  m1_emm,
                  m1_pairs,
                  hide_pairs = NA, # rows of m1_pairs to hide
                  rescale = 1, # divide y-axis by this amount
                  join_blocks = FALSE,
                  show_nest_data = FALSE,
                  block_id = NA, # this is the column containing the blocks
                  nest_id = NA, # this is the column containing the cluster
                  jitter_spread = 0.8,
                  jitter_width = 0.2,
                  jitter_type = "density", # "none", "density", "jitter"
                  palette = NULL,
                  y_label = NA,
                  y_units = NA,
                  x_axis_labels = NA,
                  font_size = 12
){
  
  # correct m1_pairs if its a list
  if(!is.null(names(m1_pairs[[1]]))){
    m1_pairs <- combine_contrasts(m1_pairs)
  }

  ptm <- get_ptm_parameters(m1, m1_pairs)
  if(is.na(block_id)){
      ptm$blocked <- FALSE
  }else{
       ptm$blocked <- TRUE
  }
  ptm$show_nest <- show_nest_data
  if(!is.na(nest_id)){ptm$nest_id <- nest_id}
  
  gg_data <- create_plot_data(m1, ptm)
  gg_emm <- create_emm_data(m1_emm, ptm)
  ptm$plot_factor_levels <- gg_emm[, plot_factor] |> as.character()
  gg_pairs <- create_pairs_data(m1_pairs, hide_pairs, ptm)
  if(!is.na(nest_id)){
    gg_nest <- create_nest_data(m1, gg_data, ptm)
  }else{
    gg_nest <- NA
  }
  gg_pairs <- add_y_pos(gg_pairs, gg_data, gg_emm, gg_nest, ptm)
  
  if(any(is.na(x_axis_labels)) == TRUE){x_axis_labels <- levels(gg_data$plot_factor)}
  
  # rescale
  gg_data[, y := y / rescale]
  gg_emm[, mean := mean / rescale]
  gg_emm[, lo := lo / rescale]
  gg_emm[, hi := hi / rescale]
  gg_pairs[, y_pos := y_pos / rescale]
  if(!is.na(nest_id)){
    gg_nest[, nest_mean := nest_mean / rescale]
  }
  
  # if cov is offset then create proportion and rescale emm by mean of covarariate to make proportion
  if(ptm$offset == TRUE){
    gg_data[, y := y/get(ptm$covariate) * 100]
    common_scale <- mean(gg_data[, get(ptm$covariate)])
    gg_emm[, mean := mean / common_scale * 100]
    gg_emm[, lo := lo / common_scale * 100]
    gg_emm[, hi := hi / common_scale * 100]
    gg_pairs[, y_pos := y_pos / common_scale * 100]
    if(is.na(y_label)){
#      y_label <- paste0("Relative ", ptm$response_label, " (% of ", ptm$covariate)
      y_label <- paste("%", ptm$response_label)
    }
  }
  
  if(is.na(y_label)){y_label <- ptm$response_label}
  
  gg <- ggplot(data = gg_data,
               aes(x = plot_factor_id,
                   y = y))
  
  # add nested data
  if(ptm$nested == TRUE){
    # nested reps
    if(show_nest_data == TRUE){
      gg <- gg +
        geom_sina(data = gg_data,
                  aes(x = plot_factor_id,
                      y = y,
                      alpha = 1),
                  scale = "width",
                  maxwidth = jitter_width,
                  size = 2,
                  color = "gray",
                  show.legend = FALSE)
    }
  }
  
  # join blocks
  if(ptm$blocked == TRUE & join_blocks == TRUE){
    gg <- gg +
      geom_line(data = gg_data,
                aes(x = plot_factor_id,
                    y = y,
                    group = get(ptm$block_id)),
                position = position_dodge(width = jitter_width),
                color = "grey"
      )
  }
  
  # show points
  if(ptm$nested == FALSE){
    # experimental reps
    if(ptm$blocked == TRUE & join_blocks == TRUE){
      gg <- gg +
        geom_point(data = gg_data,
                   aes(x = plot_factor_id,
                       y = y,
                       color = factor_1,
                       group = get(ptm$block_id)),
                   position = position_dodge(width = jitter_width),
                   size = 3,
                   alpha = .5,
                   show.legend = FALSE)
    }else{
      if(jitter_type == "density"){
        gg <- gg +
          geom_sina(data = gg_data,
                    aes(x = plot_factor_id,
                        y = y,
                        color = factor_1),
                    scale = "width",
                    maxwidth = jitter_width,
                    size = 3,
                    show.legend = FALSE)
      }
      if(jitter_type == "jitter"){
        gg <- gg +
          geom_jitter(data = gg_data,
                      aes(x = plot_factor_id,
                          y = y,
                          color = factor_1),
                      width = jitter_width,
                      size = 3,
                      show.legend = FALSE)
      }
      
    }
  }
  
  # add nest means = experimental reps
  if(ptm$nested == TRUE){
    gg <- gg +
      geom_jitter(data = gg_nest,
                  aes(x = plot_factor_id,
                      y = nest_mean,
                      color = factor_1),
                  width = jitter_width,
                  size = 3,
                  show.legend = FALSE)
  }
  
  
  # add model means and CI
  gg <- gg +
    geom_errorbar(data = gg_emm,
                  aes(x = plot_factor_id,
                      y = mean,
                      ymin = lo,
                      ymax = hi,
                      width =.1,
                      color = factor_1),
                  show.legend = FALSE) +
    geom_point(data = gg_emm,
               aes(x = plot_factor_id,
                   y = mean,
                      color = factor_1),
               size = 4,
               show.legend = FALSE)
  
  # add some color
  # if(palette != "pal_ggplot"){
  #   gg <- gg +
  #     scale_color_manual(values = get(palette))
  # }
  if(!is.null(palette)){
    gg <- gg +
      scale_color_manual(values = palette)
  }
  
  # add p-value brackets
  gg <- gg +
    stat_pvalue_manual(gg_pairs,
                       label = "p.print",
                       y.position = "y_pos",
                       # xmin = "minx",
                       # xmax = "maxx",
                       size = 4,
                       tip.length = 0.01)
  
  # add y-axis label
  if(rescale != 1){
    if(rescale > 1000){
      rescale_str <- sci_to_10(rescale)
    }else{
      rescale_str <- as.character(rescale)
    }
    if(is.na(y_units)){
      y_units <- paste0("X", rescale_str)
    }else{
      y_units <- paste0("X", rescale_str, " ", y_units)
    }
  }
  # if % is in label then
  if(any(str_detect(y_label, "%"))){
    gg <- gg + ylab(y_label)
  }else{
    y_label <- str_replace_all(y_label, " ", "~")
    y_units <- str_replace_all(y_units, " ", "~")
    if(is.na(y_units)){
      gg <- gg +
        ylab(bquote(.(rlang::parse_expr(paste(y_label)))))
    }else{
      gg <- gg +
        ylab(bquote(.(rlang::parse_expr(paste(y_label)))
                    ~(.(rlang::parse_expr(paste(y_units))))))
    }
  }
  
  
  # x-axis tick labels
  gg <- gg + scale_x_discrete(labels = x_axis_labels)
    
  # add theme
  gg <- gg + theme_pubr(base_size = font_size) +
    theme(axis.title.x = element_blank())
  
  # gg
  
  return(gg)
}


## --------------------------------------------------------------------------------------------------------------
plot_response <- function(m1,
                          m1_emm,
                          m1_pairs,
                          hide_pairs = NA, # rows of m1_pairs to hide
                          rescale = 1, # divide y-axis by this amount
                          join_blocks = FALSE,
                          show_nest_data = FALSE,
                          block_id = NA, # this is the column containing the blocks
                          nest_id = NA, # this is the column containing the cluster
                          jitter_spread = 0.8,
                          jitter_width = 0.2,
                          jitter_type = "density", # "none", "density", "jitter"
                          palette = "pal_ggplot",
                          y_label = NA,
                          y_units = NA,
                          x_axis_labels = NA,
                          font_size = 12
){
  return(ggptm(
    m1,
    m1_emm,
    m1_pairs,
    hide_pairs, # rows of m1_pairs to hide
    rescale, # divide y-axis by this amount
    join_blocks,
    show_nest_data,
    block_id, # this is the column containing the blocks
    nest_id, # this is the column containing the cluster
    jitter_spread,
    jitter_width,
    jitter_type, # "none", "density", "jitter"
    palette,
    y_label,
    y_units,
    x_axis_labels,
    font_size    
  ))
}


## --------------------------------------------------------------------------------------------------------------
geom_ancova <- function(m1){
  geom_smooth(method = "lm",
              se = FALSE,
              mapping = aes(y = predict(m1, type = "response"))
  )
}



## --------------------------------------------------------------------------------------------------------------

plot_the_ancova_model <- function(m1, m1_emm, m1_pairs){
  
  ptm <- get_ptm_parameters(m1, m1_pairs)
  ptm$show_nest <- show_nest_data
  if(!is.na(nest_id)){ptm$nest_id <- nest_id}
  gg_data <- create_plot_data(m1, ptm)
  gg_emm <- create_emm_data(m1_emm, ptm)
  ptm$plot_factor_levels <- gg_emm[, plot_factor] |> as.character()
  gg_pairs <- create_pairs_data(m1_pairs, hide_pairs, ptm)
  if(!is.na(nest_id)){
    gg_nest <- create_nest_data(m1, gg_data, ptm)
  }else{
    gg_nest <- NA
  }
  gg_pairs <- add_y_pos(gg_pairs, gg_data, gg_emm, gg_nest, ptm)
  if(is.na(y_label)){y_label <- ptm$response_label}

  gg <- ggplot(data = gg_data,
               aes(x = get(ptm$covariate),
                   y = y,
                   color = plot_factor)) +
    geom_point(size = 2) +
    geom_ancova(m1)

  # add some color
  if(palette != "pal_ggplot"){
    gg <- gg +
      scale_color_manual(values = get(palette))
  }
  
  # add p-value brackets
  # gg <- gg +
  #   stat_pvalue_manual(gg_pairs,
  #                      label = "p.print",
  #                      y.position = "y_pos",
  #                      # xmin = "minx",
  #                      # xmax = "maxx",
  #                      size = 4,
  #                      tip.length = 0.01)
  
  # add x-axis label
  gg <- gg +
    xlab(ptm$covariate)
  
  # add y-axis label
  if(rescale != 1){
    if(rescale > 1000){
      rescale_str <- sci_to_10(rescale)
    }else{
      rescale_str <- as.character(rescale)
    }
    if(is.na(y_units)){
      y_units <- paste0("X", rescale_str)
    }else{
      y_units <- paste0("X", rescale_str, " ", y_units)
    }
  }
  # if % is in label then
  if(any(str_detect(y_label, "%"))){
    gg <- gg + ylab(y_label)
  }else{
    y_label <- str_replace(y_label, " ", "~")
    y_units <- str_replace(y_units, " ", "~")
    if(is.na(y_units)){
      gg <- gg +
        ylab(bquote(.(rlang::parse_expr(paste(y_label)))))
    }else{
      gg <- gg +
        ylab(bquote(.(rlang::parse_expr(paste(y_label)))
                    ~(.(rlang::parse_expr(paste(y_units))))))
    }
  }
  
  # add theme
  gg <- gg + theme_pubr(base_size = font_size)
  

  return(gg)
}




## ----output-as-R-file------------------------------------------------------------------------------------------
# highlight and run to put update into R folder
# knitr::purl("ggptm.Rmd")

