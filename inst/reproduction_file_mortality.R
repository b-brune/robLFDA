# Code to reproduce the mortality example:
set.seed(369258147)

library(tidyverse)
library(patchwork)
library(gghighlight)
library(shades)
library(gamlss)
library(tikzDevice)
library(magrittr)
library(glue)
library(RBF)
library(robLFDA)

theme_set(theme_minimal() + theme(panel.grid.minor = element_blank()))

data(mortality)

# Logarithmize the target function
modeling_df$log_y = log(modeling_df$y + 0.1)

plot_country <- function(code, log=FALSE) {
  if (!log) {
    modeling_df |>
      filter(Code == code, Year %% 4 == 0) |>
      ggplot(aes(Age, Total, group=Year, color=Year)) +
      geom_line() +
      ggtitle(code) +
      scale_color_gradient(limits = c(1950, 2023), low="yellow", high="red")
  } else {
    modeling_df |>
      filter(Code == code, Year %% 4 == 0) |>
      ggplot(aes(Age, log_y, group=Year, color=Year)) +
      geom_line() +
      # ggtitle(code) +
      ylab("log(Total)") +
      scale_color_gradient(limits = c(1950, 2023), low="yellow", high="red")
  }
}

p1 = (plot_country("PRT") + ggtitle(""))

# tikzDevice::tikz(glue("{figpath}/mortality_curve_example_1.tex"), width=6, height=2.5, standAlone=TRUE)
print(p1)
# dev.off()


p2 = plot_country("AUT") +
  plot_country("RUS") +
  plot_country("DEUTW") +
  plot_country("JPN") +
  plot_country("AUT", log=TRUE) +
  plot_country("RUS", log=TRUE) +
  plot_country("DEUTW", log=TRUE) +
  plot_country("JPN", log=TRUE) +
  plot_layout(guides = "collect", nrow=2) & theme(legend.position="bottom", plot.title=element_text(size=10))

# tikzDevice::tikz(glue("{figpath}/mortality_curve_examples.tex"), width=6.5, height=3.5, standAlone=TRUE)
print(p2)
# dev.off()

# Model fitting:

## -----
## Step 1: Centering the data
## -----

# rob_fit = backf.rob(
#   log_y ~ s + timepoint,
#   data=as.data.frame(mean_estimation_df),
#   windows=c(0.1, 0.15), type="Huber"
# )

predict_from_rbf_fit <- function(rbf_fit, newdata=NULL, return_functions = FALSE) {

  if (is.null(newdata)) {
    print("Returning fitted values from initial dataset")
    return(fitted(rbf_fit))
  }

  print("Interpolating...")
  all_values = tibble(s = rbf_fit$Xp[, 1],
                      timepoint = rbf_fit$Xp[, 2],
                      gs = rbf_fit$g.matrix[, 1],
                      gt = rbf_fit$g.matrix[, 2]
  )

  time_spline = all_values |> select(timepoint, gt) |> distinct()
  age_spline = all_values |> select(s, gs) |> distinct()
  # with(time_spline, plot(timepoint, gt, type="l"))
  # with(age_spline, plot(s, gs, type="l"))

  time_fn = approxfun(time_spline$timepoint, time_spline$gt)
  age_fn = approxfun(age_spline$s, age_spline$gs)

  if (return_functions) {
    return(
      list(time_fn = time_fn, age_fn = age_fn)
    )
  }

  mean_comp_time = time_fn(newdata$timepoint)
  mean_comp_age = age_fn(newdata$s)


  return(mean_comp_time + mean_comp_age + rbf_fit$alpha)
}

modeling_df$fitted_rob = predict_from_rbf_fit(rob_fit, newdata=modeling_df)
modeling_df$y_centered_rob = with(modeling_df, log_y - fitted_rob)

# Plot the fitted mean function

plot_helper = modeling_df |>
  select(s, timepoint, Age, Year) |>
  distinct()

fns = predict_from_rbf_fit(rob_fit, plot_helper, return_functions = TRUE)

plot_helper = plot_helper |>
  mutate(`$\\nu(s)$` = fns$age_fn(s),
         `$\\delta(t)$` = fns$time_fn(timepoint))

p3 = (plot_helper |>
    select(Age, `$\\nu(s)$`) |>
    distinct() |>
    ggplot(aes(Age, `$\\nu(s)$`)) +
    geom_line()) +
  (plot_helper |>
     select(Year, `$\\delta(t)$`) |>
     distinct() |>
     ggplot(aes(Year, `$\\delta(t)$`)) +
     geom_line()
  )

# tikzDevice::tikz(glue("{figpath}/mean_function_splines.tex"), width=6, height=2, standAlone=TRUE)
print(p3)
# dev.off()

## ------
## Step 2: Estimate the marginal eigenfunctions
## ------

marginal_eigenfn = estimate_phi(
  modeling_df,
  centered_y = modeling_df$y_centered_rob,
  cov_est="MRCD",
  number_of_components="pve",
  pve=0.85
)


marginal_eigenfn_mrct = estimate_phi(
  modeling_df,
  centered_y = modeling_df$y_centered_rob,
  cov_est="MRCT",
  number_of_components="pve",
  pve=0.85
)




# Visualize the eigenfunctions:
p4 = data.frame(marginal_eigenfn$phi_k) |>
  rename(`$\\phi_1$`= phi_1, `$\\phi_2$` = phi_2,
         `$\\phi_3$` = phi_3) |> #, `$\\phi_4$` = phi_4) |>
  mutate(s = sort(unique(modeling_df$s)),
         Age = sort(unique(modeling_df$Age))) |>
  pivot_longer(starts_with("$")) |>
  #  mutate(value=exp(value) - 0.1) |>
  ggplot(aes(Age, value, linetype=name)) +
  geom_line() +
  ylab("$\\phi_i(s)$") +
  guides(linetype=guide_legend("Component")) +
  scale_linetype_manual(values = c("$\\phi_1$"="solid", "$\\phi_2$"="longdash",
                                   "$\\phi_3$"="dotted", "$\\phi_4$"="dotdash"))

tikzDevice::tikz(glue("{figpath}/mortality_fpcs_log.tex"), width=6, height=2, standAlone=TRUE)
print(p4)
dev.off()

plot_mean_plus_component <- function(component, pve="") {
  data.frame(marginal_eigenfn$phi_k) |>
    bind_cols(cbind(s = modeling_df$s |> unique(),
                    m = fns$age_fn(modeling_df$s |> unique()))) |>
    mutate(Age = s * 85) |>
    filter(Age %% 2 == 0) |>
    ggplot() +
    geom_line(aes(Age, m), color="gray") +
    geom_point(aes(Age, m + !!sym(glue("phi_{component}"))), pch="+", size=2) +
    geom_point(aes(Age, m - !!sym(glue("phi_{component}"))), pch="-", size=3) +
    xlab("Age") + ylab(glue("$\\nu(s) \\pm \\phi_{component}(s)$")) +
    ggtitle(glue("PVE: {pve}"))
}

tikzDevice::tikz(glue("{figpath}/components_with_mean.tex"), standAlone=TRUE, width=6, height=2)
p5 = (plot_mean_plus_component(1, 0.596) +
  plot_mean_plus_component(2, 0.215) +
  plot_mean_plus_component(3, 0.067)) & theme(plot.title=element_text(size=11))
print(p5)
dev.off()
## ------
## Step 3: Estimate the marginal eigenfunctions
## ------

scores = estimate_scores(
  modeling_df,
  centered_y = modeling_df$y_centered_rob,
  phi = marginal_eigenfn$phi_k,
  cov_center = marginal_eigenfn$cov_center
)

score_df = design_df |> bind_cols(scores)

k = ncol(marginal_eigenfn$phi_k)

# Fit the non-parametric models using the efpca function:
time_models = lapply(1:k, \(i) {

  X = robLFDA:::data_for_time_dynamic_fitting(design_df, scores[, i])

  efpca(
    X=X,
    opt.h.cov = 0.03,
    rho.param=1e-3,
    ncov=floor(length(unique(design_df$t)) / 2), #length(unique(score_df$t)),
    max.kappa=1e3,
    prediction_grid = sort(unique(design_df$t)),
    cutoff_outliers=2,
    mean_estimator = function(x) median(x, na.rm=TRUE) # mean(x, trim=0.1, na.rm=TRUE)
    ###
  )
})


pred = predict_grid(
  design_df = design_df,
  ss = sort(unique(modeling_df$s)),
  tt = sort(unique(modeling_df$timepoint)),
  subject_ids = unique(design_df$subject_id),
  y = t(matrix(modeling_df$y_centered_rob, nrow = length(sort(unique(modeling_df$s))))),
  models=time_models,
  phi=marginal_eigenfn$phi,
  mean_fn_model=NULL,
  cov_center=marginal_eigenfn$cov_center
)

# Add the mean since predict_grid cannot yet deal with the custom predict function
pred$yhat = pred$yhat + predict_from_rbf_fit(rob_fit, newdata=pred |> rename(timepoint=t))

# Combine the most important parts of the model fitting into a list:
model = list(
  design_df=design_df,
  modeling_df = modeling_df,
  eigenfn = marginal_eigenfn$phi_k,
  lambda_k = marginal_eigenfn$lambda_k,
  scores_fitted_models = time_models,
  y_matrix = t(matrix(modeling_df$y_centered_rob, nrow = length(sort(unique(modeling_df$s)))))

)

##  ----
##  Outlier analysis
##  ----

# Function to detect the outlying score functions:
detect_outlying_scorefunctions <- function(
    model,
    alpha=0.05,
    seed = 1
) {

  set.seed(seed)

  k = length(model$scores_fitted_models)
  score_models = model$scores_fitted_models
  estimated_scores = estimate_scores(modeling_df, modeling_df$y_centered_rob, model$eigenfn, marginal_eigenfn$cov_center)
  lambda_k = model$lambda_k

  # One model/norm per subject
  norms_score = matrix(nrow=length(unique(model$design_df$subject_id)), ncol=k)
  norms_orth = matrix(nrow=length(unique(model$design_df$subject_id)), ncol=k)

  cutoffs_score_bonferroni = vector("numeric", length = k)
  cutoffs_score_separate = vector("numeric", length = k)
  cutoffs_orth_bonferroni = vector("numeric", length = k)
  cutoffs_orth_separate = vector("numeric", length = k)

  for (i in 1:k) {

    X = robLFDA:::data_for_time_dynamic_fitting(model$design_df, estimated_scores[, i])
    pred = fitted(score_models[[i]], X=X, pve=0.95)

    eigen_tmp = eigen(score_models[[i]]$cov.fun)

    ev = eigen_tmp$values/ nrow(eigen_tmp$vectors) # normalization of eigenvalues

    ev_score = ev[1:ncol(pred$xis)]
    ev_orth = ev[ev > 1e-6][(ncol(pred$xis)+1):length(ev[ev > 1e-6])]

    norms_score[, i] = sapply(1:nrow(pred$pred), \(j) mean(pred$pred[j, ]^2)) / lambda_k[i]

    norms_orth[, i] = sapply(1:nrow(pred$pred), \(j) mean(( X$x[[j]] - pred$pred[j, score_models[[i]]$grid %in% X$pp[[j]]] )^2))

    # Monte Carlo simulation of the quantiles
    qs_score = quantile(replicate(10000, sum(ev_score / lambda_k[i] * rchisq(length(ev_score), df=1))), c(1 - alpha/k, 1-alpha))
    qs_orth = quantile(replicate(10000, sum(ev_orth * rchisq(length(ev_orth), df=1))), c(1 - alpha/k, 1-alpha))

    cutoffs_score_bonferroni[i] = qs_score[1]
    cutoffs_score_separate[i] = qs_score[2]

    cutoffs_orth_bonferroni[i] = qs_orth[1]
    cutoffs_orth_separate[i] = qs_orth[2]
  }


  calc_flags <- function(norms, cutoffs_separate, cutoffs_bonferroni) {
    flag_individual = sapply(1:k, \(i) norms[, i] > cutoffs_separate[i])

    list(
      norms = norms,
      cutoffs_bonferroni=cutoffs_bonferroni,
      cutoffs_separate=cutoffs_separate,

      # individual flages
      flag_individual = flag_individual,
      # flag an individual row
      flag_bonferroni = rowSums(norms) > sum(cutoffs_bonferroni),
      # flagged at least one individual
      flag_min1 = rowSums(flag_individual) >= 1,
      # flagged both individuals separately
      flag_both = rowSums(flag_individual) == 2
    )
  }

  score_results = calc_flags(norms_score, cutoffs_score_separate, cutoffs_score_bonferroni)
  orth_results = calc_flags(norms_orth, cutoffs_orth_separate, cutoffs_orth_bonferroni)


  return(list(score = score_results, orth = orth_results))
}


flags = detect_outlying_scorefunctions(model, alpha=0.05, seed=1111)

# Countries flagged by orthogonal distance:
filter(modeling_df, subject_id %in% which(flags$orth$flag_bonferroni)) %$% unique(Code)

# Countries flagged by score distance:
filter(modeling_df, subject_id %in% which(flags$score$flag_bonferroni)) %$% unique(Code)

flags$bonferroni = flags$orth$flag_bonferroni | flags$score$flag_bonferroni

# Distance-distance plot:
dd_plot = tibble(
  score_dist = rowSums(flags$score$norms),
  orth_dist = rowSums(flags$orth$norms),
  cutoff_score = sum(flags$score$cutoffs_bonferroni),
  cutoff_orth = sum(flags$orth$cutoffs_bonferroni),
  code = arrange(modeling_df, subject_id) |> pull(Code) |> unique()
)

p6 = dd_plot |>
  ggplot(aes(score_dist, orth_dist)) +
  geom_point() +
  geom_vline(aes(xintercept=cutoff_score), linetype="dotted") +
  geom_hline(aes(yintercept=cutoff_orth), linetype="dotted") +
  gghighlight(orth_dist > cutoff_orth | score_dist > cutoff_score, label_key = code) +
  xlab("Score distance $SD_i$") + ylab("Orthogonal distance $OD_i$")

tikzDevice::tikz(glue("{figpath}/distance_distance_plot.tex"), standAlone=TRUE, width=6, height=2.5)
print(p6)
dev.off()

# Outlying curves:

flag_any = rowMeans(matrix(unlist(lapply(time_models, \(x) unlist(x[["outlyingness"]]))), ncol=k)) > 2

score_df = score_df |>
  mutate(flag1_clean = unlist(time_models[[1]]$flagged),
         flag2_clean = unlist(time_models[[2]]$flagged),
         flag3_clean = unlist(time_models[[3]]$flagged),
         #flag4_clean = unlist(result$scores_fitted_models[[4]]$flagged),
         flag_clean = flag_any)


which = c(2, 9, 25, 31, 32)

tmp_plot_df = pred |>
  left_join(score_df, by=c("subject_id" = "subject_id", "t" = "t")) |>
  mutate(flag_bonferroni = subject_id %in% flags) |>
  drop_na() |>
  select(Code, flag_bonferroni, subject_id, starts_with("xi"),
         flag1_clean, flag2_clean, flag3_clean, t) |>
  filter(subject_id %in% which) |>
  pivot_longer(starts_with("xi")) |>
  unique() |>
  mutate(type = ifelse(str_detect(name, ".x"), "smoothed", "raw")) |>
  mutate(name = sapply(name, \(x) str_split(x, "\\.")[[1]][[1]]))

labels = labeller(name = c(xi1 = "$\\xi_1(t)$", xi2 = "$\\xi_2(t)$", xi3 = "$\\xi_3(t)$"))

p7 = tmp_plot_df |>
  ggplot() +
  geom_line(aes(t * (2023 - 1953) + 1953, value),
            data = subset(tmp_plot_df, type == "smoothed")) +
  geom_point(aes(t * (2023 - 1953) + 1953, value, color=flag1_clean),
             data = filter(tmp_plot_df, (type == "raw" & name == "xi1")),
             alpha=0.5, size=0.8) +
  geom_point(aes(t * (2023 - 1953) + 1953, value, color=flag2_clean),
             data = filter(tmp_plot_df, (type == "raw" & name == "xi2")),
             alpha=0.5, size=0.8) +
  geom_point(aes(t * (2023 - 1953) + 1953, value, color=flag3_clean),
             data = filter(tmp_plot_df, (type == "raw" & name == "xi3")),
             alpha=0.5, size=0.8) +
  facet_grid(name ~ Code, scales="free_y", labeller=labels) +
  scale_color_manual(values=c("FALSE" = "black", "TRUE" = "red"), labels=c("No", "Yes")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom",
        strip.text.y = element_text(face=c("bold")),
        strip.text.x = element_text(face=c("bold"))) +
  guides(color=guide_legend("Flagged as \noutlying"))  +
  xlab("Year") + ylab("")

tikzDevice::tikz(glue("{figpath}/fitted_score_functions.tex"), standAlone=TRUE, width=6, height=3.5)
print(p7)
dev.off()

# All score functions

score_fn <- pred |> select(subject_id, t, xi1, xi2, xi3) |> unique() |>
  left_join(modeling_df |> select(subject_id, Code) |> unique(), by="subject_id") |>ungroup()

tikz(glue("{figpath}/all_score_fn.tex"), standAlone=TRUE, width=6, height=2)
p8 = (score_fn |>
    ggplot(aes(t * (2023 - 1953) + 1953, xi1, group=interaction(Code, subject_id), color=Code)) +
    geom_line() +
    gghighlight(subject_id %in% which(flags$bonferroni), label_key=Code) +
    xlab("Year") + ylab("$\\xi_1(t)$") +
    scale_color_brewer(palette="Paired")
) + (
  score_fn |>
    ggplot(aes(t * (2023 - 1953) + 1953, xi2, group=interaction(Code, subject_id), color=Code)) +
    geom_line() +
    gghighlight(subject_id %in% which(flags$bonferroni), label_key=Code) +
    xlab("Year") + ylab("$\\xi_2(t)$") +
    scale_color_brewer(palette="Paired")
) + (
  score_fn |>
    ggplot(aes(t * (2023 - 1953) + 1953, xi3, group=interaction(Code, subject_id), color=Code)) +
    geom_line() +
    gghighlight(subject_id %in% which(flags$bonferroni), label_key=Code) +
    xlab("Year") + ylab("$\\xi_3(t)$") +
    scale_color_brewer(palette="Paired")
)
print(p8)
dev.off()


##################
##################


## DDC Table

ddc_df = score_df|>
  mutate(outlyingness1 = unlist(time_models[[1]]$outlyingness),
         outlyingness2 = unlist(time_models[[2]]$outlyingness),
         outlyingness3 = unlist(time_models[[3]]$outlyingness),
         outlyingness = (outlyingness1 + outlyingness2 + outlyingness3) / 3)



# Arrange the countries by mean outlyingness

Code_ordered = ddc_df |>
  group_by(Code) |>
  summarize(ol = mean(outlyingness)) |>
  arrange(ol) |>
  pull(Code)


plot_ddc <- function(which) {
  pp = ddc_df |>
    mutate(flag_bonferroni = subject_id %in% which(flags$bonferroni),
           Code = str_replace(Code, "_", "-"),
           Code = factor(Code, levels=Code_ordered, ordered=TRUE)
    ) |>
    ggplot(aes(
      (t * (2023 - 1951)) + 1951, Code, fill = !! sym(glue("outlyingness{which}"))
    )) +
    geom_tile() +
    geom_tile(data = ddc_df |> filter(!! sym(glue("flag{which}_clean"))),
              mapping = aes((t * (2023 - 1951)) + 1951, Code),
              fill = "transparent", color="black", inherit.aes=FALSE, linewidth=0.3) +
    scale_fill_gradient2(low="white", mid="yellow", high="red", midpoint=1.5) +
    # geom_hline(yintercept = unique(plot_df$Code[plot_df$subject_id %in% which(flags$bonferroni)]),
    #            linetype = "dotted") +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(face="bold")
    ) +
    guides(fill = guide_colorbar(title="Outlyingness", title.position = "left")) +
    xlab("Year")


  countries = unique(filter(ddc_df, subject_id %in% which(flags$bonferroni))$Code)

  for (country in countries) {
    pp = pp + geom_rect(ymin=which(Code_ordered == str_replace(country, "_", "-")) - 0.5,
                        ymax=which(Code_ordered == str_replace(country, "_", "-")) + 0.5,
                        xmin=1950, xmax=2023, color="red", fill="transparent", linetype="dashed", linewidth=0.3)
  }

  return (pp)
}


for (i in c(1,2,3,"")) {
  dev.new()
  # pdf(file=glue("ddc_table_component{i}.pdf"), width=8, height=5)
  if (i != "") {
    print(plot_ddc(i) + ggtitle(glue("Outlying cells for component {i}")))
  } else {
    print(plot_ddc(i) + ggtitle(glue("Outlying cells (average)")))
  }
  # dev.off()
}


pdf(glue("{figpath}/ddc_table.pdf"), width=8, height=5)
p9 = plot_ddc("")
print(p9)
dev.off()

####
#### Reconstruction errors:
####

## Fit non-robust reference model:

mean_fn = mgcv::gam(log_y ~ s(s) + s(timepoint), data=modeling_df)
modeling_df$y_centered_classic = modeling_df$log_y - fitted(mean_fn)


marginal_eigenfn_classic = estimate_phi(
  modeling_df,
  centered_y = modeling_df$y_centered_classic,
  cov_est="classic",
  number_of_components="pve",
  pve=0.85
)

phi_k_classic = marginal_eigenfn_classic$phi_k

scores_classic = estimate_scores(modeling_df, modeling_df$y_centered_classic, phi_k_classic)

score_df2 = design_df |>
  bind_cols(scores_classic)

k = ncol(phi_k_classic)

time_models_classic = robLFDA:::fit_time_dynamics_nonrobust(
  y= t(matrix(modeling_df$y_centered_classic, nrow = length(sort(unique(modeling_df$s))))),
  phi=phi_k_classic,
  design_df=design_df,
  parametric=FALSE
)

pred_classic = robLFDA:::predict_grid_nonrobust(
  ss = unique(modeling_df$s),
  tt = unique(modeling_df$timepoint),
  subject_ids = unique(design_df$subject_id), models = time_models_classic,
  phi = phi_k_classic,
  mean_fn_model = mean_fn
)


##########
# Evaluation of errors:

eval_df = modeling_df |>
  left_join(pred |> select(subject_id, s, t, yhat),
            by=c("subject_id", "s", "timepoint" = "t")) |>
  left_join(pred_classic |> select(subject_id, s, t, yhat) |>
              rename(yhat_classic = yhat),
            by=c("subject_id", "s", "timepoint" = "t"))

export_errors = eval_df |>
  group_by(Code, timepoint) |>
  summarize(rob = mean((log_y - yhat)^2),
            classic = mean((log_y - yhat_classic)^2),
  ) |>
  group_by(Code) |>
  summarize(rob = mean(rob, trim=0.2)*100, classic=mean(classic, trim=0.2)*100) |>
  pivot_longer(c(rob, classic)) |>
  pivot_wider(names_from=Code, values_from=value)


overall_error = eval_df |>
  group_by(Code, timepoint) |>
  summarize(rob = mean((log_y - yhat)^2),
            classic = mean((log_y - yhat_classic)^2),
  ) |>
  ungroup() |>
  summarize(rob = mean(rob, trim=.2), classic=median(classic, trim=.2))




tikz(glue("{figpath}/errors_mortality.tex"), standAlone = TRUE, height=3, width=6)
p10 = eval_df |>
  group_by(Code, timepoint) |>
  summarize(rob = mean((log_y - yhat)^2),
            classic = mean((log_y - yhat_classic)^2),
  ) |>
  group_by(Code) |>
  summarize(robLFDA = mean(rob, trim=.2), classicLFDA=mean(classic, trim=.2)) |>
  pivot_longer(c(robLFDA, classicLFDA)) |>
  ggplot(aes(Code, value, color=name, fill=name)) +
  geom_bar(stat="identity", position="dodge", alpha=0.3) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  ylab("20\\% trimmed mean squared \nreconstruction error") +
  guides(fill=guide_legend("Method"), color="none") +
  geom_hline(yintercept = pull(overall_error, rob), linetype="dashed", color="red")+
  geom_hline(yintercept = pull(overall_error, classic), linetype="dashed", color="blue") +
  scale_color_manual(values=c("classicLFDA"="blue", "robLFDA"="red")) +
  scale_fill_manual(values=c("classicLFDA"="blue", "robLFDA"="red"))+
  annotate("label", x = 34.8, y = pull(overall_error, rob) - 0.003,
           label = round(pull(overall_error, rob), 3), color="red") +
  annotate("label", x = 34.8, y = pull(overall_error, classic) + 0.005,
           label = round(pull(overall_error, classic), 3), color="blue")
print(p10)
dev.off()



