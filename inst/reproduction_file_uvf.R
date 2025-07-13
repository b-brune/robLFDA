# Code to reproduce the PV example and all plots
library(tidyverse)
library(robLFDA)
library(patchwork)
library(glue)
library(RBF)

theme_set(theme_minimal() + theme(panel.grid.minor = element_blank()))

aging_setting_colors = c(
  shades::gradient("viridis", steps=5, space="Lab")
)
names(aging_setting_colors) = c("Alpin1", "Arid1", "Moderate1", "Moderate5", "Tropical2")


figpath = "C:\\Users\\bbrune\\Documents\\Papers_Repo\\roblfda\\fig"


# Load the data:
data(uvf)


# Visualize the modules
plot_module <- function(aging, module) {
  modeling_df |>
    filter(aging_name == aging, name == module) |>
    ggplot(aes(s, y, group=cum_time_h, color=cum_time_h)) +
    geom_line() +
    scale_color_gradient2(low = "black", mid = "yellow", high="red",
                          midpoint=1000, limits=c(0, 3108)) +
    ggtitle(aging) +
    guides(color = guide_colorbar(title = "Exposure \n time [h]")) +
    theme(plot.title=element_text(size=10, face="bold")) +
    xlab("Wavelength [nm]") + ylab("Intensity") #+ ylim(c(0, 28))
}



p1 = (
  plot_module("Alpin1", "INFS003_046") +
    plot_module("Arid1", "INFS003_058") +
    plot_module("Moderate1", "INFS003_055") +
    plot_module("Moderate5", "INFS003_010") +
    plot_module ("Tropical2", "INFS003_033")
) +
  plot_layout(guides = "collect") & theme(legend.position="bottom")

tikzDevice::tikz(glue("{figpath}/example_spectra.tex"), standAlone = TRUE, width=7, height=3.8)
print(p1)
dev.off()

###
### Fitting the robust model:
###

# ...
# Step 1: Mean function
# ...

# Fit like that but due to runtime loaded with data(uvf)
# mean_fn_rob = RBF::backf.rob(
#   y ~ s,
#   data=as.data.frame(modeling_df), windows=0.1,
#   type="Huber"
# )


# Visualize the mean function:
mean_fn_matrix = bind_cols(mean_fn_rob$Xp, mean_fn_rob$g.matrix + mean_fn_rob$alpha)
colnames(mean_fn_matrix) = c("s", "robust")

p2 = mean_fn_matrix |>
  select(s, `robust`) |>
  distinct() |>
  ggplot(aes(s, robust)) +
  geom_line() +
  xlab("$s$") +
  ylab("$\\nu(s)$") +
  guides(linetype=guide_legend(title=""))

tikzDevice::tikz(glue("{figpath}/mean_fn_uvvis.tex"), standAlone = TRUE, width=6, height=3)
print(p2)
dev.off()

# Center the data:
modeling_df$fitted_rob = fitted(mean_fn_rob) # Calculate fitted values for mean function.
modeling_df$y_centered_rob = with(modeling_df, y - fitted_rob) # Center the y-observations


# ...
# Step 2: Marginal eigenfunction
# ...
#
# marginal_eigenfn = estimate_phi(
#   modeling_df,
#   centered_y = modeling_df$y_centered_rob,
#   number_of_components = "pve",
#   pve = 0.95,
#   cov_est = "MRCD")

marginal_eigenfn = estimate_phi(
  modeling_df,
  centered_y = modeling_df$y_centered_rob,
  number_of_components = "pve",
  pve = 0.95,
  cov_est = "MRCT")



k = ncol(marginal_eigenfn$phi)

# Visualize the marginal eigenfunctions
p3 = as.data.frame(marginal_eigenfn$phi) |>
  rename(`$-\\phi_1$`= phi_1, `$-\\phi_2$` = phi_2) |>
  mutate(s = modeling_df$s |> unique()) |>
  pivot_longer(starts_with("$")) |>
  ggplot(aes(s, -value, color=name)) +#, linetype=method)) +
  geom_line() +
  #scale_linetype_manual(values=c("classicLFDA" = "dashed", "robLFDA"  = "solid")) +
  scale_color_manual(values = c("$-\\phi_1$" = "black", "$-\\phi_2$"="gray")) +
  guides(color=guide_legend("Component"), linetype=guide_legend("Method")) +
  xlab("Wavelength [nm]") + ylab("FPC")

tikzDevice::tikz(glue("{figpath}/pv_fpcs_raw.tex"), standAlone = TRUE, width=6, height=2.5)
print(p3)
dev.off()


plot_eigenfn_effect = function(component, title="") {
  as.data.frame(marginal_eigenfn$phi) |>
    bind_cols(distinct(mean_fn_matrix)) |>
    ggplot() +
    geom_line(aes(s, robust), color="gray") +
    geom_point(aes(s, robust + !!sym(glue::glue("phi_{component}"))), pch="+", size=2) +
    geom_point(aes(s, robust - !!sym(glue::glue("phi_{component}"))), pch="-", size=3) +
    xlab("$s$") +
    ylab(glue::glue("$\\nu(s) \\pm \\phi_{component}(s)$")) +
    ggtitle(title) +
    theme(title=element_text(size=9))
}

p4 = plot_eigenfn_effect(1, "PVE: 0.910") +
  plot_eigenfn_effect(2, "PVE: 0.064")

tikzDevice::tikz(glue("{figpath}/pv_fpcs_with_mean.tex"), standAlone = TRUE, width=6, height=2.5)
print(p4)
dev.off()


# ...
# Step 3: Score functions
# ...

# Calculate score proxies
scores = estimate_scores(
  modeling_df,
  centered_y = modeling_df$y_centered_rob,
  phi = marginal_eigenfn$phi,
  cov_center = marginal_eigenfn$cov_center
)

# Create a data frame for model fitting
score_df = design_df |>
  rename(t = timepoint) |>
  bind_cols(scores)

# Specify the parametric models
xi_formulas = list(
  xi1 ~ 1  + ramp + factor(subject_id) * I(t^2),
  xi2 ~ 1  + ramp + factor(subject_id) * I(t^2)
)

# Fit the models
time_models =  lapply(1:k, \(i) robustbase::lmrob(
  xi_formulas[[i]],
  data = score_df,setting="KS2014"))

# ...
# Obtain fitted values
# ...

pred = predict_grid(
  design_df = design_df,
  ss = sort(unique(modeling_df$s)),
  tt = sort(unique(modeling_df$timepoint)),
  subject_ids = unique(design_df$subject_id),
  models=time_models,
  phi=marginal_eigenfn$phi,
  mean_fn_model=NULL,
  cov_center=marginal_eigenfn$cov_center
)

# Hack to obtain the predicted values outside the values used for fitting from
# the estimated RBF mean function:
pred_mean = function(rbf_fit, newdata=NULL, return_functions=FALSE) {

  if (is.null(newdata)) {
    print("Returning fitted values from initial dataset")
    return(fitted(rbf_fit))
  }

  all_values = tibble(s = rbf_fit$Xp[, 1],
                      gs = rbf_fit$g.matrix[, 1])

  freq_spline = all_values |> distinct()
  interp_fn = approxfun(freq_spline$s, freq_spline$gs)

  if (return_functions) {
    return(
      interp_fn
    )
  }

  pred = interp_fn(newdata$s)

  return(pred + rbf_fit$alpha)
}

pred$mean_fn = pred_mean(mean_fn_rob, newdata=pred)
pred$yhat = pred$yhat + pred$mean_fn

# Determine outliers from time models:
score_df$downweighted <- (time_models[[1]]$rweights < 1) | (time_models[[2]]$rweights < 1)


# Calculate prediction error:
errors = modeling_df |>
  left_join(pred |> select(t, subject_id, s, yhat) |> rename(robLFDA = yhat),
            by = c("timepoint" = "t", "subject_id"="subject_id", "s"="s")) |>
  left_join(score_df |> select(subject_id, t, downweighted),
            by = c("timepoint" = "t", "subject_id"="subject_id")) |>
  group_by(subject_id, timepoint, aging_name, downweighted) |>
  summarize(robLFDA = rmse(robLFDA, y))


# Calculate the trimmed mean error:
mean(errors$robLFDA, trim=0.1)


# Plotting the fitted score functions:
group_means = pred |>
  left_join(score_df |> select(subject_id, t, aging_name), by=c("subject_id" = "subject_id", "t" = "t")) |>
  drop_na() |>
  group_by(t, aging_name) |> summarize(xi1_mean = mean(xi1), xi2_mean=mean(xi2))

plot_fitted_score_fn <- function(component) {
  pred |>
    left_join(score_df |> select(subject_id, t, aging_name), by=c("subject_id" = "subject_id", "t" = "t")) |>
    distinct() |>
    drop_na() |>
    mutate(highlight = subject_id %in% c(2, 7, 11, 14, 18)) |>
    ggplot(aes(t, !!sym(glue("xi{component}")), color=aging_name, group=subject_id)) +
    geom_line(aes(t * 3108, !!sym(glue("xi{component}"))), alpha = 0.5, size=.7) +
    geom_line(data = group_means, aes(t * 3108, !!sym(glue("xi{component}_mean")), color=aging_name), inherit.aes=FALSE, size=1) +
    guides(color=guide_legend(title = "Climate \nsetting")) +
    xlab("Exposure time [h]") + ylab(glue("$\\xi_{component}(t)$")) +
    theme(panel.grid.minor = element_blank()) +
    scale_color_manual(values=aging_setting_colors)

  }

p5 = (plot_fitted_score_fn(1) + plot_fitted_score_fn(2)) +
  plot_layout(guides="collect")

tikzDevice::tikz(glue("{figpath}/fitted_score_functions_uvvis.tex"), standAlone = TRUE, width=6, height=2.5)
print(p5)
dev.off()

###
### Fitting the non-robust reference model
###

# Estimate mean function
mean_fn_classic = mgcv::gam(y ~ s(s),
                            dat=as.data.frame(modeling_df))

modeling_df$fitted_classic = fitted(mean_fn_classic)
modeling_df$y_centered_classic = with(modeling_df, y - fitted_classic)

## -------
# Step 2:
## -------

marginal_eigenfn_classic = estimate_phi(
  modeling_df,
  centered_y = modeling_df$y_centered_classic,
  number_of_components="pve",
  pve=0.95,
  cov_est="classic")

k = ncol(marginal_eigenfn_classic$phi)

## -------
# Step 3:
## ------


scores = estimate_scores(
  modeling_df,
  centered_y = modeling_df$y_centered_classic,
  phi = marginal_eigenfn_classic$phi
)

score_df_classic = design_df |> rename(t = timepoint) |> bind_cols(scores)

time_models_classic =  lapply(1:k, \(i) lm(
  xi_formulas[[i]],
  data = score_df_classic |> ungroup())
)


pred_classic = robLFDA:::predict_grid_nonrobust(
  ss = sort(unique(modeling_df$s)),
  tt = sort(unique(modeling_df$timepoint)),
  subject_ids = unique(design_df$subject_id),
  models=time_models_classic,
  phi=marginal_eigenfn_classic$phi,
  mean_fn_model=mean_fn_classic
)

errors_classic = modeling_df |>
  left_join(pred_classic |> select(t, subject_id, s, yhat) |> rename(classicLFDA = yhat),
            by=c("timepoint"="t", "subject_id"="subject_id", "s" ="s")) |>
  left_join(score_df_classic |> select(subject_id, t),
            by = c("timepoint" = "t", "subject_id"="subject_id")) |>
  group_by(subject_id, timepoint, aging_name) |>
  summarize(classicLFDA = rmse(classicLFDA, y))

#Ratio:
mean(errors$robLFDA, trim=0.1) /  mean(errors_classic$classicLFDA, trim=0.1)

mean(errors$robLFDA, trim=0.1)
mean(errors_classic$classicLFDA, trim=0.1)

