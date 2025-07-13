###
#
#
###

library(robLFDA)
data(mortality)


##
## Perform the train-test split
##

set.seed(369258147)

codes = unique(mortality$Code)

train = vector("list", length(codes))
test = vector("list", length(codes))

for (i in seq_along(codes)) {
  tmp = full_data_df |> 
    filter(Code == codes[i])
  
  tp = unique(tmp$timepoint)
  
  test_tp = sample(tp, size = floor(length(tp) * 0.1))
  train_tp = tp[!(tp %in% test_tp)]
  
  train[[i]] = tmp |> filter(timepoint %in% train_tp) |> arrange(timepoint)
  test[[i]] = tmp |> filter(timepoint %in% test_tp) |> arrange(timepoint)
}

modeling_df = do.call("bind_rows", train)
test_df = do.call("bind_rows", test)

# Visualization of the splitted data
bind_rows(modeling_df |> mutate(dataset = "train"),
          test_df |> mutate(dataset = "test")) |>
  arrange(subject_id, timepoint) |>
  ggplot(aes(Year, Code, fill=dataset)) +
  geom_tile() +
  ggtitle("Data split")


design_df = modeling_df |>
  group_by(subject_id, Code, timepoint) |>
  summarize() |>
  ungroup()

design_df <- design_df |> rename(t = timepoint)

modeling_df$log_y = log(modeling_df$y + 0.1)


# ****************************************************************************
#
# Model fitting
# Robust model
#
# ****************************************************************************

# Step 1: Estimation of the mean function
mean_fn_rob = gamlss(log_y ~ cs(s) + cs(timepoint),
                     data=as.data.frame(modeling_df), family = NO())

modeling_df$fitted_rob = predict(mean_fn_rob, newdata=modeling_df)
modeling_df$y_centered_rob = with(modeling_df, log_y - fitted_rob)


# Step 2: Estimate the marginal covariance function
dim_s = length(unique(modeling_df$s))
y_rob = t(matrix(modeling_df$y_centered_rob, nrow=dim_s)) # Hier wird die Matrix erstellt mit der dann weitergerechnet wird

tmp = estimate_phi(
  y=y_rob,
  cov_est="MRCD",
  number_of_components="pve",
  pve=0.90
) 

phi_k_rob = tmp$phi_k
cov_center = tmp$cov_center

y_rob= t(apply(y_rob, 1, "-", cov_center))


# Step 3: Smooth the score functions

scores = estimate_scores(y_rob, phi_k_rob)
score_df = design_df |> bind_cols(scores)

k = ncol(phi_k_rob)
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


# Create "result" object:
result = list(
  mean_model = mean_fn_rob,
  eigenfn = phi_k_rob,
  lambda_k = tmp$lambda_k,
  cov_center = cov_center,
  scores_raw = scores,
  scores_fitted_models = time_models,
  y = y_rob,
  ss = sort(unique(modeling_df$s)),
  tt = sort(unique(modeling_df$timepoint)),
  design_df = design_df |> ungroup()
)

class(result) = "robLFDAmod"

pred = predict(result)


# ****************************************************************************
#
# Model fitting
# Non-robust model
#
# ****************************************************************************

y_classic = t(matrix(modeling_df$y_centered_rob, nrow=dim_s)) # Hier wird die Matrix erstellt mit der dann weitergerechnet wird

tmp_classic = estimate_phi(
  y=y_classic,
  cov_est="classic",
  number_of_components="pve",
  pve=0.90
) 

phi_k_classic = tmp_classic$phi_k

scores_classic = estimate_scores(y_classic, phi_k_classic)

score_df2 = design_df |> 
  bind_cols(scores_classic)

k = ncol(phi_k_classic)

time_models_classic = fit_time_dynamics_nonrobust(
  y=y_classic,
  phi=phi_k_classic,
  design_df=design_df,
  parametric=FALSE
)

result_classic = list(
  mean_model = mean_fn_rob,
  eigenfn = phi_k_classic,
  scores_raw = scores_classic,
  scores_fitted_models = time_models_classic,
  design_df = design_df,
  ss = sort(unique(modeling_df$s)),
  tt = sort(unique(modeling_df$timepoint))
)


class(result_classic) = "LFDAmod"

pred_classic = predict(result_classic)



# ****************************************************************************
#
# Analysis of model fitting results
#
# ****************************************************************************

## Analysis of detected outliers:

## Outlying score functions:

flags = detect_outlying_scorefunctions(result, alpha=0.1, seed=1211)

which(flags$bonferroni)
filter(modeling_df, subject_id %in% which(flags$bonferroni)) %$% unique(Code)

# Outlying observations:

flag_any = rowMeans(matrix(unlist(lapply(result$scores_fitted_models, \(x) unlist(x[["outlyingness"]]))), ncol=ncol(phi_k_rob))) > 1.64

score_df = score_df |>  
  mutate(flag1 = unlist(result$scores_fitted_models[[1]]$flagged),
         flag2 = unlist(result$scores_fitted_models[[2]]$flagged),
         flag3 = unlist(result$scores_fitted_models[[3]]$flagged),
         flag4 = unlist(result$scores_fitted_models[[4]]$flagged),
         flag = flag_any) 

flag_props = score_df |>
  group_by(subject_id) |>
  summarize(across(starts_with("flag"), ~ mean(.x), .names="{.col}_mean"))

score_df = score_df |>
  left_join(flag_props, by="subject_id") |>
  mutate(flag1_clean = ifelse(flag1_mean >= 0.9, FALSE, flag1),
         flag2_clean = ifelse(flag2_mean >= 0.9, FALSE, flag2),
         flag3_clean = ifelse(flag3_mean >= 0.9, FALSE, flag3),
         flag4_clean = ifelse(flag4_mean >= 0.9, FALSE, flag4),
         flag_clean  = ifelse(flag_mean >= 0.9, FALSE, flag)
  )



plot_component_effects(1, round(sqrt(3), 1))
plot_component_effects(2, sqrt(0.02)) 


plot_matrix(result$scores_fitted_models[[1]]$cov.fun)
plot_matrix(result$scores_fitted_models[[2]]$cov.fun)
plot_matrix(result$scores_fitted_models[[3]]$cov.fun)
plot_matrix(result$scores_fitted_models[[4]]$cov.fun)


plot_subjects(which(flags$bonferroni), 1, flags=NULL, data=pred, score_df = score_df)


modeling_df |> filter(Code %in% c("AUT", "DEUTW", "SWE",  "UKR", "HUN")) %$% unique(subject_id)


plot_subjects(c(which(flags$bonferroni), 1, 2, 10, 15, 20), 2, flags=NULL, data=pred, score_df = score_df)

plot_subjects(c(2,9,32), 1, flags=which(flags$bonferroni), data=pred, score_df = score_df)
plot_subjects(c(2,9,32), 2, flags=which(flags$bonferroni), data=pred, score_df = score_df)
plot_subjects(c(2,9,32), 3, flags=which(flags$bonferroni), data=pred, score_df = score_df)
plot_subjects(c(2,9,32), 3, flags=which(flags$bonferroni), data=pred, score_df = score_df)




plot_smoothed_score_fn(1, highlight=TRUE) +  
  plot_smoothed_score_fn(2, highlight=TRUE) +
  plot_smoothed_score_fn(3, highlight=TRUE) +
  plot_smoothed_score_fn(4, highlight=TRUE) +
  plot_layout(guides="collect") & ggtitle("")

plot_smoothed_score_fn(1, highlight=FALSE) + 
  plot_smoothed_score_fn(2, highlight=FALSE) +
  plot_smoothed_score_fn(3, highlight=FALSE) +
  plot_smoothed_score_fn(4, highlight=FALSE) +
  plot_layout(guides="collect") & ggtitle("")

#
# Calculate the errors
#

test_df = test_df |> mutate(log_y = log(y + 0.1))

eval_df = bind_rows(
  modeling_df |> mutate(data = "train"),
  test_df |> mutate(data = "test")
) |>
  left_join(pred |> select(subject_id, s, t, yhat), 
            by=c("subject_id", "s", "timepoint" = "t")) |>
  left_join(pred_classic |> select(subject_id, s, t, yhat) |> rename(yhat_classic = yhat),
            by=c("subject_id", "s", "timepoint" = "t")) 

overall_error = eval_df |>
  group_by(Code, timepoint, data) |>
  summarize(rob = mean((log_y - yhat)^2),
            classic = mean((log_y - yhat_classic)^2),
  ) |>
  filter(data == "test") |>
  ungroup() |>
  summarize(rob = mean(rob, trim=.2), classic=median(classic, trim=.2)) 

# Comparison of errors:
eval_df |>
  group_by(Code, timepoint, data) |>
  summarize(rob = mean((log_y - yhat)^2),
            classic = mean((log_y - yhat_classic)^2),
  ) |>
  group_by(Code, data) |>
  summarize(rob = mean(rob, trim=.2), classic=mean(classic, trim=.2)) |>
  pivot_longer(c(rob, classic)) |>
  filter(data == "test") |>
  ggplot(aes(Code, value, color=name, fill=name)) +
  geom_bar(stat="identity", position="dodge", alpha=0.3) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  ylab("20\\% trimmed mean squared \nreconstruction error") +
  guides(fill=guide_legend("Method"), color="none") +
  geom_hline(yintercept = pull(overall_error, rob), linetype="dashed", color="red")+
  geom_hline(yintercept = pull(overall_error, classic), linetype="dashed", color="blue") +
  scale_color_manual(values=c("classic"="blue", "rob"="red")) +  
  scale_fill_manual(values=c("classic"="blue", "rob"="red"))+
  annotate("label", x = 34.8, y = pull(overall_error, rob) - 0.003, 
           label = round(pull(overall_error, rob), 3), color="red") +
  annotate("label", x = 34.8, y = pull(overall_error, classic) + 0.005, 
           label = round(pull(overall_error, classic), 3), color="blue")

##
## Create DDC Table
##

score_df = score_df |>  
  mutate(outlyingness1 = unlist(result$scores_fitted_models[[1]]$outlyingness),
         outlyingness2 = unlist(result$scores_fitted_models[[2]]$outlyingness),
         outlyingness3 = unlist(result$scores_fitted_models[[3]]$outlyingness),
         outlyingness4 = unlist(result$scores_fitted_models[[3]]$outlyingness),
         outlyingness = (outlyingness1 + outlyingness2 + outlyingness4) / 4)

plot_df = score_df 


# Arrange the countries by mean outlyingness

Code_ordered = plot_df |>
  group_by(Code) |>
  summarize(ol = mean(outlyingness)) |>
  arrange(ol) |> 
  pull(Code)



plot_ddc <- function(which) {
  pp = plot_df |>
    mutate(flag_bonferroni = subject_id %in% which(flags$bonferroni),
           Code = str_replace(Code, "_", "-"),
           Code = factor(Code, levels=Code_ordered, ordered=TRUE)
    ) |>
    ggplot(aes(
      (t * (2023 - 1951)) + 1951, Code, fill = !! sym(glue("outlyingness{which}"))
    )) +
    geom_tile() +
    geom_tile(data = plot_df |> filter(!! sym(glue("flag{which}_clean"))),
              mapping = aes((t * (2023 - 1951)) + 1951, Code),
              fill = "transparent", color="black", inherit.aes=FALSE, linewidth=0.9) +
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
  
  
  countries = unique(filter(plot_df, subject_id %in% which(flags$bonferroni))$Code)
  
  for (country in countries) {
    pp = pp + geom_rect(ymin=which(Code_ordered == str_replace(country, "_", "-")) - 0.5, 
                        ymax=which(Code_ordered == str_replace(country, "_", "-")) + 0.5, 
                        xmin=1950, xmax=2023, color="red", fill="transparent", linetype="dashed")
  }
  
  return (pp)
}


plot_ddc("")



bind_rows(as.data.frame(phi_k_rob), as.data.frame(phi_k_classic)) |>
  mutate(method = rep(c("robust", "classic"), each = nrow(phi_k_rob)),
         index = rep(sort(unique(modeling_df$Age)), 2)) |>
  pivot_longer(starts_with("phi")) |>
  ggplot(aes(index, value, linetype=method, color=method)) +
  geom_line() +
  facet_wrap(name ~ .) +
  ggtitle("Estimated FPCs")

