## Plots for paper draft
library(tidyverse)
library(data.table)
library(patchwork)
library(tikzDevice)
library(ggh4x)
library(glue)

source("inst/evaluation_functions.R")


theme_set(theme_minimal() + theme(panel.grid.minor=element_blank()))
algo_colors = c("classicLFDA" = "red", "robLFDA" = "blue")
new_colorcode = c("robLFDA - in-sample" = "blue", "robLFDA - out-of-sample" = "lightblue",
                  "classicLFDA - in-sample" = "red", "classicLFDA - out-of-sample" = "orange")

sim_path = "C:\\Users\\bbrune\\OneDrive - TÜV AUSTRIA Data Intelligence GmbH\\Robust longitudinal functional data analysis\\Simulation results"

###
### Evaluation for amplitude outliers
###

# load(glue("{sim_path}/20241012_MRCT_experiments.Rdata"))
load(glue("{sim_path}/20241013_MRCT_experiments_preprocessed.Rdata"))

# load(glue("{sim_path}/20241012_MRCT_experiments.Rdata"))
# res1 = res
# load(glue("{sim_path}/20241013_MRCT_experiments_fullcontamination.Rdata"))
# res2 = res
# load(glue("{sim_path}/20241013_MRCT_experiments_formererrors.Rdata"))
# res3 = res
#
# res = bind_rows(res1, res2, res3)
#
# rm(res1, res2, res3)

# res = preprocess_raw_batchtools_results(res)

# res = bind_cols(res, batchtools::unwrap(res |> select(angles)))
# save(res, file = glue("{sim_path}/20241013_MRCT_experiments_preprocessed.Rdata"))


# figpath = "C:\\Users\\bbrune\\Documents\\Papers_Repo\\roblfda\\fig"

figpath = "C:\\Users\\bbrune\\Documents\\PhD_Work\\DISSERTATION\\figures\\roblfda\\"


##
## Contamination of whole subjects:
##

## RECONSTRUCTION ERROR
plot_reconstruction_error <- function(nn, nc, title=TRUE) {
  p <- res |>
    filter(score_fn_type == "fourier", whole_subject, outlier_type == "amplitude",
           n == nn, number_of_curves == nc, opt.h.cov == 0.1 | is.na(opt.h.cov),
           cutoff_score_outliers == 2 | is.na(cutoff_score_outliers)) |>
    rename(`in-sample` = rmse_osFALSE_obsTRUE,
           `out-of-sample` = rmse_osFALSE_obsFALSE,
           `$c_1$` = outlying_subject_proportion,
           Algorithm = algorithm) |>
    pivot_longer(c(`in-sample`, `out-of-sample`), names_to="Observed") |>
    mutate(
      boxplot_labels = glue::glue("{Algorithm} - {Observed}"),
      outlier_size = ifelse(outlier_size != 1, outlier_size, "none"),
      outlier_size = factor(outlier_size, levels = c("none", 3 ,5, 10, 20))
    ) |>
    ggplot(aes(outlier_size, value,
               group=interaction(outlier_size, Algorithm, boxplot_labels,
                                 number_of_curves, n, Observed),
               color=boxplot_labels)) +
    geom_boxplot() +
    facet_grid(. ~ `$c_1$`, labeller = "label_both") +
    scale_color_manual(values = new_colorcode) +
    xlab("$b$") + ylab("$\\widehat{RMSE}(Y_i(\\cdot, t))$") +
    theme(legend.position = "bottom") +
    guides(color=guide_legend(title="", nrow=2)) +
    ylim(0, NA)
  if (title) {
    p <- p + ggtitle(glue("Reconstruction error, OS1, n={nn}, $n_i$={nc}"))
  }

  return(p)
}



tikzDevice::tikz(glue("{figpath}\\rmse_y_n100_sparse.tex"), width=7, height=3.5, standAlone=TRUE, pointsize=11)
plot_reconstruction_error(nn=100, nc="15-30", title = FALSE)
dev.off()

plot_reconstruction_error(nn=100, nc="5-15", title = TRUE)
plot_reconstruction_error(nn=50, nc="15-30", title = TRUE) + ylim(0, 3)
plot_reconstruction_error(nn=50, nc="5-15", title=TRUE) + ylim(0, 5.5)


## Effect of increasing from n=50 to n=100
tikzDevice::tikz(glue("{figpath}\\supp_rmse_y_comparison_n.tex"), width=7, height=3, standAlone=TRUE, pointsize=11)
res |>
  filter(score_fn_type == "fourier", whole_subject, outlier_type == "amplitude",
         number_of_curves == "15-30", opt.h.cov == 0.1 | is.na(opt.h.cov),
         cutoff_score_outliers == 2 | is.na(cutoff_score_outliers),
         algorithm == "robLFDA") |>
  rename(`in-sample` = rmse_osFALSE_obsTRUE,
         `out-of-sample` = rmse_osFALSE_obsFALSE,
         `$c_1$` = outlying_subject_proportion,
         Algorithm = algorithm) |>
  pivot_longer(c(`in-sample`, `out-of-sample`), names_to="Observed") |>
  mutate(
    n = as.factor(n),
    boxplot_labels = glue::glue("{Algorithm} - {Observed}"),
    outlier_size = ifelse(outlier_size != 1, outlier_size, "none"),
    outlier_size = factor(outlier_size, levels = c("none", 3 ,5, 10, 20))
  ) |>
  ggplot(aes(outlier_size, value,
             group=interaction(outlier_size, Algorithm, boxplot_labels,
                               number_of_curves, n, Observed),
             color=boxplot_labels, linetype=n)) +
  geom_boxplot(outlier.size=0.7) +
  facet_grid(. ~ `$c_1$`, labeller = "label_both") +
  scale_color_manual(values = new_colorcode) +
  xlab("$b$") + ylab("$\\widehat{RMSE}(Y_i(\\cdot, t))$") +
  theme(legend.position = "bottom") +
  guides(color=guide_legend(title="", nrow=1))
dev.off()

tikzDevice::tikz(glue("{figpath}\\supp_rmse_y_comparison_ni.tex"), width=7, height=3, standAlone=TRUE, pointsize=11)
res |>
  filter(score_fn_type == "fourier", whole_subject, outlier_type == "amplitude",
         n==100, opt.h.cov == 0.1 | is.na(opt.h.cov),
         cutoff_score_outliers == 2 | is.na(cutoff_score_outliers),
         algorithm == "robLFDA") |>
  rename(`in-sample` = rmse_osFALSE_obsTRUE,
         `out-of-sample` = rmse_osFALSE_obsFALSE,
         `$c_1$` = outlying_subject_proportion,
         Algorithm = algorithm) |>
  pivot_longer(c(`in-sample`, `out-of-sample`), names_to="Observed") |>
  mutate(
    n = as.factor(n),
    boxplot_labels = glue::glue("{Algorithm} - {Observed}"),
    outlier_size = ifelse(outlier_size != 1, outlier_size, "none"),
    outlier_size = factor(outlier_size, levels = c("none", 3 ,5, 10, 20)),
    `$n_i$` = factor(number_of_curves, levels = c("5-15", "15-30"), ordered=TRUE)
  ) |>
  ggplot(aes(outlier_size, value,
             group=interaction(outlier_size, Algorithm, boxplot_labels,
                               number_of_curves, n, Observed),
             color=boxplot_labels, linetype=`$n_i$`)) +
  geom_boxplot(outlier.size=0.7) +
  facet_grid(. ~ `$c_1$`, labeller = "label_both") +
  scale_color_manual(values = new_colorcode) +
  xlab("$b$") + ylab("$\\widehat{RMSE}(Y_i(\\cdot, t))$") +
  theme(legend.position = "bottom") +
  guides(color=guide_legend(title="", nrow=1)) +
  scale_linetype_manual(values = c("5-15" = "solid", "15-30"="dashed"))
dev.off()


## SCORE FN RECONSTRUCTION ERROR

plot_score_fn_reconstruction_error <- function(nn, nc, title=TRUE) {
  p <- res |>
    filter(score_fn_type == "fourier", whole_subject, outlier_type == "amplitude",
           n == nn, number_of_curves == nc, opt.h.cov == 0.1 | is.na(opt.h.cov),
           cutoff_score_outliers == 2 | is.na(cutoff_score_outliers),
           outlying_subject_proportion==0.1) |>
    pivot_longer(c(inlying.1, inlying.2), names_to="Component") |>
    mutate(Component = ifelse(Component == "inlying.1", "k=1", "k=2")) |>
    rename(Algorithm = algorithm) |>
    mutate(outlier_size = ifelse(outlier_size != 1, outlier_size, "none"),
           outlier_size = factor(outlier_size, levels = c("none", 3 ,5, 10, 20))) |>
    ggplot(aes(outlier_size, value,
               group=interaction(outlier_size, Algorithm, n),
               color=Algorithm)) +
    geom_boxplot() +
    facet_grid(. ~ Component, labeller = "label_both") +
    scale_color_manual(values = algo_colors) +
    # ggtitle("Reconstruction error of score functions",
    #         subtitle="Whole subject contaminated\n Outlying subject proportion: 0.1") +
    theme(legend.position="bottom") +
    xlab("$b$") + ylab("RMSE") +
    ylim(0, NA)

  if (title) {
    p <- p + ggtitle(glue("Reconstruction error for score functions, OS1, n={nn}, $n_i$={nc}"))
  }

  return(p)
}


tikzDevice::tikz(glue("{figpath}\\supp_rmse_xi_n100_sparse.tex"), width=7, height=3, standAlone=TRUE, pointsize=11)
plot_score_fn_reconstruction_error(nn=100, nc="15-30", title = FALSE)
dev.off()



plot_score_fn_reconstruction_error(nn=100, nc="5-15", title = TRUE)
plot_score_fn_reconstruction_error(nn=50, nc="15-30", title = TRUE)
plot_score_fn_reconstruction_error(nn=50, nc="5-15", title=TRUE) + ylim(c(0, 4))


res |>
  filter(score_fn_type == "fourier", whole_subject, outlier_type == "amplitude",
          number_of_curves == "15-30", opt.h.cov == 0.1 | is.na(opt.h.cov),
         cutoff_score_outliers == 2 | is.na(cutoff_score_outliers),
         outlying_subject_proportion==0.1, algorithm=="robLFDA") |>
  pivot_longer(c(inlying.1, inlying.2), names_to="Component") |>
  mutate(Component = ifelse(Component == "inlying.1", "k=1", "k=2")) |>
  rename(Algorithm = algorithm) |>
  mutate(outlier_size = ifelse(outlier_size != 1, outlier_size, "none"),
         outlier_size = factor(outlier_size, levels = c("none", 3 ,5, 10, 20))) |>
  ggplot(aes(outlier_size, value,
             group=interaction(outlier_size, Algorithm, n),
             color=Algorithm, linetype=factor(n))) +
  geom_boxplot() +
  facet_grid(. ~ Component, labeller = "label_both") +
  scale_color_manual(values = algo_colors) +
  # ggtitle("Reconstruction error of score functions",
  #         subtitle="Whole subject contaminated\n Outlying subject proportion: 0.1") +
  theme(legend.position="bottom") +
  xlab("$b$") + ylab("RMSE") +
  ylim(0, NA)


res |>
  filter(score_fn_type == "fourier", whole_subject, outlier_type == "amplitude",
         n==100, opt.h.cov == 0.1 | is.na(opt.h.cov),
         cutoff_score_outliers == 2 | is.na(cutoff_score_outliers),
         outlying_subject_proportion==0.1, algorithm=="robLFDA") |>
  pivot_longer(c(inlying.1, inlying.2), names_to="Component") |>
  mutate(Component = ifelse(Component == "inlying.1", "k=1", "k=2")) |>
  rename(Algorithm = algorithm) |>
  mutate(outlier_size = ifelse(outlier_size != 1, outlier_size, "none"),
         outlier_size = factor(outlier_size, levels = c("none", 3 ,5, 10, 20))) |>
  ggplot(aes(outlier_size, value,
             group=interaction(outlier_size,number_of_curves, Algorithm, n),
             linetype=number_of_curves)) +
  geom_boxplot() +
  facet_grid(. ~ Component, labeller = "label_both") +
  # ggtitle("Reconstruction error of score functions",
  #         subtitle="Whole subject contaminated\n Outlying subject proportion: 0.1") +
  theme(legend.position="bottom") +
  xlab("$b$") + ylab("RMSE") +
  ylim(0, NA)



## Angles:

plot_angles <- function(nn, nc, title=TRUE) {
  p <- res |>
    filter(score_fn_type == "fourier", whole_subject, outlier_type == "amplitude",
           n == nn, number_of_curves == nc, opt.h.cov == 0.1 | is.na(opt.h.cov),
           cutoff_score_outliers == 2 | is.na(cutoff_score_outliers)) |>
    rename(`in-sample` = rmse_osFALSE_obsTRUE,
           `out-of-sample` = rmse_osFALSE_obsFALSE,
           `$c_1$` = outlying_subject_proportion,
           Algorithm = algorithm,
           `$k=1$` = angles.1,
           `$k=2$` = angles.2
    ) |>
    pivot_longer(c(`$k=1$`, `$k=2$`), names_to="Component") |>
    mutate(
      outlier_size = ifelse(outlier_size != 1, outlier_size, "none"),
      outlier_size = factor(outlier_size, levels = c("none", 3 ,5, 10, 20))
    ) |>
    ggplot(aes(outlier_size, value, color=Algorithm,
               group=interaction(
                 outlier_size, Algorithm,
                 number_of_curves, n))) +
    geom_boxplot() +
    facet_grid(Component ~ `$c_1$`, labeller = "label_both") +
    scale_color_manual(values = algo_colors) +
    xlab("$b$") + ylab("RMSE") +
    theme(legend.position = "bottom") +
    guides(color=guide_legend(title="", nrow=1))

  if (title) {
    p <- p + ggtitle(glue("Angles $\\phi_k$, OS1, n={nn}, $n_i$={nc}"))
  }

  return(p)
}

plot_angles(nn=100, nc="15-30", title = TRUE)
plot_angles(nn=100, nc="5-15", title = TRUE)
plot_angles(nn=50, nc="15-30", title = TRUE)
plot_angles(nn=50, nc="5-15", title=TRUE)


## Outlier detection analysis

# Analysis of misclassification rates: Bonferroni type procedure, score distance
eval_sd <- res |> filter(algorithm=="robLFDA", whole_subject)
eval_sd <- eval_sd |>
  mutate(
    TPR = unlist(lapply(lapply(eval_sd$confmat_outlying_subjects, "[[", "bonferroni"), \(x) TPR(x))),
    FPR = unlist(lapply(lapply(eval_sd$confmat_outlying_subjects, "[[", "bonferroni"), \(x) FPR(x))),
    F1 = unlist(lapply(lapply(eval_sd$confmat_outlying_subjects, "[[", "bonferroni"), \(x) F1_score(x)))
  )

eval_od <- res |> filter(algorithm=="robLFDA", whole_subject)
eval_od <- eval_od |>
  mutate(
    TPR = unlist(lapply(lapply(eval_od$confmat_outlying_subjects_orth, "[[", "bonferroni"), \(x) TPR(x))),
    FPR = unlist(lapply(lapply(eval_od$confmat_outlying_subjects_orth, "[[", "bonferroni"), \(x) FPR(x))),
    F1 = unlist(lapply(lapply(eval_od$confmat_outlying_subjects_orth, "[[", "bonferroni"), \(x) F1_score(x)))
  )


plot_outlier_detection <- function(nn, nc, dataset, title=TRUE, which = "SD") {
  p <- dataset |>
    filter(
      n == nn,
      number_of_curves == nc,
      score_fn_type == "fourier",
      whole_subject,
      outlier_size > 1
    ) |>
    pivot_longer(c("TPR", "FPR", "F1")) |>
    mutate(name = factor(name, levels = c("TPR", "FPR", "F1"))) |>
    mutate(outlier_size = ifelse(outlier_size != 1, outlier_size, "none"),
           outlier_size = factor(outlier_size, levels = c("none", 3 ,5, 10, 20))) |>
    mutate(line = ifelse(name == "FPR", 0.05, NA)) |>
    ggplot(aes(outlier_size, value, color=factor(outlying_subject_proportion),
               group=interaction(outlier_size, outlying_subject_proportion))) +
    geom_boxplot(outlier.size=0.5) +
    geom_hline(aes(yintercept = line), linetype="dotted") +
    facet_wrap(. ~ name) + #, scale="free_y") +
    # ggtitle("Outlying subject detection") +
    ylab("Score") +
    #geom_hline(yintercept=0.05, lty="dotted") +
    xlab("$b$") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    scale_color_manual(values = c("0.05" = "lightblue", "0.1" = "blue", "0.2" = "darkblue")) +
    guides(color=guide_legend(title="$c_1$"))

  if (title) {
    p <- p + ggtitle(glue("Outlier detection based on {which}, OS1, n={nn}, $n_i$={nc}"))
  }

  return(p)
}



tikzDevice::tikz(glue("{figpath}\\outlier_analysis_sd_n100_sparse.tex"), width=7, height=3, standAlone=TRUE, pointsize=11)
plot_outlier_detection(100, "15-30", dataset=eval_sd, which="SD", title=FALSE)
dev.off()

plot_outlier_detection(100, "5-15", dataset=eval_sd, which="SD")
plot_outlier_detection(50, "15-30", dataset=eval_sd, which="SD")
plot_outlier_detection(50, "5-15", dataset=eval_sd, which="SD")


tikzDevice::tikz(glue("{figpath}\\supp_outlier_analysis_od_n100_sparse.tex"), width=7, height=3, standAlone=TRUE, pointsize=11)
plot_outlier_detection(100, "15-30", dataset=eval_od, which="OD", title=FALSE)
dev.off()

plot_outlier_detection(100, "5-15", dataset=eval_od, which="OD")
plot_outlier_detection(50, "15-30", dataset=eval_od, which="OD")
plot_outlier_detection(50, "5-15", dataset=eval_od, which="OD")






## Tables with outlier detection performances

ws |>
  filter(score_fn_type == "fourier", whole_subject, outlier_size > 1) |>
  group_by(outlier_size, outlying_subject_proportion) |>
  summarize(across(c(TPR, FPR, F1), ~ mean(.x, na.rm=TRUE), .names="{.col}_mean"),
            across(c("TPR", "FPR", "F1"), ~ sd(.x, na.rm=TRUE), .names="{.col}_sd")
  )

ws |>
  filter(score_fn_type == "fourier", !whole_subject, outlier_size > 1) |>
  group_by(outlier_size, outlying_subject_proportion) |>
  summarize(across(c(TPR, FPR, F1), ~ mean(.x, na.rm=TRUE), .names="{.col}_mean"),
            across(c("TPR", "FPR", "F1"), ~ sd(.x, na.rm=TRUE), .names="{.col}_sd")
  )


######
#### Partial contamination
######

visualize_partial_contamination <- function(nn, nc, c1, data = res, title=TRUE) {
  p <- data |>
    filter(score_fn_type == "fourier", !whole_subject, outlier_type == "amplitude",
           n == nn, number_of_curves == nc, outlying_subject_proportion  == c1,
           cutoff_score_outliers == 2 | is.na(cutoff_score_outliers),
           !(algorithm == "robLFDA" & mean_estimator %in% c("mean / none", "trimmed_mean")),
           opt.h.cov==0.1 | is.na(opt.h.cov),
           outliers_per_individual == 0.1) |>
    pivot_longer(c(rmse_ocFALSE_osFALSE_obsFALSE,
                   rmse_ocFALSE_osFALSE_obsTRUE,
                   rmse_ocFALSE_osTRUE_obsFALSE,
                   rmse_ocFALSE_osTRUE_obsTRUE), names_to="RMSE_type") |>
    rename(Algorithm = algorithm, `Outlying subject proportion` = outlying_subject_proportion) |>
    mutate(
      Observed = sapply(RMSE_type, \(x) ifelse(grepl("obsTRUE", x, fixed = TRUE), "in-sample", "out-of-sample")),
      boxplot_labels = glue::glue("{Algorithm} - {Observed}"),
      `Partially contaminated` = sapply(RMSE_type, \(x) grepl("osTRUE", x, fixed = TRUE)),
      outlier_size = ifelse(outlier_size != 1, outlier_size, "none"),
      outlier_size = factor(outlier_size, levels = c("none", 3 ,5, 10, 20))) |>
    ggplot(aes(outlier_size, value,
               group=interaction(outlier_size, Algorithm, mean_estimator, n, Observed),
               color=boxplot_labels)) +
    geom_boxplot() +
    facet_wrap(. ~ `Partially contaminated`, labeller = "label_both", scales="free_y") +
    # ggtitle("Comparison of RMSEs, in-and-out-of-sample",
    #         subtitle="$n = 100, ~ n_i\\in\\{15, ...., 30\\}$ \n Proportion of outlying subjects: 0.1 \n Proportion of outliers per individual: 0.1") +
    scale_color_manual(values = new_colorcode) +
    xlab("$b$") + ylab("$\\widehat{RMSE}(Y_i(\\cdot,t))$") +
    theme(legend.position = "bottom") +
    guides(color=guide_legend(title="", nrow=2))

  if (title) {
    p <- p + ggtitle(glue("Reconstruction errors, OS2, n={nn}, $n_i$={nc}, $\\c_1$={c1}"))
  }

  return(p)
}


create_figure_title <- function(figpath, main, n, nc, c1) {
  glue("{figpath}/{main}_n{n}_nc{nc}_c1{c1}.tex")
}

tikzDevice::tikz(glue("{figpath}\\partial_rmse_y_n100_sparse.tex"), width=7, height=3.5, standAlone=TRUE, pointsize=11)
visualize_partial_contamination(nn=100, nc="15-30", c1=0.1, data=res, title=FALSE)
dev.off()

tikzDevice::tikz(glue("{figpath}\\supp_partial_rmse_y_n100_sparse_differentc2.tex"), width=7, height=3, standAlone=TRUE, pointsize=11)
res |>
  filter(score_fn_type == "fourier", !whole_subject, outlier_type == "amplitude",
         n == 100, number_of_curves == "15-30", outlying_subject_proportion  == c1,
         cutoff_score_outliers == 2 | is.na(cutoff_score_outliers), algorithm == "robLFDA",
         opt.h.cov==0.1 | is.na(opt.h.cov)) |>
  pivot_longer(c(rmse_ocFALSE_osFALSE_obsFALSE,
                 rmse_ocFALSE_osFALSE_obsTRUE,
                 rmse_ocFALSE_osTRUE_obsFALSE,
                 rmse_ocFALSE_osTRUE_obsTRUE), names_to="RMSE_type") |>
  rename(Algorithm = algorithm, `Outlying subject proportion` = outlying_subject_proportion) |>
  mutate(
    Observed = sapply(RMSE_type, \(x) ifelse(grepl("obsTRUE", x, fixed = TRUE), "in-sample", "out-of-sample")),
    boxplot_labels = glue::glue("{Algorithm} - {Observed}"),
    `Partially contaminated` = sapply(RMSE_type, \(x) grepl("osTRUE", x, fixed = TRUE)),
    outlier_size = ifelse(outlier_size != 1, outlier_size, "none"),
    outlier_size = factor(outlier_size, levels = c("none", 3 ,5, 10, 20)),
    `$c_2$`= factor(outliers_per_individual)) |>
  ggplot(aes(outlier_size, value,
             group=interaction(outlier_size, Algorithm, outliers_per_individual, mean_estimator, n, Observed),
             color=boxplot_labels, linetype=`$c_2$`)) +
  geom_boxplot(outlier.size=0.7) +
  facet_wrap(. ~ `Partially contaminated`, labeller = "label_both", scales="free_y") +
  # ggtitle("Comparison of RMSEs, in-and-out-of-sample",
  #         subtitle="$n = 100, ~ n_i\\in\\{15, ...., 30\\}$ \n Proportion of outlying subjects: 0.1 \n Proportion of outliers per individual: 0.1") +
  scale_color_manual(values = new_colorcode) +
  xlab("$b$") + ylab("$\\widehat{RMSE}(Y_i(\\cdot,t))$") +
  theme(legend.position = "bottom") +
  guides(color=guide_legend(title="", nrow=1))
dev.off()

for (c1 in c(0.05, 0.1, 0.2)) {
  print(visualize_partial_contamination(nn=100, nc="15-30", c1=c1, data=res))
  print(visualize_partial_contamination(nn=100, nc="5-15", c1=c1, data=res))
  print(visualize_partial_contamination(nn=50, nc="15-30", c1=c1, data=res))
  print(visualize_partial_contamination(nn=50, nc="5-15", c1=c1, data=res))
}




# visualize partiob
visualize_partial_contamination_scorefn <- function(nn, nc, c1, data = res, title=TRUE) {
  p <- data |>
    filter(score_fn_type == "fourier", !whole_subject, outlier_type == "amplitude",
           n == nn, number_of_curves == nc, outlying_subject_proportion  == c1,
           cutoff_score_outliers == 2 | is.na(cutoff_score_outliers),
           !(algorithm == "robLFDA" & mean_estimator %in% c("mean / none", "trimmed_mean")),
           opt.h.cov==0.1 | is.na(opt.h.cov),
           outliers_per_individual == 0.1) |>
    mutate(outlier_size = ifelse(outlier_size != 1, outlier_size, "none"),
           outlier_size = factor(outlier_size, levels = c("none", 3 ,5, 10, 20))) |>
    pivot_longer(c(inlying.1, inlying.2, outlying.1, outlying.2), names_to="Components") |>
    mutate(Component = ifelse(Components == "outlying.1" | Components == "inlying.1", "k=1", "k=2"),
           `Partially contaminated` = sapply(Components, function(x) grepl("out", x))) |>
    rename(Algorithm = algorithm) |>
    ggplot(aes(outlier_size, value,
               group=interaction(outlier_size, Algorithm, number_of_curves, n),
               color=Algorithm)) +
    geom_boxplot() +
    facet_wrap(Component ~ `Partially contaminated`, labeller = "label_both",
               scales = "free") +
    scale_color_manual(values = algo_colors) +
    # ggtitle("Reconstruction error of score functions",
    #         subtitle="Partial contamination\n Outlying subject proportion: 0.1") +
    theme(legend.position="bottom") +
    xlab("$b$") + ylab("$\\widehat{RMSE}(\\xi_{ik},\\hat{\\xi}_{ik})$")

  if (title) {
    p <- p + ggtitle(glue("Reconstruction errors (score functions), OS2, n={nn}, $n_i$={nc}, $\\c_1$={c1}"))
  }

  return(p)
}

tikzDevice::tikz(glue("{figpath}\\partial_rmse_xi_n100_sparse.tex"), width=7, height=6, standAlone=TRUE, pointsize=11)
visualize_partial_contamination_scorefn(100, "15-30", 0.1, data=res, title=FALSE)
dev.off()


visualize_partial_contamination_scorefn(100, "5-15", 0.1, data=res)
visualize_partial_contamination_scorefn(50, "15-30", 0.1, data=res)
visualize_partial_contamination_scorefn(50, "5-15", 0.1, data=res)

visualize_partial_contamination_scorefn(100, "15-30", 0.2, data=res)
visualize_partial_contamination_scorefn(100, "5-15", 0.2, data=res)
visualize_partial_contamination_scorefn(50, "15-30", 0.2, data=res)
visualize_partial_contamination_scorefn(50, "5-15", 0.2, data=res)

visualize_partial_contamination_scorefn(100, "15-30", 0.05, data=res)
visualize_partial_contamination_scorefn(100, "5-15", 0.05, data=res)
visualize_partial_contamination_scorefn(50, "15-30", 0.05, data=res)
visualize_partial_contamination_scorefn(50, "5-15", 0.05, data=res)



### Outlier detection analysis

ws <- res |> filter(algorithm=="robLFDA", !whole_subject)


ws <- ws |>
  mutate(
    TPR = unlist(lapply(ws$confmat_outlying_curves_comb, \(x) TPR(x))),
    FPR = unlist(lapply(ws$confmat_outlying_curves_comb, \(x) FPR(x))),
    F1 = unlist(lapply(ws$confmat_outlying_curves_comb, \(x) F1_score(x))),
    precision = unlist(lapply(ws$confmat_outlying_curves_comb, \(x) precision(x))),
    recall = unlist(lapply(ws$confmat_outlying_curves_comb, \(x) recall(x)))
  )



scales <- list(
  # Here you have to specify all the scales, one for each facet row in your case
  scale_y_continuous(limits = c(0, 1)),
  scale_y_continuous(limits = c(0, 0.05)),
  scale_y_continuous(limits = c(0, 1))
)

tikzDevice::tikz(glue("{figpath}/partial_outlier_analysis_n100_sparse.tex"), standAlone = TRUE,  width=7, height=3, pointsize=11)
ws |>
  filter(score_fn_type == "fourier", !whole_subject, outlier_size > 1, algorithm == "robLFDA", n == 100,
         number_of_curves == "15-30",
         !(algorithm == "robLFDA" & mean_estimator %in% c("mean / none", "trimmed_mean")),
         outlying_subject_proportion == 0.1
         ) |>
  #select(-TPR, -FPR, - F1) |>
  #rename(TPR = TPR_comb, FPR = FPR_comb, F1 = F1_comb) |>
  pivot_longer(c(TPR, FPR, F1)) |>
  mutate(line = ifelse(name == "FPR", 0.05, NA)) |>
  mutate(name = factor(name, levels = c("TPR", "FPR", "F1"))) |>
  ggplot(aes(factor(outlier_size), color=factor(outliers_per_individual), y = value,
             group=interaction(outlier_size, outliers_per_individual, covariance_estimator)),) +
  geom_boxplot(outlier.size=0.7) +
  geom_hline(aes(yintercept = line), linetype="dotted") +
  facet_wrap(. ~ name) + #, scales="free_y", labeller=scales)+
  # ggtitle("Outlying curve detection")  +
  scale_color_manual(values = c("0.05" = "lightblue", "0.1" = "blue", "0.2" = "darkblue")) +
  guides(color=guide_legend(title="$c_2$")) +
  theme(legend.position="bottom") +
  ylab("Score") + xlab("$b$")
dev.off()





###
### Evaluation for structural outliers
###

# RMSE (curves)
tikzDevice::tikz(glue("{figpath}/rmse_y_os3_n100_sparse.tex"), standAlone = TRUE,  width=4, height=3)
res |>
  filter(score_fn_type == "fourier", whole_subject, outlier_type == "shape_score",
         n == 100, number_of_curves == "15-30", opt.h.cov == 0.1 | is.na(opt.h.cov),
         cutoff_score_outliers == 2 | is.na(cutoff_score_outliers)) |>
  rename(`in-sample` = rmse_osFALSE_obsTRUE,
         `out-of-sample` = rmse_osFALSE_obsFALSE,
         `Outlying subject proportion` = outlying_subject_proportion,
         Algorithm = algorithm) |>
  pivot_longer(c(`in-sample`, `out-of-sample`), names_to="Observed") |>
  mutate(
    boxplot_labels = glue("{Algorithm} - {Observed}"),
    outlier_size = ifelse(outlier_size != 1, outlier_size, "none"),
    outlier_size = factor(outlier_size, levels = c("none", 3 ,5, 10, 20))) |>
  ggplot(aes(factor(`Outlying subject proportion`), value,
             group=interaction(`Outlying subject proportion`, outlier_size, Algorithm, number_of_curves, n, Observed),
             color=boxplot_labels)) +
  geom_boxplot(outlier.size=0.5) +
  scale_color_manual(values = new_colorcode) +
  xlab("$c_1$") + ylab("$\\widehat{RMSE}(Y_i(\\cdot,t))$") +
  guides(color=guide_legend("", nrow=2)) +
  theme(legend.position = "bottom", legend.box="vertical", legend.margin=margin())
dev.off()

# Reconstruction errors for the score functions

res |>
  filter(score_fn_type == "fourier", whole_subject, outlier_type == "shape_score",
         n == 100, number_of_curves == "15-30", opt.h.cov == 0.1 | is.na(opt.h.cov),
         cutoff_score_outliers == 2 | is.na(cutoff_score_outliers),
         outlying_subject_proportion==0.1) |>
  pivot_longer(c(inlying.1, inlying.2), names_to="Component") |>
  mutate(Component = ifelse(Component == "inlying.1", "k=1", "k=2")) |>
  rename(Algorithm = algorithm) |>
  mutate(outlier_size = ifelse(outlier_size != 1, outlier_size, "none"),
         outlier_size = factor(outlier_size, levels = c("none", 3 ,5, 10, 20))) |>
  ggplot(aes(Component, value,
             group=interaction(outlier_size, Algorithm, number_of_curves, n, Component),
             color=Algorithm)) +
  geom_boxplot() +
  scale_color_manual(values = algo_colors) +
  theme(legend.position="bottom") +
  xlab("Component") + ylab("RMSE")


# Outlier detection performance
ws <- res |> filter(algorithm=="robLFDA")
ws <- ws |>
  mutate(
    TPR = unlist(lapply(lapply(ws$confmat_outlying_subjects_orth, "[[", "bonferroni"), \(x) TPR(x))),
    FPR = unlist(lapply(lapply(ws$confmat_outlying_subjects_orth, "[[", "bonferroni"), \(x) FPR(x))),
    F1 = unlist(lapply(lapply(ws$confmat_outlying_subjects_orth, "[[", "bonferroni"), \(x) F1_score(x)))
  )
sd <- ws |>
  mutate(
    TPR = unlist(lapply(lapply(ws$confmat_outlying_subjects, "[[", "bonferroni"), \(x) TPR(x))),
    FPR = unlist(lapply(lapply(ws$confmat_outlying_subjects, "[[", "bonferroni"), \(x) FPR(x))),
    F1 = unlist(lapply(lapply(ws$confmat_outlying_subjects, "[[", "bonferroni"), \(x) F1_score(x)))
  )




scales <- list(
  # Here you have to specify all the scales, one for each facet row in your case
  scale_y_continuous(limits = c(0, 1)),
  scale_y_continuous(limits = c(0, 1)),
  scale_y_continuous(limits = c(0, 1))
)


tikzDevice::tikz(glue("{figpath}/outlier_detection_os3.tex"), standAlone = TRUE,  width=7, height=3, pointsize=11)
(sd |>
    filter(score_fn_type == "fourier", whole_subject, n==100, number_of_curves == "15-30", opt.h.cov==0.1, cutoff_score_outliers==2,
           outlier_type=="shape_score") |>
    pivot_longer(c("TPR", "FPR", "F1")) |>
    mutate(name = factor(name, levels = c("TPR", "FPR", "F1"))) |>
    mutate(line = ifelse(name == "FPR", 0.05, NA)) |>
    mutate(outlier_size = ifelse(outlier_size != 1, outlier_size, "none"),
           outlier_size = factor(outlier_size, levels = c("none", 3 ,5, 10, 20))) |>
    ggplot(aes(factor(outlying_subject_proportion), value, color=factor(outlying_subject_proportion),
               group=interaction(name, outlying_subject_proportion))) +
    geom_boxplot(outlier.size=0.5) +
    ylab("Score") +
    xlab("") +
    theme_minimal() +
    facet_wrap(. ~ name) +
    theme(legend.position = "bottom", legend.box="vertical", legend.margin=margin()) +
    scale_color_manual(values = c("0.05" = "lightblue", "0.1" = "blue", "0.2" = "darkblue")) +
    guides(color=guide_legend(title="$c_1$")) +
    geom_hline(aes(yintercept=line), linetype="dotted") +
    ggtitle("Outlier detection based on $SD_i$"))  +
  (ws |>
  filter(score_fn_type == "fourier", whole_subject, n==100, number_of_curves == "15-30", opt.h.cov==0.1, cutoff_score_outliers==2,
         outlier_type=="shape_score") |>
  pivot_longer(c("TPR", "FPR", "F1")) |>
  mutate(name = factor(name, levels = c("TPR", "FPR", "F1"))) |>
  mutate(line = ifelse(name == "FPR", 0.05, NA)) |>
  mutate(outlier_size = ifelse(outlier_size != 1, outlier_size, "none"),
         outlier_size = factor(outlier_size, levels = c("none", 3 ,5, 10, 20))) |>
  ggplot(aes(factor(outlying_subject_proportion), value, color=factor(outlying_subject_proportion),
             group=interaction(name, outlying_subject_proportion))) +
  geom_boxplot(outlier.size=0.5) +
  ylab("Score") +
  xlab("") +
  theme_minimal() +
  facet_wrap(. ~ name) +
  theme(legend.position = "bottom", legend.box="vertical", legend.margin=margin()) +
  scale_color_manual(values = c("0.05" = "lightblue", "0.1" = "blue", "0.2" = "darkblue")) +
  guides(color=guide_legend(title="$c_1$")) +
  geom_hline(aes(yintercept=line), linetype="dotted") +
  ggtitle("Outlier detection based on $OD_i$")
  )
dev.off()

# tikzDevice::tikz(glue("{figpath}/results_shape.tex"), standAlone = TRUE,  width=6, height=3)
# (p1 + p2)
# dev.off()



res |>
    filter(score_fn_type == "fourier", whole_subject, outlier_type == "shape_score",
           n == 100, number_of_curves == "15-30", opt.h.cov == 0.1 | is.na(opt.h.cov),
           cutoff_score_outliers == 2 | is.na(cutoff_score_outliers)) |>
    rename(`in-sample` = rmse_osFALSE_obsTRUE,
           `out-of-sample` = rmse_osFALSE_obsFALSE,
           `$c_1$` = outlying_subject_proportion,
           Algorithm = algorithm,
           `$k=1$` = angles.1,
           `$k=2$` = angles.2
    ) |>
    pivot_longer(c(`$k=1$`, `$k=2$`), names_to="Component") |>
    mutate(
      outlier_size = ifelse(outlier_size != 1, outlier_size, "none"),
      outlier_size = factor(outlier_size, levels = c("none", 3 ,5, 10, 20))
    ) |>
    ggplot(aes(outlier_size, value, color=Algorithm,
               group=interaction(
                 outlier_size, Algorithm,
                 number_of_curves, n))) +
    geom_boxplot() +
    facet_grid(Component ~ `$c_1$`, labeller = "label_both") +
    scale_color_manual(values = algo_colors) +
    xlab("$b$") + ylab("RMSE") +
    theme(legend.position = "bottom") +
    guides(color=guide_legend(title="", nrow=1))
