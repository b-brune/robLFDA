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


###
### Evaluation for amplitude outliers
###

load("inst/simulation_amplitude_outliers.Rdata")

##
## Contamination of whole subjects:
##

figpath = "C:\\Users\\bbrune\\Documents\\Papers_Repo\\roblfda\\fig"

# In- and out-of-sample RMSE (observations)
tikzDevice::tikz(glue("{figpath}/new_RMSE_whole_subject.tex"), standAlone = TRUE,  width=6, height=3)
res |>
  filter(score_fn_type == "fourier", whole_subject, outlier_type == "amplitude",
         n == 100, number_of_curves == "15-30", opt.h.cov == 0.1 | is.na(opt.h.cov),
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
             group=interaction(outlier_size, Algorithm, number_of_curves, n, Observed),
             color=boxplot_labels)) +
  geom_boxplot() +
  facet_grid(. ~ `$c_1$`, labeller = "label_both") +
  scale_color_manual(values = new_colorcode) +
  xlab("$b$") + ylab("RMSE") +
  theme(legend.position = "bottom") +
  guides(color=guide_legend(title="", nrow=2))
  dev.off()


######

# In- and out-of-sample RMSE (score functions)
tikzDevice::tikz(glue("{figpath}/new_reconstruction_errors_whole_subject.tex"), standAlone = TRUE,  width=6, height=3)
res |>
  filter(score_fn_type == "fourier", whole_subject, outlier_type == "amplitude",
           n == 100, number_of_curves == "15-30", opt.h.cov == 0.1 | is.na(opt.h.cov),
           cutoff_score_outliers == 2 | is.na(cutoff_score_outliers),
         outlying_subject_proportion==0.1) |>
  pivot_longer(c(inlying.1, inlying.2), names_to="Component") |>
  mutate(Component = ifelse(Component == "inlying.1", "k=1", "k=2")) |>
  rename(Algorithm = algorithm) |>
  mutate(outlier_size = ifelse(outlier_size != 1, outlier_size, "none"),
         outlier_size = factor(outlier_size, levels = c("none", 3 ,5, 10, 20))) |>
  ggplot(aes(outlier_size, value,
             group=interaction(outlier_size, Algorithm, number_of_curves, n),
             color=Algorithm)) +
  geom_boxplot() +
  facet_grid(. ~ Component, labeller = "label_both") +
  scale_color_manual(values = algo_colors) +
  # ggtitle("Reconstruction error of score functions",
  #         subtitle="Whole subject contaminated\n Outlying subject proportion: 0.1") +
  theme(legend.position="bottom") +
  xlab("$b$") + ylab("RMSE")
dev.off()


# Analysis of misclassification rates: Bonferroni type procedure, score distance

ws <- res |> filter(algorithm=="robLFDA", whole_subject)
ws <- ws |>
  mutate(
    TPR = unlist(lapply(ws$bonferroni, \(x) TPR(x))),
    FPR = unlist(lapply(ws$bonferroni, \(x) FPR(x))),
    F1 = unlist(lapply(ws$bonferroni, \(x) F1_score(x)))
  )


# TPR, FPR, F1 score for whole subject with amplitude outliers
tikzDevice::tikz(glue("{figpath}/new_outlying_subject_detection.tex"), standAlone = TRUE,  width=6, height=2.5)
ws |>
  filter(score_fn_type == "fourier", whole_subject, outlier_size > 1) |>
  pivot_longer(c("TPR", "FPR", "F1")) |>
  mutate(name = factor(name, levels = c("TPR", "FPR", "F1"))) |>
  mutate(outlier_size = ifelse(outlier_size != 1, outlier_size, "none"),
         outlier_size = factor(outlier_size, levels = c("none", 3 ,5, 10, 20))) |>
  ggplot(aes(outlier_size, value, color=factor(outlying_subject_proportion),
             group=interaction(outlier_size, outlying_subject_proportion))) +
  geom_boxplot(outlier.size=0.5) +
  facet_wrap(. ~ name, scale="free_y") +
  # ggtitle("Outlying subject detection") +
  ylab("Score") +
  #geom_hline(yintercept=0.05, lty="dotted") +
  xlab("$b$") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("0.05" = "lightblue", "0.1" = "blue", "0.2" = "darkblue")) +
  guides(color=guide_legend(title="$c_1$"))
dev.off()

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


##
## Partial contamination
##
tikzDevice::tikz(glue("{figpath}/new_RMSE_non_whole_subject.tex"), standAlone = TRUE,  width=6, height=3)
res |>
  filter(score_fn_type == "fourier", !whole_subject, outlier_type == "amplitude",
         n == 100, number_of_curves == "15-30", outlying_subject_proportion  == 0.1,
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
  xlab("$b$") + ylab("RMSE") +
  theme(legend.position = "bottom") +
  guides(color=guide_legend(title="", nrow=2))
dev.off()

tikzDevice::tikz(glue("{figpath}/new_reconstruction_errors_partial.tex"), standAlone = TRUE,  width=6, height=5)
res |>
  filter(score_fn_type == "fourier", !whole_subject, outlier_type == "amplitude",
         n == 100, number_of_curves == "15-30",
         cutoff_score_outliers == 2 | is.na(cutoff_score_outliers),
         outlying_subject_proportion == 0.1,
        outliers_per_individual == 0.1,
         !(algorithm == "robLFDA" & mean_estimator %in% c("mean / none", "trimmed_mean"))) |>
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
  xlab("$b$") + ylab("RMSE")
dev.off()


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

tikzDevice::tikz(glue("{figpath}/new_outlying_curve_detection.tex"), standAlone = TRUE,  width=6, height=2.5)
ws |>
  filter(score_fn_type == "fourier", !whole_subject, outlier_size > 1, algorithm == "robLFDA", n == 100, number_of_curves == "15-30",
         !(algorithm == "robLFDA" & mean_estimator %in% c("mean / none", "trimmed_mean")),
         outlying_subject_proportion == 0.1
         ) |>
  #select(-TPR, -FPR, - F1) |>
  #rename(TPR = TPR_comb, FPR = FPR_comb, F1 = F1_comb) |>
  pivot_longer(c(TPR, FPR, F1)) |>
  mutate(name = factor(name, levels = c("TPR", "FPR", "F1"))) |>
  ggplot(aes(factor(outlier_size), color=factor(outliers_per_individual), y = value,
             group=interaction(outlier_size, outliers_per_individual)),) +
  geom_boxplot(outlier.size=0.5) +
  facet_wrap(cutoff_score_outliers ~ name, scales="free_y", labeller=)+
  # ggtitle("Outlying curve detection")  +
  scale_color_manual(values = c("0.05" = "lightblue", "0.1" = "blue", "0.2" = "darkblue")) +
  guides(color=guide_legend(title="$c_2$")) +
  theme(legend.position="bottom") +
  ylab("Score") + xlab("$b$") +
  facetted_pos_scales(y = scales)
dev.off()


###
### Evaluation for amplitude outliers
###

load("inst/simulation_structural_outliers.Rdata")

# RMSE (curves)
tikzDevice::tikz(glue("{figpath}/results_meaningless_scorefn_rmse.tex"), standAlone = TRUE,  width=3.7, height=3.7)
p1 = res |>
  filter(score_fn_type == "fourier", whole_subject, outlier_type == "shape_score",
         n == 100, number_of_curves == "15-30", opt.h.cov == 0.05 | is.na(opt.h.cov),
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
  geom_boxplot() +
  scale_color_manual(values = new_colorcode) +
  xlab("Outlying subject proportion") + ylab("RMSE") +
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
    TPR = unlist(lapply(ws$confmat_outlying_subjects_orth_bonferroni, \(x) TPR(x))),
    FPR = unlist(lapply(ws$confmat_outlying_subjects_orth_bonferroni, \(x) FPR(x))),
    F1 = unlist(lapply(ws$confmat_outlying_subjects_orth_bonferroni, \(x) F1_score(x)))
  )

tikzDevice::tikz(glue("{figpath}/results_meaningless_scorefn_scores.tex"), standAlone = TRUE,  width=3, height=2)
p2 = ws |>
  filter(score_fn_type == "fourier", whole_subject, n==100, number_of_curves == "15-30") |>
  pivot_longer(c("TPR", "FPR", "F1")) |>
  mutate(name = factor(name, levels = c("TPR", "FPR", "F1"))) |>
  mutate(outlier_size = ifelse(outlier_size != 1, outlier_size, "none"),
         outlier_size = factor(outlier_size, levels = c("none", 3 ,5, 10, 20))) |>
  ggplot(aes(name, value, color=factor(outlying_subject_proportion),
             group=interaction(name, outlier_size, outlying_subject_proportion, number_of_curves))) +
  geom_boxplot(outlier.size=0.5) +
  ylab("Score") +
  xlab("") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.box="vertical", legend.margin=margin()) +
  scale_color_manual(values = c("0.05" = "lightblue", "0.1" = "blue", "0.2" = "darkblue")) +
  guides(color=guide_legend(title="$c_1$"))
dev.off()

tikzDevice::tikz(glue("{figpath}/results_shape.tex"), standAlone = TRUE,  width=6, height=3)
(p1 + p2)
dev.off()


