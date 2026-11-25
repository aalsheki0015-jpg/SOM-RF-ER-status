# =============================================================================
#  Hybrid Self-Organizing Map (SOM) + Random Forest Framework for
#  ER-Status Classification from DCE-MRI Radiomic Features
#  Dataset: Duke-Breast-Cancer-MRI (TCIA) | Saha et al., Br J Cancer 2018
#  Environment: Posit Cloud / RStudio
#  v3 - real data (synthetic block removed); 14 figures saved at 300 dpi
# =============================================================================

# ---- 1. Setup ---------------------------------------------------------
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/jammy/latest"))

required_pkgs <- c("kohonen", "randomForest", "ggplot2", "dplyr", "pROC",
                   "gridExtra", "readxl")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs) > 0) install.packages(new_pkgs)

library(kohonen)
library(randomForest)
library(ggplot2)
library(dplyr)
library(pROC)
library(gridExtra)
library(readxl)

dir.create("figures", showWarnings = FALSE)

# Helper: render a plotting block to screen AND to a 300-dpi PNG
save_fig <- function(filename, code, width = 3000, height = 2200, res = 300, show = TRUE) {
  expr <- substitute(code)
  env  <- parent.frame()
  if (isTRUE(show)) {
    tryCatch(eval(expr, envir = env),
             error = function(e) message("(screen preview skipped: ", conditionMessage(e), ")"))
  }
  png(file.path("figures", filename), width = width, height = height, res = res)
  eval(expr, envir = env)
  dev.off()
  message("Saved: figures/", filename)
}

# Global settings
SEED     <- 42
SOM_X    <- 12      # 12 x 12 grid (~ 5*sqrt(N) units for N = 922)
SOM_Y    <- 12
SOM_RLEN <- 200
K_SOM    <- 8       # number of SOM super-clusters (revisit after inspecting Fig. 6)
COR_CUT  <- 0.95    # correlation filter threshold
CLIP_Z   <- 5       # z-score clipping

# ---- 2. Real Dataset: Duke-Breast-Cancer-MRI ---------------------------
# Imaging file: sheet "Imaging Features", column 1 = Patient ID, 529 features
# Clinical file: sheet "Data", 3 header rows, column 1 = Patient ID,
#                column 24 = ER status {0 = neg, 1 = pos}
img <- as.data.frame(read_excel("Imaging_Features.xlsx", sheet = "Imaging Features"))
names(img)[1] <- "PatientID"

clin <- read_excel("Clinical_and_Other_Features.xlsx", sheet = "Data",
                   col_names = FALSE, skip = 3, col_types = "text")
clin_df <- data.frame(PatientID = clin[[1]],
                      ER = suppressWarnings(as.integer(clin[[24]])),
                      stringsAsFactors = FALSE)
stopifnot(all(clin_df$ER %in% c(0, 1)))

dat <- merge(img, clin_df, by = "PatientID")
cat("=== Merged data ===\n")
cat("Patients:", nrow(dat), "| Imaging features:", ncol(img) - 1, "\n")

y <- factor(ifelse(dat$ER == 1, "Positive", "Negative"),
            levels = c("Negative", "Positive"))
print(table(ER = y))

# Feature matrix with syntactically safe names (lookup table saved)
orig_names <- setdiff(names(img), "PatientID")
safe_names <- make.names(orig_names, unique = TRUE)
write.csv(data.frame(SafeName = safe_names, OriginalName = orig_names),
          "figures/feature_name_lookup.csv", row.names = FALSE)

X_raw <- suppressWarnings(sapply(dat[, orig_names], as.numeric))
colnames(X_raw) <- safe_names
cat("Missing cells:", sum(is.na(X_raw)), "\n")

# ---- Stratified 80/20 split (done BEFORE preprocessing to avoid leakage) --
strat_split <- function(y, p = 0.8) {
  sort(unlist(lapply(split(seq_along(y), y),
                     function(i) sample(i, floor(p * length(i))))))
}
set.seed(SEED)
train_idx <- strat_split(y)

# ---- Preprocessing (parameters estimated on the training set only) ------
# (a) median imputation
med <- apply(X_raw[train_idx, ], 2, median, na.rm = TRUE)
X_imp <- X_raw
for (j in seq_len(ncol(X_imp))) X_imp[is.na(X_imp[, j]), j] <- med[j]

# (b) remove near-constant features
sd_tr <- apply(X_imp[train_idx, ], 2, sd)
keep  <- which(is.finite(sd_tr) & sd_tr > 1e-8)
X_imp <- X_imp[, keep]

# (c) remove one of each feature pair with |r| > COR_CUT
cm <- abs(cor(X_imp[train_idx, ]))
cm[upper.tri(cm, diag = TRUE)] <- 0
drop_cor <- which(apply(cm, 1, max) > COR_CUT)
X <- if (length(drop_cor) > 0) X_imp[, -drop_cor] else X_imp
feature_cols <- colnames(X)
cat("Features retained after filtering:", length(feature_cols), "\n")

Xdf <- as.data.frame(X)

# ---- 3. SOM Training (descriptive map, all patients, unsupervised) ------
cat("\n=== SOM Training ===\n")
clip <- function(Z) pmin(pmax(Z, -CLIP_Z), CLIP_Z)
X_scaled <- clip(scale(X))

set.seed(SEED)
som_grid  <- somgrid(xdim = SOM_X, ydim = SOM_Y, topo = "hexagonal")
som_model <- som(X_scaled, grid = som_grid, rlen = SOM_RLEN, alpha = c(0.05, 0.01))
n_nodes   <- nrow(getCodes(som_model))
unit_id   <- som_model$unit.classif

# Display features: the 7 features varying most across the SOM codebook
code_var      <- apply(getCodes(som_model), 2, var)
display_feats <- names(sort(code_var, decreasing = TRUE))[1:7]
short <- function(s, n = 30) ifelse(nchar(s) > n, paste0(substr(s, 1, n), "..."), s)

# ---- FIGURE 1 (2x3): univariate distributions --------------------------
hist_cols <- c("#0072B2", "#009E73", "#E69F00", "#56B4E9", "#D55E00")
save_fig("fig01_distributions.png", {
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
  for (i in 1:5) {
    f <- display_feats[i]
    hist(X[, f], col = hist_cols[i], border = "white",
         main = paste0(LETTERS[i], ". ", short(f, 24)), xlab = "Feature value",
         cex.main = 0.9)
  }
  barplot(table(y), col = c("#999999", "#CC79A7"), border = "white",
          main = "F. ER Status", ylab = "Patients", cex.main = 0.9)
  par(mfrow = c(1, 1))
})

# ---- FIGURE 2: correlation heatmap -------------------------------------
cor_mat <- cor(X[, display_feats])
dimnames(cor_mat) <- list(short(display_feats, 22), short(display_feats, 22))
cor_df  <- as.data.frame(as.table(cor_mat))

p_cor <- ggplot(cor_df, aes(Var1, Var2, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", Freq)), size = 3.2) +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00",
                       midpoint = 0, limits = c(-1, 1)) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Figure 2. Correlation Matrix of SOM-Dominant Features",
       x = "", y = "", fill = "r")

save_fig("fig02_correlation.png", { print(p_cor) })

# ---- FIGURE 3 (1x2): SOM convergence + codebook vectors ----------------
save_fig("fig03_som_training.png", {
  par(mfrow = c(1, 2))
  plot(som_model, type = "changes", main = "A. SOM Training Convergence")
  plot(som_model, type = "codes", codeRendering = "lines",
       main = "B. SOM Codebook Vectors")
  par(mfrow = c(1, 1))
}, width = 3000, height = 1600)

# ---- FIGURE 4 (1x2): node counts + U-Matrix -----------------------------
save_fig("fig04_counts_umatrix.png", {
  par(mfrow = c(1, 2))
  plot(som_model, type = "counts",          main = "A. Patients per Node")
  plot(som_model, type = "dist.neighbours", main = "B. U-Matrix")
  par(mfrow = c(1, 1))
}, width = 3000, height = 1600)

# ---- FIGURE 5 (2x4): 7 component planes + ER-positive fraction ----------
pos_frac <- sapply(1:n_nodes, function(i) {
  idx <- which(unit_id == i)
  if (length(idx) == 0) NA else mean(y[idx] == "Positive")
})

save_fig("fig05_component_planes.png", {
  par(mfrow = c(2, 4), mar = c(3, 3, 3, 4))
  for (f in display_feats) {
    plot(som_model, type = "property",
         property = getCodes(som_model)[, f], main = short(f, 26), cex.main = 0.8)
  }
  plot(som_model, type = "property", property = pos_frac,
       main = "ER-Positive Fraction")
  par(mfrow = c(1, 1))
}, width = 3800, height = 2200)

# ---- FIGURE 6 (1x2): super-clusters + patients by true ER status --------
som_cluster <- cutree(hclust(dist(getCodes(som_model))), k = K_SOM)
clu_pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
             "#0072B2", "#D55E00", "#CC79A7", "#999999",
             "#882255", "#44AA99")[1:K_SOM]
er_cols <- c("Negative" = "#D55E00", "Positive" = "#0072B2")

save_fig("fig06_superclusters.png", {
  par(mfrow = c(1, 2))
  plot(som_model, type = "mapping", bgcol = clu_pal[som_cluster], pch = NA,
       main = sprintf("A. SOM Super-Clusters (k = %d)", K_SOM))
  add.cluster.boundaries(som_model, som_cluster)
  plot(som_model, type = "mapping", col = er_cols[as.character(y)],
       pch = 19, cex = 0.5, main = "B. Patients Mapped by ER Status")
  legend("bottomleft", legend = names(er_cols), col = er_cols, pch = 19,
         bty = "n", cex = 0.8)
  par(mfrow = c(1, 1))
}, width = 3000, height = 1600)

# Cluster x ER association
clu_tab <- table(Cluster = som_cluster[unit_id], ER = y)
cat("\n--- SOM cluster x ER ---\n"); print(clu_tab)
print(suppressWarnings(chisq.test(clu_tab)))

# ---- 4. Baseline Random Forest -----------------------------------------
cat("\n=== Baseline Random Forest ===\n")
x_train <- Xdf[train_idx, ];  y_train <- y[train_idx]
x_test  <- Xdf[-train_idx, ]; y_test  <- y[-train_idx]

set.seed(SEED)
rf_model <- randomForest(x = x_train, y = y_train, ntree = 300, importance = TRUE)
print(rf_model)

predictions <- predict(rf_model, newdata = x_test)
conf_mat <- table(Predicted = predictions, Actual = y_test)
accuracy <- sum(diag(conf_mat)) / sum(conf_mat)

cohen_kappa <- function(cm) {
  po <- sum(diag(cm)) / sum(cm)
  pe <- sum(rowSums(cm) * colSums(cm)) / sum(cm)^2
  (po - pe) / (1 - pe)
}
bal_acc <- function(cm) mean(diag(cm) / colSums(cm))   # mean per-class recall
kappa_base <- cohen_kappa(conf_mat)
bacc_base  <- bal_acc(conf_mat)

prob_test <- predict(rf_model, newdata = x_test, type = "prob")
roc_base  <- roc(response = y_test, predictor = prob_test[, "Positive"],
                 levels = c("Negative", "Positive"), direction = "<", quiet = TRUE)
auc_base  <- as.numeric(auc(roc_base))

cat("\n--- Baseline Confusion Matrix ---\n"); print(conf_mat)
cat(sprintf("Baseline Accuracy: %.2f%% | Balanced Acc: %.3f | Kappa: %.3f | AUC: %.3f\n",
            accuracy * 100, bacc_base, kappa_base, auc_base))
cat(sprintf("Majority-class (no-information) rate: %.2f%%\n",
            100 * max(table(y_test)) / length(y_test)))

# ---- FIGURE 7: RF error convergence -------------------------------------
save_fig("fig07_rf_error.png", {
  matplot(rf_model$err.rate, type = "l", lty = 1, lwd = 2,
          col = c("black", "#D55E00", "#0072B2"),
          xlab = "Number of Trees", ylab = "Error Rate",
          main = "Figure 7. RF Error Convergence")
  legend("topright", legend = colnames(rf_model$err.rate),
         col = c("black", "#D55E00", "#0072B2"), lty = 1, lwd = 2, cex = 0.8)
})

# ---- FIGURE 8 (1x2): MDA + MDG variable importance (top 20) -------------
imp <- as.data.frame(importance(rf_model))
imp$Feature <- rownames(imp)
write.csv(imp[order(-imp$MeanDecreaseGini), ],
          "figures/rf_feature_importance.csv", row.names = FALSE)

top_mda <- imp[order(-imp$MeanDecreaseAccuracy), ][1:20, ]
top_mdg <- imp[order(-imp$MeanDecreaseGini), ][1:20, ]

p_mda <- ggplot(top_mda, aes(x = reorder(short(Feature), MeanDecreaseAccuracy),
                             y = MeanDecreaseAccuracy)) +
  geom_col(fill = "#D55E00", width = 0.6) + coord_flip() +
  theme_minimal(base_size = 10) +
  labs(title = "A. Mean Decrease Accuracy (top 20)", x = "", y = "MDA")

p_mdg <- ggplot(top_mdg, aes(x = reorder(short(Feature), MeanDecreaseGini),
                             y = MeanDecreaseGini)) +
  geom_col(fill = "#0072B2", width = 0.6) + coord_flip() +
  theme_minimal(base_size = 10) +
  labs(title = "B. Mean Decrease Gini (top 20)", x = "", y = "MDG")

save_fig("fig08_importance.png", {
  grid.arrange(p_mda, p_mdg, ncol = 2,
               top = "Figure 8. Baseline RF Variable Importance")
}, width = 3600, height = 2000)

# ---- FIGURE 9 (1x2): confusion heatmap + precision/recall/F1 ------------
cm_df <- as.data.frame(conf_mat)
p_cm <- ggplot(cm_df, aes(Actual, Predicted, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), size = 6, fontface = "bold") +
  scale_fill_gradient(low = "white", high = "#0072B2") +
  theme_minimal(base_size = 12) +
  labs(title = "A. Confusion Matrix (Baseline RF)")

prec <- diag(conf_mat) / rowSums(conf_mat)
rec  <- diag(conf_mat) / colSums(conf_mat)
f1   <- 2 * prec * rec / (prec + rec)
metrics_df <- data.frame(
  Class  = rep(rownames(conf_mat), 3),
  Metric = rep(c("Precision", "Recall", "F1"), each = nrow(conf_mat)),
  Value  = c(prec, rec, f1)
)
p_met <- ggplot(metrics_df, aes(Class, Value, fill = Metric)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("#0072B2", "#E69F00", "#009E73")) +
  ylim(0, 1) + theme_minimal(base_size = 12) +
  labs(title = "B. Per-Class Performance", y = "Score")

save_fig("fig09_confusion_metrics.png", {
  grid.arrange(p_cm, p_met, ncol = 2,
               top = "Figure 9. Baseline RF Classification Performance")
}, width = 3000, height = 1700)

# ---- FIGURE 10: ROC curve (binary target) -------------------------------
save_fig("fig10_roc.png", {
  plot.roc(roc_base, col = "#0072B2", lwd = 2, legacy.axes = TRUE,
           main = "Figure 10. ROC Curve for ER Status (Baseline RF)")
  legend("bottomright", bty = "n", lwd = 2, col = "#0072B2",
         legend = sprintf("ER-positive vs negative (AUC = %.3f)", auc_base))
})

# ---- 5. HYBRID MODEL: SOM-augmented Random Forest -----------------------
# The SOM used for the hybrid features is fitted on the TRAINING set only;
# test patients are projected onto it (no information leakage).
# SOM unit is encoded by its lattice coordinates (randomForest cannot handle
# a 144-level factor), plus the super-cluster label.
cat("\n=== Hybrid SOM + Random Forest ===\n")

make_som_feats <- function(X_tr, X_te, k = K_SOM) {
  mu  <- colMeans(X_tr)
  sdv <- apply(X_tr, 2, sd); sdv[!is.finite(sdv) | sdv == 0] <- 1
  Z_tr <- clip(scale(X_tr, center = mu, scale = sdv))
  Z_te <- clip(scale(X_te, center = mu, scale = sdv))
  sm <- som(Z_tr, grid = somgrid(xdim = SOM_X, ydim = SOM_Y, topo = "hexagonal"),
            rlen = SOM_RLEN, alpha = c(0.05, 0.01))
  cl   <- cutree(hclust(dist(getCodes(sm))), k = k)
  u_tr <- sm$unit.classif
  u_te <- kohonen::map(sm, newdata = Z_te)$unit.classif
  feats <- function(u) data.frame(
    SOM_X       = sm$grid$pts[u, 1],
    SOM_Y       = sm$grid$pts[u, 2],
    SOM_Cluster = factor(cl[u], levels = 1:k)
  )
  list(train = feats(u_tr), test = feats(u_te))
}
som_vars <- c("SOM_X", "SOM_Y", "SOM_Cluster")

set.seed(SEED)
sf <- make_som_feats(X[train_idx, ], X[-train_idx, ])
train_hybrid <- cbind(x_train, sf$train)
test_hybrid  <- cbind(x_test,  sf$test)

set.seed(SEED)
rf_hybrid <- randomForest(x = train_hybrid, y = y_train, ntree = 300, importance = TRUE)
pred_hybrid <- predict(rf_hybrid, newdata = test_hybrid)
conf_hybrid <- table(Predicted = pred_hybrid, Actual = y_test)
acc_hybrid   <- sum(diag(conf_hybrid)) / sum(conf_hybrid)
kappa_hybrid <- cohen_kappa(conf_hybrid)
bacc_hybrid  <- bal_acc(conf_hybrid)
prob_hybrid  <- predict(rf_hybrid, newdata = test_hybrid, type = "prob")
roc_hybrid   <- roc(response = y_test, predictor = prob_hybrid[, "Positive"],
                    levels = c("Negative", "Positive"), direction = "<", quiet = TRUE)
auc_hybrid   <- as.numeric(auc(roc_hybrid))

cat("\n--- Hybrid Confusion Matrix ---\n"); print(conf_hybrid)
cat(sprintf("Hybrid Accuracy: %.2f%% | Balanced Acc: %.3f | Kappa: %.3f | AUC: %.3f\n",
            acc_hybrid * 100, bacc_hybrid, kappa_hybrid, auc_hybrid))

results_tab <- data.frame(
  Model       = c("Baseline RF", "Hybrid SOM+RF"),
  Accuracy    = c(accuracy, acc_hybrid),
  BalancedAcc = c(bacc_base, bacc_hybrid),
  Kappa       = c(kappa_base, kappa_hybrid),
  AUC         = c(auc_base, auc_hybrid)
)
cat("\n=== Model Comparison ===\n"); print(results_tab)
cat("DeLong test (baseline vs hybrid AUC):\n")
print(roc.test(roc_base, roc_hybrid, method = "delong"))
write.csv(results_tab, "figures/model_comparison.csv", row.names = FALSE)

# ---- FIGURE 11 (1x2): baseline vs hybrid + hybrid importance ------------
comp_df <- data.frame(
  Model  = rep(c("Baseline RF", "Hybrid SOM+RF"), 4),
  Metric = rep(c("Accuracy", "Balanced Acc", "Kappa", "AUC"), each = 2),
  Value  = c(accuracy, acc_hybrid, bacc_base, bacc_hybrid,
             kappa_base, kappa_hybrid, auc_base, auc_hybrid)
)
p_comp <- ggplot(comp_df, aes(Metric, Value, fill = Model)) +
  geom_col(position = "dodge", width = 0.6) +
  scale_fill_manual(values = c("#999999", "#0072B2")) +
  ylim(0, 1) + theme_minimal(base_size = 12) +
  labs(title = "A. Baseline vs Hybrid Performance", y = "Score")

imp_h <- as.data.frame(importance(rf_hybrid))
imp_h$Feature <- rownames(imp_h)
imp_h$Type <- ifelse(imp_h$Feature %in% som_vars, "SOM-derived", "Radiomic")
top_h <- imp_h[order(-imp_h$MeanDecreaseGini), ][1:20, ]
top_h <- unique(rbind(top_h, imp_h[imp_h$Feature %in% som_vars, ]))  # always show SOM vars

p_imp_h <- ggplot(top_h, aes(x = reorder(short(Feature), MeanDecreaseGini),
                             y = MeanDecreaseGini, fill = Type)) +
  geom_col(width = 0.6) + coord_flip() +
  scale_fill_manual(values = c("Radiomic" = "#0072B2", "SOM-derived" = "#D55E00")) +
  theme_minimal(base_size = 10) +
  labs(title = "B. Hybrid RF Variable Importance", x = "", y = "MDG")

save_fig("fig11_hybrid_comparison.png", {
  grid.arrange(p_comp, p_imp_h, ncol = 2,
               top = "Figure 11. Added Value of SOM Topological Features")
}, width = 3600, height = 2000)

# ---- FIGURE 12 (1x3): observed vs predicted ER on the SOM lattice -------
# Honest predictions for every patient: out-of-bag for training patients,
# hold-out predictions for test patients.
pred_all <- y
pred_all[train_idx]  <- rf_hybrid$predicted
pred_all[-train_idx] <- pred_hybrid

pred_pos_frac <- sapply(1:n_nodes, function(i) {
  idx <- which(unit_id == i)
  if (length(idx) == 0) NA else mean(pred_all[idx] == "Positive")
})
node_agree <- sapply(1:n_nodes, function(i) {
  idx <- which(unit_id == i)
  if (length(idx) == 0) NA else mean(pred_all[idx] == y[idx])
})

save_fig("fig12_som_validation.png", {
  par(mfrow = c(1, 3))
  plot(som_model, type = "property", property = pos_frac,
       main = "A. Observed ER-Positive Fraction")
  plot(som_model, type = "property", property = pred_pos_frac,
       main = "B. Hybrid RF Predicted ER-Positive")
  plot(som_model, type = "property", property = node_agree,
       main = "C. Node-wise Agreement (OOB + test)")
  par(mfrow = c(1, 1))
}, width = 3900, height = 1500)

# ---- FIGURE 13 (2x2): manual partial dependence plots -------------------
pdp_one <- function(model, data, var, class = "Positive", grid_n = 40) {
  grid <- seq(quantile(data[[var]], 0.02, na.rm = TRUE),
              quantile(data[[var]], 0.98, na.rm = TRUE),
              length.out = grid_n)
  pd <- sapply(grid, function(g) {
    tmp <- data
    tmp[[var]] <- g
    mean(predict(model, newdata = tmp, type = "prob")[, class])
  })
  data.frame(x = grid, y = pd)
}

num_imp <- imp_h[!imp_h$Feature %in% som_vars, ]
top4 <- num_imp$Feature[order(-num_imp$MeanDecreaseGini)][1:4]
print(top4)

pdp_list <- lapply(top4, function(v) {
  pd <- pdp_one(rf_hybrid, train_hybrid, v)
  ggplot(pd, aes(x, y)) +
    geom_line(color = "#0072B2", linewidth = 1.2) +
    geom_rug(data = train_hybrid, aes(x = .data[[v]]),
             inherit.aes = FALSE, alpha = 0.15) +
    coord_cartesian(xlim = range(pd$x)) +
    theme_minimal(base_size = 10) +
    labs(title = short(v, 40), x = "Feature value", y = "Mean P(ER-positive)")
})

save_fig("fig13_pdp.png", {
  grid.arrange(grobs = pdp_list, ncol = 2,
               top = "Figure 13. Partial Dependence on P(ER-positive) - Hybrid SOM+RF")
})

# ---- FIGURE 14: stability over 20 stratified random splits --------------
# SOM features are re-derived inside every split (train-only fit).
set.seed(123)
stab <- data.frame(Rep = integer(), Model = character(),
                   Acc = numeric(), AUC = numeric())
for (r in 1:20) {
  idx  <- strat_split(y)
  rf_b <- randomForest(x = Xdf[idx, ], y = y[idx], ntree = 200)
  sfr  <- make_som_feats(X[idx, ], X[-idx, ])
  tr_h <- cbind(Xdf[idx, ],  sfr$train)
  te_h <- cbind(Xdf[-idx, ], sfr$test)
  rf_h <- randomForest(x = tr_h, y = y[idx], ntree = 200)
  
  get_auc <- function(m, nd) as.numeric(auc(roc(
    response = y[-idx], predictor = predict(m, nd, type = "prob")[, "Positive"],
    levels = c("Negative", "Positive"), direction = "<", quiet = TRUE)))
  
  stab <- rbind(stab,
                data.frame(Rep = r, Model = "Baseline RF",
                           Acc = mean(predict(rf_b, Xdf[-idx, ]) == y[-idx]),
                           AUC = get_auc(rf_b, Xdf[-idx, ])),
                data.frame(Rep = r, Model = "Hybrid SOM+RF",
                           Acc = mean(predict(rf_h, te_h) == y[-idx]),
                           AUC = get_auc(rf_h, te_h)))
  cat("split", r, "done\n")
}
write.csv(stab, "figures/stability_20_splits.csv", row.names = FALSE)

mk_box <- function(var, ttl) {
  ggplot(stab, aes(Model, .data[[var]], fill = Model)) +
    geom_boxplot(width = 0.5, outlier.shape = 21) +
    geom_jitter(width = 0.12, alpha = 0.4, size = 1.5) +
    scale_fill_manual(values = c("#999999", "#0072B2")) +
    theme_minimal(base_size = 12) + theme(legend.position = "none") +
    labs(title = ttl, y = var, x = "")
}
save_fig("fig14_stability.png", {
  grid.arrange(mk_box("Acc", "A. Test Accuracy"), mk_box("AUC", "B. Test AUC"),
               ncol = 2, top = "Figure 14. Stability over 20 Stratified Splits")
}, width = 3000, height = 1700)

cat("\nMean over 20 splits:\n")
print(aggregate(cbind(Acc, AUC) ~ Model, data = stab, FUN = mean))
cat("Paired Wilcoxon test on AUC (hybrid vs baseline):\n")
print(wilcox.test(stab$AUC[stab$Model == "Hybrid SOM+RF"],
                  stab$AUC[stab$Model == "Baseline RF"], paired = TRUE))

cat("\n=== Files written to ./figures/ ===\n")
print(list.files("figures"))
cat("\nPipeline complete. 14 figures + result tables saved.\n")
sessionInfo()
# =============================================================================
# 6. POST HOC ANALYSIS: class-imbalance handling
#    (a) balanced Random Forest: each tree is grown on an equal number of
#        ER-negative and ER-positive patients (stratified down-sampling)
#    (b) decision threshold tuned on OUT-OF-BAG training probabilities
#        (Youden's J); the test set is never used for tuning
# =============================================================================
cat("\n=== Post hoc: balanced RF + OOB-tuned threshold ===\n")

fit_bal <- function(x, y, ntree = 300) {
  n_min <- min(table(y))
  randomForest(x = x, y = y, ntree = ntree,
               strata = y, sampsize = c(n_min, n_min))
}

oob_threshold <- function(rf, y) {
  p    <- rf$votes[, "Positive"]            # out-of-bag vote fractions
  grid <- seq(0.05, 0.95, by = 0.01)
  J <- sapply(grid, function(t) {
    pr <- p >= t
    mean(pr[y == "Positive"]) + mean(!pr[y == "Negative"]) - 1
  })
  grid[which.max(J)]
}

eval_model <- function(p, ytrue, thr) {
  pred <- factor(ifelse(p >= thr, "Positive", "Negative"),
                 levels = c("Negative", "Positive"))
  cm   <- table(Predicted = pred, Actual = ytrue)
  sens <- cm["Positive", "Positive"] / sum(cm[, "Positive"])
  spec <- cm["Negative", "Negative"] / sum(cm[, "Negative"])
  r <- roc(response = ytrue, predictor = p,
           levels = c("Negative", "Positive"), direction = "<", quiet = TRUE)
  list(cm = cm, roc = r,
       row = data.frame(Threshold = thr,
                        Accuracy = sum(diag(cm)) / sum(cm),
                        Sensitivity = sens, Specificity = spec,
                        BalancedAcc = (sens + spec) / 2,
                        Kappa = cohen_kappa(cm),
                        AUC = as.numeric(auc(r))))
}

# ---- Single hold-out split (same split as the main analysis) -------------
set.seed(SEED)
rf_bal   <- fit_bal(x_train, y_train)
set.seed(SEED)
rf_bal_h <- fit_bal(train_hybrid, y_train)

thr_b <- oob_threshold(rf_bal,   y_train)
thr_h <- oob_threshold(rf_bal_h, y_train)
cat("OOB-tuned thresholds | balanced RF:", thr_b, "| balanced hybrid:", thr_h, "\n")

p_bal   <- predict(rf_bal,   x_test,      type = "prob")[, "Positive"]
p_bal_h <- predict(rf_bal_h, test_hybrid, type = "prob")[, "Positive"]

e1 <- eval_model(prob_test[, "Positive"], y_test, 0.5)    # standard RF (main analysis)
e2 <- eval_model(p_bal,   y_test, 0.5)                    # balanced RF, default threshold
e3 <- eval_model(p_bal,   y_test, thr_b)                  # balanced RF, tuned threshold
e4 <- eval_model(p_bal_h, y_test, thr_h)                  # balanced hybrid, tuned threshold

mod_names <- c("Standard RF", "Balanced RF (0.5)",
               "Balanced RF (tuned)", "Balanced Hybrid (tuned)")
bal_tab <- cbind(Model = mod_names, rbind(e1$row, e2$row, e3$row, e4$row))
cat("\n--- Hold-out results ---\n"); print(bal_tab)
cat("\nConfusion matrix: balanced RF (tuned)\n");     print(e3$cm)
cat("\nConfusion matrix: balanced hybrid (tuned)\n"); print(e4$cm)
cat("\nDeLong: standard RF vs balanced RF\n")
print(roc.test(e1$roc, e3$roc, method = "delong"))
cat("\nDeLong: balanced RF vs balanced hybrid\n")
print(roc.test(e3$roc, e4$roc, method = "delong"))
write.csv(bal_tab, "figures/balanced_model_comparison.csv", row.names = FALSE)

# ---- 20 stratified splits ------------------------------------------------
set.seed(321)
stab2 <- data.frame()
for (r in 1:20) {
  idx <- strat_split(y)
  ytr <- y[idx]; yte <- y[-idx]
  
  rf_s <- randomForest(x = Xdf[idx, ], y = ytr, ntree = 200)
  rf_b <- fit_bal(Xdf[idx, ], ytr, ntree = 200)
  sfr  <- make_som_feats(X[idx, ], X[-idx, ])
  tr_h <- cbind(Xdf[idx, ],  sfr$train)
  te_h <- cbind(Xdf[-idx, ], sfr$test)
  rf_h <- fit_bal(tr_h, ytr, ntree = 200)
  
  rs <- eval_model(predict(rf_s, Xdf[-idx, ], type = "prob")[, "Positive"], yte, 0.5)$row
  rb <- eval_model(predict(rf_b, Xdf[-idx, ], type = "prob")[, "Positive"], yte,
                   oob_threshold(rf_b, ytr))$row
  rh <- eval_model(predict(rf_h, te_h, type = "prob")[, "Positive"], yte,
                   oob_threshold(rf_h, ytr))$row
  
  stab2 <- rbind(stab2,
                 cbind(Rep = r, Model = "Standard RF",     rs),
                 cbind(Rep = r, Model = "Balanced RF",     rb),
                 cbind(Rep = r, Model = "Balanced Hybrid", rh))
  cat("split", r, "done\n")
}
stab2$Model <- factor(stab2$Model,
                      levels = c("Standard RF", "Balanced RF", "Balanced Hybrid"))
write.csv(stab2, "figures/balanced_stability_20_splits.csv", row.names = FALSE)

cat("\n--- Mean over 20 splits ---\n")
print(aggregate(cbind(Accuracy, Sensitivity, Specificity, BalancedAcc, AUC) ~ Model,
                data = stab2, FUN = mean))
cat("\n--- SD over 20 splits ---\n")
print(aggregate(cbind(Accuracy, Sensitivity, Specificity, BalancedAcc, AUC) ~ Model,
                data = stab2, FUN = sd))

getv <- function(m, v) stab2[stab2$Model == m, v]
cat("\nWilcoxon (paired) balanced accuracy: Balanced RF vs Standard RF\n")
print(wilcox.test(getv("Balanced RF", "BalancedAcc"),
                  getv("Standard RF", "BalancedAcc"), paired = TRUE))
cat("\nWilcoxon (paired) AUC: Balanced RF vs Standard RF\n")
print(wilcox.test(getv("Balanced RF", "AUC"),
                  getv("Standard RF", "AUC"), paired = TRUE))
cat("\nWilcoxon (paired) balanced accuracy: Balanced Hybrid vs Balanced RF\n")
print(wilcox.test(getv("Balanced Hybrid", "BalancedAcc"),
                  getv("Balanced RF", "BalancedAcc"), paired = TRUE))
cat("\nWilcoxon (paired) AUC: Balanced Hybrid vs Balanced RF\n")
print(wilcox.test(getv("Balanced Hybrid", "AUC"),
                  getv("Balanced RF", "AUC"), paired = TRUE))

# ---- FIGURE 15 (1x3): effect of imbalance handling -----------------------
long_a <- do.call(rbind, lapply(c("Sensitivity", "Specificity", "BalancedAcc", "AUC"),
                                function(m) data.frame(Model = bal_tab$Model, Metric = m, Value = bal_tab[[m]])))
long_a$Model <- factor(long_a$Model, levels = mod_names)

p15a <- ggplot(long_a, aes(Metric, Value, fill = Model)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("#999999", "#56B4E9", "#0072B2", "#D55E00")) +
  ylim(0, 1) + theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.text.x = element_text(angle = 20, hjust = 1)) +
  guides(fill = guide_legend(nrow = 2)) +
  labs(title = "A. Hold-out Performance", x = "", y = "Score")

mk_box2 <- function(var, ttl, ylab) {
  ggplot(stab2, aes(Model, .data[[var]], fill = Model)) +
    geom_boxplot(width = 0.5, outlier.shape = 21) +
    geom_jitter(width = 0.12, alpha = 0.4, size = 1.5) +
    scale_fill_manual(values = c("#999999", "#0072B2", "#D55E00")) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 20, hjust = 1)) +
    labs(title = ttl, x = "", y = ylab)
}

save_fig("fig15_imbalance_handling.png", {
  grid.arrange(p15a,
               mk_box2("BalancedAcc", "B. Balanced Accuracy (20 splits)", "Balanced accuracy"),
               mk_box2("AUC", "C. AUC (20 splits)", "AUC"),
               ncol = 3,
               top = "Figure 15. Effect of Class-Imbalance Handling")
}, width = 4500, height = 1800)

cat("\nPost hoc analysis complete.\n")