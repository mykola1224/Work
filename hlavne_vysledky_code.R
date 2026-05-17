

rm(list = ls())
options(stringsAsFactors = FALSE)


data_path <- "C:/Users/ubcma/Desktop/data.RData"
out_dir <- "C:/Users/ubcma/Desktop/output_glm_grouped"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

LARGE_CLAIM_SHARE <- 0.035
SMALL_CLAIM_PROB  <- 1 - LARGE_CLAIM_SHARE   # 0.965

# =========================
# 1) BALÍKY
# =========================
required_packages <- c("dplyr", "MASS", "splines", "ggplot2")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(dplyr)
library(MASS)
library(splines)
library(ggplot2)

# =========================
# 2) NAČÍTANIE DÁT
# =========================
load(data_path)

train_freq_raw <- freMTPL2freq
test_freq_raw  <- freMTPLfreq

train_sev_raw  <- freMTPL2sev
test_sev_raw   <- freMTPLsev

# =========================
# 3) POMOCNÉ FUNKCIE
# =========================
poisson_deviance <- function(y, mu) {
  eps <- 1e-15
  mu <- pmax(mu, eps)
  term <- ifelse(y == 0, 0, y * log(y / mu))
  2 * sum(term - (y - mu))
}

gamma_deviance <- function(y, mu) {
  eps <- 1e-15
  y <- pmax(y, eps)
  mu <- pmax(mu, eps)
  2 * sum((y - mu) / mu - log(y / mu))
}

build_term <- function(var_name, data, df = 4, force_linear = c("BonusMalus")) {
  if (!(var_name %in% names(data))) return(NULL)
  
  x <- data[[var_name]]
  if (!is.numeric(x)) return(var_name)
  
  x <- x[is.finite(x) & !is.na(x)]
  n_unique <- length(unique(x))
  
  if (n_unique < 2) return(NULL)
  if (var_name %in% force_linear) return(var_name)
  if (n_unique < (df + 2)) return(var_name)
  
  paste0("splines::ns(", var_name, ", df = ", df, ")")
}

safe_formula_with_offset <- function(response, offset_term, rhs_terms) {
  rhs_terms <- rhs_terms[!is.na(rhs_terms) & nzchar(rhs_terms)]
  
  if (length(rhs_terms) == 0) {
    as.formula(paste(response, "~", offset_term))
  } else {
    as.formula(paste(response, "~", offset_term, "+", paste(rhs_terms, collapse = " + ")))
  }
}

safe_formula <- function(response, rhs_terms) {
  rhs_terms <- rhs_terms[!is.na(rhs_terms) & nzchar(rhs_terms)]
  
  if (length(rhs_terms) == 0) {
    as.formula(paste(response, "~ 1"))
  } else {
    as.formula(paste(response, "~", paste(rhs_terms, collapse = " + ")))
  }
}

align_factor_levels <- function(train_df, test_df, factor_vars) {
  for (v in factor_vars) {
    train_df[[v]] <- as.factor(train_df[[v]])
    test_df[[v]]  <- factor(test_df[[v]], levels = levels(train_df[[v]]))
  }
  list(train = train_df, test = test_df)
}

make_quantile_breaks <- function(x, n_groups = 20) {
  x <- as.numeric(x)
  x <- x[is.finite(x) & !is.na(x)]
  
  if (length(x) == 0) {
    return(c(-Inf, Inf))
  }
  
  probs <- seq(0, 1, length.out = n_groups + 1)
  br <- quantile(x, probs = probs, na.rm = TRUE, type = 7)
  br <- sort(unique(as.numeric(br)))
  
  if (length(br) < 2) {
    return(c(-Inf, Inf))
  }
  
  br[1] <- -Inf
  br[length(br)] <- Inf
  br
}

apply_quantile_groups <- function(x, breaks, prefix) {
  if (is.null(x)) stop(paste("Premenná", prefix, "chýba."))
  
  x <- as.numeric(x)
  if (length(x) == 0) stop(paste("Premenná", prefix, "má dĺžku 0."))
  
  if (length(breaks) < 2) {
    return(factor(rep(paste0(prefix, "_G1"), length(x))))
  }
  
  g <- cut(x, breaks = breaks, include.lowest = TRUE, right = TRUE)
  
  if (length(g) != length(x)) stop(paste("Zoskupenie zlyhalo pre", prefix))
  
  nlev <- nlevels(g)
  if (nlev == 0) stop(paste("Nevytvorili sa žiadne skupiny pre", prefix))
  
  levels(g) <- paste0(prefix, "_G", seq_len(nlev))
  g
}

safe_existing_vars <- function(train_df, test_df, vars) {
  vars[vars %in% names(train_df) & vars %in% names(test_df)]
}

eval_frequency <- function(df) {
  actual_claims <- df$ClaimNb
  pred_claims   <- df$pred_claims
  actual_freq   <- df$ClaimNb / df$Exposure
  pred_freq     <- df$pred_claims / df$Exposure
  
  data.frame(
    n = nrow(df),
    total_exposure = sum(df$Exposure, na.rm = TRUE),
    total_actual_claims = sum(actual_claims, na.rm = TRUE),
    total_predicted_claims = sum(pred_claims, na.rm = TRUE),
    actual_mean_frequency = sum(actual_claims, na.rm = TRUE) / sum(df$Exposure, na.rm = TRUE),
    predicted_mean_frequency = sum(pred_claims, na.rm = TRUE) / sum(df$Exposure, na.rm = TRUE),
    MAE_claims = mean(abs(actual_claims - pred_claims), na.rm = TRUE),
    RMSE_claims = sqrt(mean((actual_claims - pred_claims)^2, na.rm = TRUE)),
    MAE_frequency = mean(abs(actual_freq - pred_freq), na.rm = TRUE),
    RMSE_frequency = sqrt(mean((actual_freq - pred_freq)^2, na.rm = TRUE)),
    Poisson_deviance = poisson_deviance(actual_claims, pred_claims)
  )
}

eval_severity <- function(df) {
  data.frame(
    n = nrow(df),
    actual_mean_amount = mean(df$ClaimAmount, na.rm = TRUE),
    predicted_mean_amount = mean(df$pred_amount, na.rm = TRUE),
    MAE = mean(abs(df$ClaimAmount - df$pred_amount), na.rm = TRUE),
    RMSE = sqrt(mean((df$ClaimAmount - df$pred_amount)^2, na.rm = TRUE)),
    Gamma_deviance = gamma_deviance(df$ClaimAmount, df$pred_amount)
  )
}

# =========================
# 4) PRÍPRAVA DÁT
# =========================
prepare_freq_data <- function(df) {
  df <- as.data.frame(df)
  
  factor_vars_all  <- c("Area", "VehBrand", "VehGas", "Region")
  numeric_vars_all <- c("VehPower", "VehAge", "DrivAge", "BonusMalus", "Density")
  
  factor_vars  <- intersect(names(df), factor_vars_all)
  numeric_vars <- intersect(names(df), numeric_vars_all)
  
  for (v in factor_vars)  df[[v]] <- as.factor(df[[v]])
  for (v in numeric_vars) df[[v]] <- as.numeric(df[[v]])
  
  df <- df %>%
    filter(
      !is.na(Exposure), is.finite(Exposure), Exposure > 0,
      !is.na(ClaimNb), is.finite(ClaimNb)
    )
  
  df$claim_freq <- df$ClaimNb / df$Exposure
  
  list(
    data = df,
    factor_vars = factor_vars,
    numeric_vars = numeric_vars
  )
}

prepare_sev_data <- function(df) {
  df <- as.data.frame(df)
  
  factor_vars_all  <- c("Area", "VehBrand", "VehGas", "Region")
  numeric_vars_all <- c("VehPower", "VehAge", "DrivAge", "BonusMalus", "Density")
  
  factor_vars  <- intersect(names(df), factor_vars_all)
  numeric_vars <- intersect(names(df), numeric_vars_all)
  
  for (v in factor_vars)  df[[v]] <- as.factor(df[[v]])
  for (v in numeric_vars) df[[v]] <- as.numeric(df[[v]])
  
  df <- df %>%
    filter(!is.na(ClaimAmount), is.finite(ClaimAmount), ClaimAmount > 0)
  
  list(
    data = df,
    factor_vars = factor_vars,
    numeric_vars = numeric_vars
  )
}

# =========================
# 5) HRANICA ŠKODY 96,5 / 3,5
# =========================
train_sev_valid_for_threshold <- as.data.frame(train_sev_raw) %>%
  filter(!is.na(ClaimAmount), is.finite(ClaimAmount), ClaimAmount > 0)

large_claim_threshold <- as.numeric(
  quantile(
    train_sev_valid_for_threshold$ClaimAmount,
    probs = SMALL_CLAIM_PROB,
    na.rm = TRUE
  )
)

threshold_table <- data.frame(
  typ_hranice = "96_5_percent_trening",
  hodnota_hranice = large_claim_threshold
)

write.csv(
  threshold_table,
  file.path(out_dir, "claim_threshold_96_5.csv"),
  row.names = FALSE
)

cat("\nHranica škody (96,5 % / 3,5 %):", round(large_claim_threshold, 2), "\n")

# =========================
# 6) DÁTA PRE FREKVENCIU
# =========================
train_freq_obj <- prepare_freq_data(train_freq_raw)
test_freq_obj  <- prepare_freq_data(test_freq_raw)

train_freq <- train_freq_obj$data
test_freq  <- test_freq_obj$data

group_vars <- c("DrivAge", "VehAge", "VehPower")
group_vars <- Reduce(intersect, list(group_vars, names(train_freq), names(test_freq)))

freq_group_breaks <- list()

for (v in group_vars) {
  freq_group_breaks[[v]] <- make_quantile_breaks(train_freq[[v]], n_groups = 20)
  
  train_freq[[paste0(v, "_grp")]] <- apply_quantile_groups(
    train_freq[[v]],
    freq_group_breaks[[v]],
    v
  )
  
  test_freq[[paste0(v, "_grp")]] <- apply_quantile_groups(
    test_freq[[v]],
    freq_group_breaks[[v]],
    v
  )
}

grouped_factor_vars <- paste0(group_vars, "_grp")

common_factor_vars  <- intersect(train_freq_obj$factor_vars, test_freq_obj$factor_vars)
common_numeric_vars <- intersect(train_freq_obj$numeric_vars, test_freq_obj$numeric_vars)

common_numeric_vars <- setdiff(common_numeric_vars, group_vars)
common_factor_vars  <- unique(c(common_factor_vars, grouped_factor_vars))

keep_vars <- unique(c("Exposure", "ClaimNb", common_factor_vars, common_numeric_vars))
keep_vars <- safe_existing_vars(train_freq, test_freq, keep_vars)

train_freq <- train_freq[, keep_vars, drop = FALSE]
test_freq  <- test_freq[, keep_vars, drop = FALSE]

common_factor_vars  <- intersect(common_factor_vars, keep_vars)
common_numeric_vars <- intersect(common_numeric_vars, keep_vars)

aligned_freq <- align_factor_levels(train_freq, test_freq, common_factor_vars)
train_freq <- aligned_freq$train
test_freq  <- aligned_freq$test

for (v in common_factor_vars) {
  test_freq <- test_freq %>% filter(!is.na(.data[[v]]))
}

# =========================
# 7) MODELY FREKVENCIE
# =========================
rhs_linear <- c(common_numeric_vars, common_factor_vars)
formula_linear <- safe_formula_with_offset("ClaimNb", "offset(log(Exposure))", rhs_linear)

model_pois_linear <- glm(
  formula_linear,
  family = poisson(link = "log"),
  data = train_freq
)

rhs_flex <- c(
  sapply(common_numeric_vars, build_term, data = train_freq, USE.NAMES = FALSE),
  common_factor_vars
)
rhs_flex <- rhs_flex[!is.na(rhs_flex) & nzchar(rhs_flex)]

formula_flex <- safe_formula_with_offset("ClaimNb", "offset(log(Exposure))", rhs_flex)

model_pois_flex <- glm(
  formula_flex,
  family = poisson(link = "log"),
  data = train_freq
)

best_poisson <- model_pois_linear
best_poisson_name <- "Poisson_linear"

if (AIC(model_pois_flex) < AIC(best_poisson)) {
  best_poisson <- model_pois_flex
  best_poisson_name <- "Poisson_flexible"
}

dispersion_ratio <- sum(residuals(best_poisson, type = "pearson")^2) / best_poisson$df.residual

model_nb <- NULL
if (dispersion_ratio > 1.5) {
  model_nb <- tryCatch(
    MASS::glm.nb(formula(best_poisson), data = train_freq),
    error = function(e) NULL
  )
  
  if (!is.null(model_nb) && AIC(model_nb) < AIC(best_poisson)) {
    best_freq_model <- model_nb
    best_freq_model_name <- "Negative_Binomial"
  } else {
    best_freq_model <- best_poisson
    best_freq_model_name <- best_poisson_name
  }
} else {
  best_freq_model <- best_poisson
  best_freq_model_name <- best_poisson_name
}

train_freq$pred_claims <- predict(best_freq_model, newdata = train_freq, type = "response")
train_freq$pred_freq   <- train_freq$pred_claims / train_freq$Exposure

test_freq$pred_claims <- predict(best_freq_model, newdata = test_freq, type = "response")
test_freq$pred_freq   <- test_freq$pred_claims / test_freq$Exposure

train_freq_metrics <- eval_frequency(train_freq)
test_freq_metrics  <- eval_frequency(test_freq)

cat("\n================ FREKVENCIA ================\n")
cat("Zvolený model frekvencie:", best_freq_model_name, "\n")
cat("Pomer disperzie:", round(dispersion_ratio, 4), "\n")
cat("\nMETRIKY FREKVENCIE - TRÉNOVACIA MNOŽINA\n")
print(train_freq_metrics)
cat("\nMETRIKY FREKVENCIE - TESTOVACIA MNOŽINA\n")
print(test_freq_metrics)

aic_table_freq <- data.frame(
  model = c("Poisson_linear", "Poisson_flexible"),
  AIC = c(AIC(model_pois_linear), AIC(model_pois_flex))
)

if (!is.null(model_nb)) {
  aic_table_freq <- rbind(
    aic_table_freq,
    data.frame(model = "Negative_Binomial", AIC = AIC(model_nb))
  )
}

write.csv(
  aic_table_freq,
  file.path(out_dir, "frequency_model_AIC.csv"),
  row.names = FALSE
)

write.csv(
  summary(best_freq_model)$coefficients,
  file.path(out_dir, "best_frequency_model_coefficients.csv")
)

capture.output(
  summary(best_freq_model),
  file = file.path(out_dir, "best_frequency_model_summary.txt")
)

freq_metrics_table <- bind_rows(
  cbind(mnozina = "trening", model = "frekvencia", train_freq_metrics),
  cbind(mnozina = "test", model = "frekvencia", test_freq_metrics)
)

write.csv(
  freq_metrics_table,
  file.path(out_dir, "frequency_metrics_train_test.csv"),
  row.names = FALSE
)

# -------------------------
# GRAFY FREKVENCIE
# -------------------------
train_freq_deciles <- train_freq %>%
  mutate(pred_decile = dplyr::ntile(pred_freq, 10)) %>%
  group_by(pred_decile) %>%
  summarise(
    expozicia = sum(Exposure, na.rm = TRUE),
    skutocny_pocet_skod = sum(ClaimNb, na.rm = TRUE),
    predikovany_pocet_skod = sum(pred_claims, na.rm = TRUE),
    skutocna_frekvencia = skutocny_pocet_skod / expozicia,
    predikovana_frekvencia = predikovany_pocet_skod / expozicia,
    .groups = "drop"
  )

test_freq_deciles <- test_freq %>%
  mutate(pred_decile = dplyr::ntile(pred_freq, 10)) %>%
  group_by(pred_decile) %>%
  summarise(
    expozicia = sum(Exposure, na.rm = TRUE),
    skutocny_pocet_skod = sum(ClaimNb, na.rm = TRUE),
    predikovany_pocet_skod = sum(pred_claims, na.rm = TRUE),
    skutocna_frekvencia = skutocny_pocet_skod / expozicia,
    predikovana_frekvencia = predikovany_pocet_skod / expozicia,
    .groups = "drop"
  )

write.csv(
  train_freq_deciles,
  file.path(out_dir, "train_frequency_deciles.csv"),
  row.names = FALSE
)

write.csv(
  test_freq_deciles,
  file.path(out_dir, "test_frequency_deciles.csv"),
  row.names = FALSE
)

p_freq_train <- ggplot(train_freq_deciles, aes(x = pred_decile)) +
  geom_line(aes(y = skutocna_frekvencia, color = "Skutočná"), linewidth = 1) +
  geom_point(aes(y = skutocna_frekvencia, color = "Skutočná"), size = 2) +
  geom_line(aes(y = predikovana_frekvencia, color = "Predikovaná"), linewidth = 1) +
  geom_point(aes(y = predikovana_frekvencia, color = "Predikovaná"), size = 2) +
  labs(
    title = "Trénovacia množina: skutočná a predikovaná frekvencia podľa decilov",
    x = "Decil predikovanej frekvencie",
    y = "Frekvencia škôd",
    color = ""
  ) +
  theme_minimal()

ggsave(
  file.path(out_dir, "train_frequency_actual_vs_predicted_deciles.png"),
  plot = p_freq_train,
  width = 10,
  height = 6,
  dpi = 150
)

p_freq_test <- ggplot(test_freq_deciles, aes(x = pred_decile)) +
  geom_line(aes(y = skutocna_frekvencia, color = "Skutočná"), linewidth = 1) +
  geom_point(aes(y = skutocna_frekvencia, color = "Skutočná"), size = 2) +
  geom_line(aes(y = predikovana_frekvencia, color = "Predikovaná"), linewidth = 1) +
  geom_point(aes(y = predikovana_frekvencia, color = "Predikovaná"), size = 2) +
  labs(
    title = "Testovacia množina: skutočná a predikovaná frekvencia podľa decilov",
    x = "Decil predikovanej frekvencie",
    y = "Frekvencia škôd",
    color = ""
  ) +
  theme_minimal()

ggsave(
  file.path(out_dir, "test_frequency_actual_vs_predicted_deciles.png"),
  plot = p_freq_test,
  width = 10,
  height = 6,
  dpi = 150
)

# =========================
# 8) ROZDELENIE ŠKÔD 96,5 / 3,5
# =========================
train_sev_small_raw <- as.data.frame(train_sev_raw) %>%
  filter(
    !is.na(ClaimAmount), is.finite(ClaimAmount), ClaimAmount > 0,
    ClaimAmount <= large_claim_threshold
  )

train_sev_large_raw <- as.data.frame(train_sev_raw) %>%
  filter(
    !is.na(ClaimAmount), is.finite(ClaimAmount), ClaimAmount > 0,
    ClaimAmount > large_claim_threshold
  )

test_sev_small_raw <- as.data.frame(test_sev_raw) %>%
  filter(
    !is.na(ClaimAmount), is.finite(ClaimAmount), ClaimAmount > 0,
    ClaimAmount <= large_claim_threshold
  )

test_sev_large_raw <- as.data.frame(test_sev_raw) %>%
  filter(
    !is.na(ClaimAmount), is.finite(ClaimAmount), ClaimAmount > 0,
    ClaimAmount > large_claim_threshold
  )

split_summary <- data.frame(
  mnozina = c("trening_mensie_96_5", "trening_velke_3_5", "test_mensie_96_5", "test_velke_3_5"),
  n = c(
    nrow(train_sev_small_raw), nrow(train_sev_large_raw),
    nrow(test_sev_small_raw), nrow(test_sev_large_raw)
  ),
  celkova_suma = c(
    sum(train_sev_small_raw$ClaimAmount),
    sum(train_sev_large_raw$ClaimAmount),
    sum(test_sev_small_raw$ClaimAmount),
    sum(test_sev_large_raw$ClaimAmount)
  ),
  priemerna_skoda = c(
    mean(train_sev_small_raw$ClaimAmount),
    mean(train_sev_large_raw$ClaimAmount),
    mean(test_sev_small_raw$ClaimAmount),
    mean(test_sev_large_raw$ClaimAmount)
  )
)

write.csv(
  split_summary,
  file.path(out_dir, "claim_split_summary_96_5_3_5.csv"),
  row.names = FALSE
)

write.csv(
  train_sev_large_raw,
  file.path(out_dir, "train_large_claims_3_5.csv"),
  row.names = FALSE
)

write.csv(
  test_sev_large_raw,
  file.path(out_dir, "test_large_claims_3_5.csv"),
  row.names = FALSE
)

# =========================
# 9) MODEL ZÁVAŽNOSTI LEN PRE MENŠIE ŠKODY 96,5 %
# =========================
train_sev_obj <- prepare_sev_data(train_sev_small_raw)
test_sev_obj  <- prepare_sev_data(test_sev_small_raw)

train_sev <- train_sev_obj$data
test_sev  <- test_sev_obj$data

group_vars_sev <- c("DrivAge", "VehAge", "VehPower")
group_vars_sev <- Reduce(intersect, list(group_vars_sev, names(train_sev), names(test_sev)))

sev_group_breaks <- list()

for (v in group_vars_sev) {
  sev_group_breaks[[v]] <- make_quantile_breaks(train_sev[[v]], n_groups = 20)
  
  train_sev[[paste0(v, "_grp")]] <- apply_quantile_groups(
    train_sev[[v]],
    sev_group_breaks[[v]],
    v
  )
  
  test_sev[[paste0(v, "_grp")]] <- apply_quantile_groups(
    test_sev[[v]],
    sev_group_breaks[[v]],
    v
  )
}

grouped_factor_vars_sev <- paste0(group_vars_sev, "_grp")

common_factor_vars_sev  <- intersect(train_sev_obj$factor_vars, test_sev_obj$factor_vars)
common_numeric_vars_sev <- intersect(train_sev_obj$numeric_vars, test_sev_obj$numeric_vars)

common_numeric_vars_sev <- setdiff(common_numeric_vars_sev, group_vars_sev)
common_factor_vars_sev  <- unique(c(common_factor_vars_sev, grouped_factor_vars_sev))

keep_vars_sev <- unique(c("ClaimAmount", common_factor_vars_sev, common_numeric_vars_sev))
keep_vars_sev <- safe_existing_vars(train_sev, test_sev, keep_vars_sev)

train_sev <- train_sev[, keep_vars_sev, drop = FALSE]
test_sev  <- test_sev[, keep_vars_sev, drop = FALSE]

common_factor_vars_sev  <- intersect(common_factor_vars_sev, keep_vars_sev)
common_numeric_vars_sev <- intersect(common_numeric_vars_sev, keep_vars_sev)

aligned_sev <- align_factor_levels(train_sev, test_sev, common_factor_vars_sev)
train_sev <- aligned_sev$train
test_sev  <- aligned_sev$test

for (v in common_factor_vars_sev) {
  test_sev <- test_sev %>% filter(!is.na(.data[[v]]))
}

rhs_sev <- c(
  sapply(common_numeric_vars_sev, build_term, data = train_sev, USE.NAMES = FALSE),
  common_factor_vars_sev
)
rhs_sev <- rhs_sev[!is.na(rhs_sev) & nzchar(rhs_sev)]

formula_sev <- safe_formula("ClaimAmount", rhs_sev)

model_sev <- glm(
  formula_sev,
  family = Gamma(link = "log"),
  data = train_sev
)

train_sev$pred_amount <- predict(model_sev, newdata = train_sev, type = "response")
test_sev$pred_amount  <- predict(model_sev, newdata = test_sev, type = "response")

train_sev_metrics <- eval_severity(train_sev)
test_sev_metrics  <- eval_severity(test_sev)

cat("\n================ ZÁVAŽNOSŤ (96,5 % MENŠÍCH ŠKÔD) ================\n")
cat("\nMETRIKY ZÁVAŽNOSTI - TRÉNOVACIA MNOŽINA\n")
print(train_sev_metrics)
cat("\nMETRIKY ZÁVAŽNOSTI - TESTOVACIA MNOŽINA\n")
print(test_sev_metrics)

write.csv(
  summary(model_sev)$coefficients,
  file.path(out_dir, "severity_model_coefficients_small_96_5.csv")
)

capture.output(
  summary(model_sev),
  file = file.path(out_dir, "severity_model_summary_small_96_5.txt")
)

sev_metrics_table <- bind_rows(
  cbind(mnozina = "trening_mensie_96_5", model = "zavaznost", train_sev_metrics),
  cbind(mnozina = "test_mensie_96_5", model = "zavaznost", test_sev_metrics)
)

write.csv(
  sev_metrics_table,
  file.path(out_dir, "severity_metrics_train_test_small_96_5.csv"),
  row.names = FALSE
)

# -------------------------
# GRAFY ZÁVAŽNOSTI
# -------------------------
train_sev_deciles <- train_sev %>%
  mutate(pred_decile = dplyr::ntile(pred_amount, 10)) %>%
  group_by(pred_decile) %>%
  summarise(
    skutocna_priemerna_skoda = mean(ClaimAmount, na.rm = TRUE),
    predikovana_priemerna_skoda = mean(pred_amount, na.rm = TRUE),
    .groups = "drop"
  )

test_sev_deciles <- test_sev %>%
  mutate(pred_decile = dplyr::ntile(pred_amount, 10)) %>%
  group_by(pred_decile) %>%
  summarise(
    skutocna_priemerna_skoda = mean(ClaimAmount, na.rm = TRUE),
    predikovana_priemerna_skoda = mean(pred_amount, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  train_sev_deciles,
  file.path(out_dir, "train_severity_deciles.csv"),
  row.names = FALSE
)

write.csv(
  test_sev_deciles,
  file.path(out_dir, "test_severity_deciles.csv"),
  row.names = FALSE
)

p_sev_train <- ggplot(train_sev_deciles, aes(x = pred_decile)) +
  geom_line(aes(y = skutocna_priemerna_skoda, color = "Skutočná"), linewidth = 1) +
  geom_point(aes(y = skutocna_priemerna_skoda, color = "Skutočná"), size = 2) +
  geom_line(aes(y = predikovana_priemerna_skoda, color = "Predikovaná"), linewidth = 1) +
  geom_point(aes(y = predikovana_priemerna_skoda, color = "Predikovaná"), size = 2) +
  labs(
    title = "Trénovacia množina: skutočná a predikovaná závažnosť podľa decilov",
    x = "Decil predikovanej závažnosti",
    y = "Priemerná výška škody",
    color = ""
  ) +
  theme_minimal()

ggsave(
  file.path(out_dir, "train_severity_actual_vs_predicted_deciles.png"),
  plot = p_sev_train,
  width = 10,
  height = 6,
  dpi = 150
)

p_sev_test <- ggplot(test_sev_deciles, aes(x = pred_decile)) +
  geom_line(aes(y = skutocna_priemerna_skoda, color = "Skutočná"), linewidth = 1) +
  geom_point(aes(y = skutocna_priemerna_skoda, color = "Skutočná"), size = 2) +
  geom_line(aes(y = predikovana_priemerna_skoda, color = "Predikovaná"), linewidth = 1) +
  geom_point(aes(y = predikovana_priemerna_skoda, color = "Predikovaná"), size = 2) +
  labs(
    title = "Testovacia množina: skutočná a predikovaná závažnosť podľa decilov",
    x = "Decil predikovanej závažnosti",
    y = "Priemerná výška škody",
    color = ""
  ) +
  theme_minimal()

ggsave(
  file.path(out_dir, "test_severity_actual_vs_predicted_deciles.png"),
  plot = p_sev_test,
  width = 10,
  height = 6,
  dpi = 150
)

p_sev_scatter <- ggplot(test_sev, aes(x = ClaimAmount, y = pred_amount)) +
  geom_point(alpha = 0.3) +
  labs(
    title = "Testovacia množina: skutočná a predikovaná výška škody",
    x = "Skutočná výška škody",
    y = "Predikovaná výška škody"
  ) +
  theme_minimal()

ggsave(
  file.path(out_dir, "test_severity_actual_vs_predicted_scatter.png"),
  plot = p_sev_scatter,
  width = 10,
  height = 6,
  dpi = 150
)

# =========================
# 10) 5-NÁSOBNÁ KRÍŽOVÁ VALIDÁCIA PRE FREKVENCIU
# =========================
set.seed(123)

cv_freq_obj <- prepare_freq_data(train_freq_raw)
cv_freq <- cv_freq_obj$data

if ("IDpol" %in% names(cv_freq)) {
  unit_ids <- unique(cv_freq$IDpol)
  fold_id_assignment <- sample(rep(1:5, length.out = length(unit_ids)))
  fold_map <- data.frame(IDpol = unit_ids, fold = fold_id_assignment)
  cv_freq <- cv_freq %>% left_join(fold_map, by = "IDpol")
} else if ("PolicyID" %in% names(cv_freq)) {
  unit_ids <- unique(cv_freq$PolicyID)
  fold_id_assignment <- sample(rep(1:5, length.out = length(unit_ids)))
  fold_map <- data.frame(PolicyID = unit_ids, fold = fold_id_assignment)
  cv_freq <- cv_freq %>% left_join(fold_map, by = "PolicyID")
} else {
  cv_freq$fold <- sample(rep(1:5, length.out = nrow(cv_freq)))
}

cv_results <- list()
cv_predictions <- list()

for (k in 1:5) {
  train_fold <- cv_freq %>% filter(fold != k)
  valid_fold <- cv_freq %>% filter(fold == k)
  
  fold_group_vars <- c("DrivAge", "VehAge", "VehPower")
  fold_group_vars <- Reduce(intersect, list(fold_group_vars, names(train_fold), names(valid_fold)))
  
  fold_group_breaks <- list()
  
  for (v in fold_group_vars) {
    fold_group_breaks[[v]] <- make_quantile_breaks(train_fold[[v]], n_groups = 20)
    
    train_fold[[paste0(v, "_grp")]] <- apply_quantile_groups(
      train_fold[[v]],
      fold_group_breaks[[v]],
      v
    )
    
    valid_fold[[paste0(v, "_grp")]] <- apply_quantile_groups(
      valid_fold[[v]],
      fold_group_breaks[[v]],
      v
    )
  }
  
  factor_vars_all  <- c("Area", "VehBrand", "VehGas", "Region", paste0(fold_group_vars, "_grp"))
  numeric_vars_all <- c("BonusMalus", "Density")
  
  factor_vars_train  <- intersect(names(train_fold), factor_vars_all)
  numeric_vars_train <- intersect(names(train_fold), numeric_vars_all)
  
  factor_vars_valid  <- intersect(names(valid_fold), factor_vars_all)
  numeric_vars_valid <- intersect(names(valid_fold), numeric_vars_all)
  
  common_factor_vars_cv  <- intersect(factor_vars_train, factor_vars_valid)
  common_numeric_vars_cv <- intersect(numeric_vars_train, numeric_vars_valid)
  
  keep_vars_cv <- unique(c("Exposure", "ClaimNb", common_factor_vars_cv, common_numeric_vars_cv, "fold"))
  keep_vars_cv <- keep_vars_cv[keep_vars_cv %in% names(train_fold) & keep_vars_cv %in% names(valid_fold)]
  
  if ("IDpol" %in% names(train_fold) && "IDpol" %in% names(valid_fold)) {
    keep_vars_cv <- unique(c(keep_vars_cv, "IDpol"))
  }
  if ("PolicyID" %in% names(train_fold) && "PolicyID" %in% names(valid_fold)) {
    keep_vars_cv <- unique(c(keep_vars_cv, "PolicyID"))
  }
  
  train_fold <- train_fold[, keep_vars_cv, drop = FALSE]
  valid_fold <- valid_fold[, keep_vars_cv, drop = FALSE]
  
  aligned <- align_factor_levels(train_fold, valid_fold, common_factor_vars_cv)
  train_fold <- aligned$train
  valid_fold <- aligned$test
  
  for (v in common_factor_vars_cv) {
    valid_fold <- valid_fold %>% filter(!is.na(.data[[v]]))
  }
  
  rhs_linear_cv <- c(common_numeric_vars_cv, common_factor_vars_cv)
  formula_linear_cv <- safe_formula_with_offset("ClaimNb", "offset(log(Exposure))", rhs_linear_cv)
  
  model_pois_linear_cv <- glm(
    formula_linear_cv,
    family = poisson(link = "log"),
    data = train_fold
  )
  
  rhs_flex_cv <- c(
    sapply(common_numeric_vars_cv, build_term, data = train_fold, USE.NAMES = FALSE),
    common_factor_vars_cv
  )
  rhs_flex_cv <- rhs_flex_cv[!is.na(rhs_flex_cv) & nzchar(rhs_flex_cv)]
  
  formula_flex_cv <- safe_formula_with_offset("ClaimNb", "offset(log(Exposure))", rhs_flex_cv)
  
  model_pois_flex_cv <- glm(
    formula_flex_cv,
    family = poisson(link = "log"),
    data = train_fold
  )
  
  best_poisson_cv <- model_pois_linear_cv
  best_poisson_cv_name <- "Poisson_linear"
  
  if (AIC(model_pois_flex_cv) < AIC(best_poisson_cv)) {
    best_poisson_cv <- model_pois_flex_cv
    best_poisson_cv_name <- "Poisson_flexible"
  }
  
  dispersion_ratio_cv <- sum(residuals(best_poisson_cv, type = "pearson")^2) / best_poisson_cv$df.residual
  
  if (dispersion_ratio_cv > 1.5) {
    model_nb_cv <- tryCatch(
      MASS::glm.nb(formula(best_poisson_cv), data = train_fold),
      error = function(e) NULL
    )
    
    if (!is.null(model_nb_cv) && AIC(model_nb_cv) < AIC(best_poisson_cv)) {
      best_model_cv <- model_nb_cv
      best_model_cv_name <- "Negative_Binomial"
    } else {
      best_model_cv <- best_poisson_cv
      best_model_cv_name <- best_poisson_cv_name
    }
  } else {
    best_model_cv <- best_poisson_cv
    best_model_cv_name <- best_poisson_cv_name
  }
  
  valid_fold$pred_claims <- predict(best_model_cv, newdata = valid_fold, type = "response")
  valid_fold$pred_freq   <- valid_fold$pred_claims / valid_fold$Exposure
  
  fold_metrics <- eval_frequency(valid_fold)
  fold_metrics$fold <- k
  fold_metrics$model_used <- best_model_cv_name
  fold_metrics$dispersion_ratio <- dispersion_ratio_cv
  
  cv_results[[k]] <- fold_metrics
  cv_predictions[[k]] <- valid_fold
}

cv_results_table <- dplyr::bind_rows(cv_results)
cv_predictions_table <- dplyr::bind_rows(cv_predictions)

cv_summary_table <- data.frame(
  metrika = c(
    "actual_mean_frequency",
    "predicted_mean_frequency",
    "MAE_claims",
    "RMSE_claims",
    "MAE_frequency",
    "RMSE_frequency",
    "Poisson_deviance",
    "dispersion_ratio"
  ),
  priemer = c(
    mean(cv_results_table$actual_mean_frequency, na.rm = TRUE),
    mean(cv_results_table$predicted_mean_frequency, na.rm = TRUE),
    mean(cv_results_table$MAE_claims, na.rm = TRUE),
    mean(cv_results_table$RMSE_claims, na.rm = TRUE),
    mean(cv_results_table$MAE_frequency, na.rm = TRUE),
    mean(cv_results_table$RMSE_frequency, na.rm = TRUE),
    mean(cv_results_table$Poisson_deviance, na.rm = TRUE),
    mean(cv_results_table$dispersion_ratio, na.rm = TRUE)
  ),
  smerodajna_odchylka = c(
    sd(cv_results_table$actual_mean_frequency, na.rm = TRUE),
    sd(cv_results_table$predicted_mean_frequency, na.rm = TRUE),
    sd(cv_results_table$MAE_claims, na.rm = TRUE),
    sd(cv_results_table$RMSE_claims, na.rm = TRUE),
    sd(cv_results_table$MAE_frequency, na.rm = TRUE),
    sd(cv_results_table$RMSE_frequency, na.rm = TRUE),
    sd(cv_results_table$Poisson_deviance, na.rm = TRUE),
    sd(cv_results_table$dispersion_ratio, na.rm = TRUE)
  )
)

write.csv(
  cv_results_table,
  file.path(out_dir, "frequency_cv_5fold_results.csv"),
  row.names = FALSE
)

write.csv(
  cv_summary_table,
  file.path(out_dir, "frequency_cv_5fold_summary.csv"),
  row.names = FALSE
)

write.csv(
  cv_predictions_table,
  file.path(out_dir, "frequency_cv_5fold_predictions.csv"),
  row.names = FALSE
)

# =========================
# 11) ULOŽENIE INFORMÁCIÍ O RELÁCII
# =========================
capture.output(
  sessionInfo(),
  file = file.path(out_dir, "sessionInfo.txt")
)

# =========================
# 12) ZÁVEREČNÝ VÝPIS
# =========================
cat("\n============================\n")
cat("Výpočet bol úspešne dokončený.\n")
cat("Model frekvencie:", best_freq_model_name, "\n")
cat("Pomer disperzie:", round(dispersion_ratio, 4), "\n")
cat("Hranica škody (96,5 % / 3,5 %):", round(large_claim_threshold, 2), "\n")
cat("Výsledky sú uložené v priečinku:\n")
cat(normalizePath(out_dir, winslash = "/", mustWork = FALSE), "\n")
cat("============================\n")