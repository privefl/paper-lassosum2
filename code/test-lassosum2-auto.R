library(bigsnpr)
ukb <- snp_attach("data/UKBB_imp_HM3.rds")

load("data/ind_gwas_val_test.RData")

(NCORES <- parallelly::availableCores() - 1L)

chr <- 22
ind.chr <- which(ukb$map$chromosome == chr)
corr0 <- readRDS(paste0("data/corr/chr", chr, ".rds"))
corr <- runonce::save_run(as_SFBM(corr0, "tmp-data/corr_chr22"),
                          "tmp-data/corr_chr22.rds")
stopifnot(length(ind.chr) == ncol(corr))
dim(corr)  # 16410 x 16410

# big_copy(ukb$genotypes, ind.col = ind.chr,
#          backingfile = "tmp-data/dosage_chr22")$save()
G <- big_attach("tmp-data/dosage_chr22.rds")
y <- snp_simuPheno(G, h2 = 0.2, M = 2000, ncores = NCORES)$pheno

ind.gwas <- sample(ind.gwas, 20e3)

# GWAS to get sumstats
gwas <- big_univLinReg(G, y[ind.gwas], ind.train = ind.gwas, ncores = NCORES)
library(dplyr)
df_beta <- gwas %>%
  transmute(beta = estim, beta_se = std.err, n_eff = length(ind.gwas))

# lassosum2
beta_lassosum <- snp_lassosum2(corr, df_beta, ncores = NCORES)
(params <- attr(beta_lassosum, "grid_param"))

# validation
ind <- which(rowSums(beta_lassosum != 0) > 0)
pred_lassosum <- big_prodMat(G, beta_lassosum[ind, ], ind.col = ind, ncores = NCORES)
params$score <- big_univLinReg(as_FBM(pred_lassosum[ind.val, ]), y[ind.val])$score

# pseudo-validation
scale <- with(df_beta, sqrt(n_eff * beta_se^2 + beta^2))
beta_hat <- df_beta$beta / scale

fdr <- fdrtool::fdrtool(beta_hat, statistic = "correlation", plot = FALSE)
beta_hat_shrunk <- round(beta_hat * (1 - fdr$lfdr), 16)

params$auto_score <- apply(beta_lassosum, 2, function(beta) {
  cat(".")
  beta <- beta / scale
  bRb <- crossprod(beta, bigsparser::sp_prodVec(corr, beta))
  crossprod(beta, beta_hat_shrunk) / sqrt(bRb)
})

library(ggplot2)
qplot(auto_score, score, color = s, data = params) +
  theme_bw(15) +
  scale_color_viridis_c() +
  labs(x = "Score from pseudo-validation", y = "Score from validation")


pval <- predict(gwas, log10 = FALSE)
fdr2 <- fdrtool::fdrtool(pval, statistic = "pvalue", plot = FALSE)
beta_hat_shrunk2 <- beta_hat * (1 - fdr2$lfdr)

params$auto_score2 <- apply(beta_lassosum, 2, function(beta) {
  cat(".")
  beta <- beta / scale
  bRb <- crossprod(beta, bigsparser::sp_prodVec(corr, beta))
  crossprod(beta, beta_hat_shrunk2) / sqrt(bRb)
})

plot_grid(
  qplot(auto_score, score, color = s, data = params) +
    theme_bw(15) +
    scale_color_viridis_c() +
    labs(x = "Score from pseudo-validation (using correlations)",
         y = "Score from validation"),
  qplot(auto_score2, score, color = s, data = params) +
    theme_bw(15) +
    scale_color_viridis_c() +
    labs(x = "Score from pseudo-validation (using p-values)",
         y = "Score from validation"),
  scale = 0.95, labels = c("A", "B"), label_size = 16, ncol = 1
)
# ggsave("figures/pseudoval.pdf", width = 8, height = 9)
